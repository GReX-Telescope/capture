// -*- c++ -*-mp1*2**cand1_widths[i]
#include <iostream>
#include <algorithm>
using std::cout;
using std::cerr;
using std::endl;
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <unistd.h>
#include <netdb.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <time.h>

#include "dada_cuda.h"
#include "dada_client.h"
#include "dada_def.h"
#include "dada_hdu.h"
#include "multilog.h"
#include "ipcio.h"
#include "ipcbuf.h"
#include "dada_affinity.h"
#include "ascii_header.h"

#include <thrust/fill.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/functional.h>
#include <thrust/transform.h>

#include "fitsio.h"

#define NCHAN 2048
#define NANT 1

// these are in bytes
#define BLOCKSIZE 409600000
#define GULP 51200000
#define NSAMPS_GULP 12500
#define NGULPS 8

// binning for single-point flagging
#define bf 16
#define bt 25

void dsaX_dbgpu_cleanup (dada_hdu_t * hdu_in, dada_hdu_t * hdu_out, multilog_t * log);
int dada_bind_thread_to_core (int core);

// functor to do the scaling
__device__ float *s1, *s2, *s3;

struct da_functor
{

  int nant;
  da_functor(int _nant) : nant(_nant) {}

  __device__
  int operator()(const int x, const int y) const {

    int i = (int)(y % (2048)); // CHANGE FOR NANT
    //printf("x = %f",x);
    /*std::cout << "s1 = " << s1[i] <<std::endl;
    std::cout << "s2 = " << s2[i] <<std::endl;
    std::cout << "s3 = " << s3[i] <<std::endl;
    std::cout << "x s1/s2 + s3 = " << x*s1[i]/s2[i]+s3[i] <<std::endl;
    std::cout << "x s1/s2 + s3 r = " << __float2int_rn(x*s1[i]/s2[i]+s3[i]) <<std::endl;*/
    return __float2int_rn(x*s1[i]/s2[i]+s3[i]);
    
  }
};

    
int dada_bind_thread_to_core(int core)
{

  cpu_set_t set;
  pid_t tpid;

  CPU_ZERO(&set);
    CPU_SET(core, &set);
      tpid = syscall(SYS_gettid);

  if (sched_setaffinity(tpid, sizeof(cpu_set_t), &set) < 0) {
      fprintf(stderr, "failed to set cpu affinity: %s", strerror(errno));
          return -1;
	    }

  CPU_ZERO(&set);
    if ( sched_getaffinity(tpid, sizeof(cpu_set_t), &set) < 0 ) {
        fprintf(stderr, "failed to get cpu affinity: %s", strerror(errno));
	    return -1;
	      }

  return 0;
}

void usage()
{
  fprintf (stdout,
	   "dsaX_spectrometer_reorder [options]\n"
	   " -c core   bind process to CPU core\n"
	   " -d        dump 1024 spectra from each input block\n"
	   " -h        print usage\n");
}

int main (int argc, char *argv[]) {

  cudaSetDevice(1);
  
  /* DADA Header plus Data Unit */
  dada_hdu_t* hdu_in = 0;
  dada_hdu_t* hdu_out = 0;

  /* DADA Logger */
  multilog_t* log = 0;

  int core = -1;

  // input data block HDU key
  key_t in_key = 0x0000dada;

  // output data block HDU key
  key_t out_key = 0x0000eada;

  // online update of BF and BT
  int BF = bf;
  int BT = bt;
  int tf, tt;
  FILE *fFT;
  
  int arg = 0;
  int dump = 0;
  int snapchoice = -1;

  while ((arg=getopt(argc,argv,"c:s:dh")) != -1)
    {
      switch (arg)
	{
	case 'c':
	  if (optarg)
	    {
	      core = atoi(optarg);
	      break;
	    }
	  else
	    {
	      fprintf (stderr, "ERROR: -c flag requires argument\n");
	      return EXIT_FAILURE;
	    }

	case 's':
	  if (optarg)
	    {
	      snapchoice = atoi(optarg);
	      break;
	    }
	  else
	    {
	      fprintf (stderr, "ERROR: -s flag requires argument\n");
	      return EXIT_FAILURE;
	    }

	case 'd':
	  dump=1;
	case 'h':
	  usage();
	  return EXIT_SUCCESS;
	}
    }
  
  // DADA stuff

  log = multilog_open ("dsaX_spectrometer_reorder", 0);

  multilog_add (log, stderr);

  multilog (log, LOG_INFO, "dsaX_spectrometer_reorder: creating in hdu\n");

  // open connection to the in/read DB
  hdu_in  = dada_hdu_create (log);
  dada_hdu_set_key (hdu_in, in_key);
  if (dada_hdu_connect (hdu_in) < 0) {
    fprintf (stderr, "dsaX_spectrometer_reorder: could not connect to dada buffer\n");
    return EXIT_FAILURE;
  }
  if (dada_hdu_lock_read (hdu_in) < 0) {
    fprintf (stderr, "dsaX_spectrometer_reorder: could not lock to dada buffer\n");
    return EXIT_FAILURE;
  }

  // open connection to the out/write DB
  hdu_out = dada_hdu_create (log);
  dada_hdu_set_key (hdu_out, out_key);
  if (dada_hdu_connect (hdu_out) < 0)
    {
      dsaX_dbgpu_cleanup (hdu_in, hdu_out, log);
      fprintf (stderr, "dsaX_spectrometer_reorder: could not connect to eada buffer\n");
      return EXIT_FAILURE;
    }
  if (dada_hdu_lock_write(hdu_out) < 0)
    {
      dsaX_dbgpu_cleanup (hdu_in, hdu_out, log);
      fprintf (stderr, "dsaX_spectrometer_reorder: could not lock4 to eada buffer\n");
      return EXIT_FAILURE;
    }

  if (core >= 0)
    {
      fprintf(stderr, "binding to core %d\n", core);
      if (dada_bind_thread_to_core(core) < 0)
	fprintf(stderr, "dsaX_spectrometer_reorder: failed to bind to core %d\n", core);
    }

  bool observation_complete=0;

  // more DADA stuff
  
  uint64_t header_size = 0;

  // read the header from the input HDU
  char * header_in = ipcbuf_get_next_read (hdu_in->header_block, &header_size);
  if (!header_in)
    {
      multilog(log ,LOG_ERR, "main: could not read next header\n");
      dsaX_dbgpu_cleanup (hdu_in, hdu_out, log);
      return EXIT_FAILURE;
    }

  // now write the output DADA header
  char * header_out = ipcbuf_get_next_write (hdu_out->header_block);
  if (!header_out)
    {
      multilog(log, LOG_ERR, "could not get next header block [output]\n");
      dsaX_dbgpu_cleanup (hdu_in, hdu_out, log);
      return EXIT_FAILURE;
    }

  // read the number of stations from the header
  int nant;
  if (ascii_header_get (header_in, "NANT", "%d", &nant) != 1)
    {
      nant = NANT;
      multilog(log, LOG_WARNING, "Header with no NANT. Setting to %d\n", NANT);
    }

  fprintf(stderr, "dsaX_spectrometer_reorder: Have NANT %d\n", nant);
  fprintf(stderr, "dsaX_spectrometer_reorder: Have snapchoice %d\n", snapchoice);
  
  
  // copy the in header to the out header
  memcpy (header_out, header_in, header_size);

  // need to change some DADA parameters
  if (ascii_header_set (header_out, "NANT", "%d", 1) < 0)
    multilog(log, LOG_WARNING, "failed to set NANT 1 in header_out\n");

  // mark the input header as cleared
  if (ipcbuf_mark_cleared (hdu_in->header_block) < 0)
    {
      multilog (log, LOG_ERR, "could not mark header block cleared [input]\n");
      dsaX_dbgpu_cleanup (hdu_in, hdu_out, log);
      return EXIT_FAILURE;
    }

  // mark the output header buffer as filled
  if (ipcbuf_mark_filled (hdu_out->header_block, header_size) < 0)
    {
      multilog (log, LOG_ERR, "could not mark header block filled [output]\n");
      dsaX_dbgpu_cleanup (hdu_in, hdu_out, log);
      return EXIT_FAILURE;
    }

  //dada_cuda_dbregister (hdu_in);

  // to scale data
  thrust::host_vector<float> mult(2048*nant), sc(nant*2048), mlt(2048*nant);
  thrust::host_vector<float> mmult(2048*nant), msc(nant*2048), mmlt(2048*nant);
  thrust::device_vector<float> d_mult(2048*nant), d_sc(nant*2048), d_mlt(2048*nant);
  float *s_mult, *s_mlt, *s_sc;
  s_mult = thrust::raw_pointer_cast(&d_mult[0]);
  s_mlt = thrust::raw_pointer_cast(&d_mlt[0]);
  s_sc = thrust::raw_pointer_cast(&d_sc[0]);
  float bpscl = NSAMPS_GULP*64;

  for (int i=0;i<2048;i++) {
    for (int j=0;j<nant;j++) {
      sc[i+j*2048] = 0.;
      mult[i+j*2048] = bpscl;
      if (i<330) mult[i+j*2048]=0.;
//if (i>1354 && i<1367) mult[i+j*2048] = 0.;
//if (i>1982 && i<1986) mult[i+j*2048] = 0.;
//if (i>1965 && i<1974) mult[i+j*2048] = 0.;
//if (i>1920 && i<1935) mult[i+j*2048] = 0.;
/*if (i>1748 && i<1764) mult[i+j*2048] = 0.;
if (i>1704 && i<1724) mult[i+j*2048] = 0.; 
if (i>1673 && i<1691) mult[i+j*2048] = 0.;
if (i>1629 && i<1651) mult[i+j*2048] = 0.;
if (i>1359 && i<1593) mult[i+j*2048] = 0.;
if (i>1056 && i<1243) mult[i+j*2048] = 0.;
if (i>988 && i<1028) mult[i+j*2048] = 0.;
if (i>541 && i<543) mult[i+j*2048] = 0.;
    if (i>2033) mult[i+j*2048]=0.;*/
/*if (i>1354 && i<1367) mult[i+j*2048] = 0.;      
if (i>1528 && i<1547) mult[i+j*2048] = 0.;
if (i>1570 && i<1589) mult[i+j*2048] = 0.;
if (i>1207 && i<1220) mult[i+j*2048] = 0.;
if (i>515 && i<552) mult[i+j*2048] = 0.;
      if (i>1208 && i<1232) mult[i+j*2048] = 0.;
      if (i>1466 && i<1482) mult[i+j*2048] = 0.;
      if (i>1777 && i<1794) mult[i+j*2048] = 0.;*/
      //if (i>1888) mult[i+j*2048] = 0.;
      if (mult[i+j*2048]==0.) mlt[i+j*2048]=64.;
      else mlt[i+j*2048] = 0.;
    }
  }
   
  
  /*cudaMemcpyToSymbol(s1,&s_mult,sizeof(float *));
  cudaMemcpyToSymbol(s2,&s_sc,sizeof(float *));
  cudaMemcpyToSymbol(s3,&s_mlt,sizeof(float *));*/

  // file for logging flagged spectra
  fitsfile *fptr;
  char fitsnam[100];
  int status=0;
  int rownum = 1;
  time_t rawtime;
  struct tm *info;
  time(&rawtime);
  info = localtime(&rawtime);
  //double MJD = (double)(57754.+info->tm_yday+(info->tm_hour+8.)/24.+info->tm_min/(24.*60.)+info->tm_sec/(24.*60.*60.));
  //char MJD[100];
  long double MJD;
   if (ascii_header_get (header_in, "MJD_START", "%Lf", &MJD) != 1)
    {
      MJD = (long double)(57754.+info->tm_yday+(info->tm_hour+8.)/24.+info->tm_min/(24.*60.)+info->tm_sec/(24.*60.*60.));
      //MJD = (char*)(&MJD_double);
      multilog(log, LOG_WARNING, "Header with no MJD_START. Setting to %s\n", &MJD);
    }
float MJD_f = (float)MJD;
//printf("reorder: %0.30Lf\n",MJD);
 
  sprintf(fitsnam,"/home/user/data/slog_%0.3Lf.fits",MJD);
  char *ttype[] = {"Spectra","Perc_ts","Perc_samp"};
  char *tform[] = {"2048E", "E", "E"}; // EDIT FOR NANT
  char *tunit[] = {"\0", "\0", "\0"};
  char extname[] = "spec_log";
  fits_create_file(&fptr, fitsnam, &status);
  if (status) cerr << "create_file FITS error " << status << endl;
  fits_create_tbl(fptr, BINARY_TBL, 0, 3, ttype, tform, tunit, extname, &status);
  if (status) cerr << "create_tbl FITS error " << status << endl;
  fits_write_key(fptr, TFLOAT, "MJD", &MJD_f, "Start MJD", &status);
  float mytsamp = NSAMPS_GULP*6.5536e-5*2;
  fits_write_key(fptr, TFLOAT, "TSAMP", &mytsamp, "Sample time (s)", &status);

  if (status) cerr << "FITS error " << status << endl;
  else
    cout << "Opened FITS file " << fitsnam << endl;
  fits_close_file(fptr, &status);
  float out_bp[nant*2048], out_pts[1], out_psamp[1];

  // setu up vars
  uint64_t block_size = ipcbuf_get_bufsz ((ipcbuf_t *) hdu_in->data_block);
  uint64_t bytes_to_read;
  uint64_t block_id;
  char *   block;
  int bytes_to_write = BLOCKSIZE;
  uint64_t bytes_written=0;
  uint64_t written=0;
  int ibyte, idx, idxo, idxs, idxg;
  
  // allocate memory to output
  unsigned short * out_data;
  thrust::host_vector<int> h_gulpdata(2048*nant*NSAMPS_GULP);
  thrust::device_vector<int> d_gulpdata(2048*nant*NSAMPS_GULP);
  thrust::device_vector<int> d_idx(2048*nant*NSAMPS_GULP);
  thrust::sequence(d_idx.begin(),d_idx.end());
  out_data = (unsigned short *)malloc(sizeof(unsigned short)*BLOCKSIZE/2);
  char *h_indata;
  uint64_t  bytes_read = 0, clipped, clipped_ts;
  int bp[2048*nant], old_bp[2048*nant], ts[NSAMPS_GULP];
  int ts_sum, ds_sum;
  int thresh_ts = (int)((2048*nant*64*BT)+(4.*sqrt(2048.*nant*BT)*10.));
  float thresh_diff = 10*(0.32/sqrt(NSAMPS_GULP*1.));
  //unsigned short repval = (unsigned short)(av_bp);
  //if (snapchoice!=-1)
  //  repval = (unsigned short)(av_bp);
  int clipthresh = (int)(64.*nant*BT*BF+50000.*sqrt(1.*nant*BT*BF))*100000;
  int clipthresh_ss = (int)(64.*nant+45.*sqrt(1.*nant));
  if (snapchoice!=-1)
    clipthresh = 104000000*200*200000;
  unsigned short tmp;
  int intct = 0;
  int started_recording = 0;
  char cmd[200];
  uint64_t specnum = 0;
  int nints=0;
  
  multilog(log, LOG_INFO, "main: starting observation\n");
  
  while (!observation_complete) {

    if (nints > 4395) {

      rownum = 1;
      nints=0;
      time(&rawtime);
      info = localtime(&rawtime);
      MJD = (double)(57754.+info->tm_yday+(info->tm_hour+8.)/24.+info->tm_min/(24.*60.)+info->tm_sec/(24.*60.*60.));
      sprintf(fitsnam,"/home/user/data/slog_%.3Lf.fits",MJD);
      char *ttype[] = {"Spectra","Perc_ts","Perc_samp"};
      char *tform[] = {"2048E", "E", "E"}; // EDIT FOR NANT
      char *tunit[] = {"\0", "\0", "\0"};
      char extname[] = "spec_log";
      fits_create_file(&fptr, fitsnam, &status);
      if (status) cerr << "create_file FITS error " << status << endl;
      fits_create_tbl(fptr, BINARY_TBL, 0, 3, ttype, tform, tunit, extname, &status);
      if (status) cerr << "create_tbl FITS error " << status << endl;
      fits_write_key(fptr, TFLOAT, "MJD", &MJD_f, "Start MJD", &status);
      fits_write_key(fptr, TFLOAT, "TSAMP", &mytsamp, "Sample time (s)", &status);

      if (status) cerr << "FITS error " << status << endl;
      else
	cout << "Opened FITS file " << fitsnam << endl;
      fits_close_file(fptr, &status);

    }

    // deal with BF and BT
    if ((fFT=fopen("/home/user/runtime/fFT.dat","r"))) {
      fscanf(fFT,"%d %d\n",&tf,&tt);
      if (tf!=BF) {
	multilog(log, LOG_INFO, "main: new BF %d\n", tf);
	BF=tf;
      }
      if (tt!=BT) {
	multilog(log, LOG_INFO, "main: new BT %d\n", tt);
	BT=tt;
      }
      fclose(fFT);
    } 
    for (int gulp=0;gulp<NGULPS;gulp++) {
//cout << "gulp = " << gulp << std::endl;
      // read a DADA block
      h_indata = ipcio_open_block_read (hdu_in->data_block, &bytes_read, &block_id);
//cout << "h_indata size" << strlen(h_indata) << std::endl;
      //multilog(log, LOG_INFO, "main: read block which contains %lld bytes\n", bytes_read);
      
      /* Assuming mean of 64, each 8x20-sample point will have stddev of 10. 
	 Then, over NSAMPS samples, the stddev is 10/sqrt(NSAMPS*nant). Can 
	 flag every 10-sigma channel, and every 6-sigma burst */

      // deal with zero-ing stuff and setting up bandpasses.
      thrust::fill(h_gulpdata.begin(),h_gulpdata.end(),0);
      clipped = 0;
      clipped_ts = 0;
      for (int i=0;i<NSAMPS_GULP;i++) ts[i] = 0;
      for (int i=0;i<2048*nant;i++) {
	if (started_recording) old_bp[i]=bp[i];
	else old_bp[i] = 0;
	bp[i] = 0;
      }
    int sum_bp = 0;
    int av_bp_n = 0;
    int av_bp;
      // unpack data into h_gulpdata, and find current bp
      for (int k=0;k<NSAMPS_GULP;k++) {
	
	for (int snp=0;snp<nant;snp++) {
	  
	  for (int i=0;i<512;i++) {
	    for (int j=0;j<4;j++) {
	    
	      idx = k*nant*4096+snp*4096+i*8+j*2;
	      idxg = k*nant*2048+snp*2048+i*4+j;
	      idxs = i*4+j+snp*2048;
	      tmp=0;
	      tmp |= (unsigned short)(h_indata[idx]) << 8;
	      //cout << "tmp1 " << tmp << std::endl;
	      tmp |= (unsigned short)(h_indata[idx+1]);
	      //cout << "tmp2 " << tmp << std::endl;
	      //cout << "h_indata[idx] = " << (unsigned short)h_indata[idx] << std::endl;
              //cout << "h_indata[idx+1] = " << (unsigned short)h_indata[idx+1] << std::endl;
	      h_gulpdata[idxg] = (int)tmp;
	      //cout << "h_gulpdata[idxg] = " << (int)tmp << std::endl;
	      bp[idxs] += (int)tmp;
	      sum_bp += (int)tmp;
	      av_bp_n+=1;
//if(i==247 && k == 1) cout << int(tmp) << std::endl;
	    }
	    
	  }
	
	}
	
      }

av_bp = sum_bp/av_bp_n;
unsigned short repval = 64.;
int chans_clipped = 0;
int bp_clipped = 0;
int spec_clipped = 0;
//cout << "av_bp " << av_bp <<std::endl;
      // set up scaling by bandpass
      for (int i=0;i<2048*nant;i++) {
	mmult[i] = mult[i];
	mmlt[i] = mlt[i];
	msc[i] = bp[i]*1.;
	if (old_bp[i]==0) {
	  if (chans_clipped < 820) {
	  	mmult[i] = 0.;
	  	mmlt[i] = 64.;
	  	bp_clipped += 1;
	  	chans_clipped +=1;
	  }
	}
	else if ((((float)(bp[i]-old_bp[i]))/((float)(old_bp[i]))>thresh_diff) || (((float)(old_bp[i]-bp[i]))/((float)(old_bp[i]))>thresh_diff)) {
		if (chans_clipped < 820) {
			//mult[i] = 0;
			//mlt[i] = 64.;
	  		mmult[i] = 0.;
	  		mmlt[i] = 64.;
	  		chans_clipped += 1.;
			bp_clipped += 1.;
		}
	}
      }
      thrust::copy(mmult.begin(),mmult.end(),d_mult.begin());
      thrust::copy(mmlt.begin(),mmlt.end(),d_mlt.begin());
      thrust::copy(msc.begin(),msc.end(),d_sc.begin());
      s_mult = thrust::raw_pointer_cast(&d_mult[0]);
      s_mlt = thrust::raw_pointer_cast(&d_mlt[0]);
      s_sc = thrust::raw_pointer_cast(&d_sc[0]);
      cudaMemcpyToSymbol(s1,&s_mult,sizeof(float *));
      cudaMemcpyToSymbol(s2,&s_sc,sizeof(float *));
      cudaMemcpyToSymbol(s3,&s_mlt,sizeof(float *));
	//std::cout << "s1 = " << mmult[1000] << std::endl;
	//std::cout << "s2 = " << msc[1000] << std::endl;
	//std::cout << "s3 = " << mmlt[1000] << std::endl;
	//std::cout << "x = " << h_gulpdata[1000] << std::endl;
	//std::cout << "x*s1/s2+s3 = " << h_gulpdata[1000]*mmult[1000]/msc[1000]+mlt[1000] << std::endl;
      // do bandpass scaling of data
      thrust::copy(h_gulpdata.begin(),h_gulpdata.end(),d_gulpdata.begin());
      thrust::transform(d_gulpdata.begin(),d_gulpdata.end(),d_idx.begin(),d_gulpdata.begin(),da_functor(nant));
      thrust::copy(d_gulpdata.begin(),d_gulpdata.end(),h_gulpdata.begin());

      // reduce to out_data, and find ts
      for (int k=0;k<NSAMPS_GULP;k++) {
	for (int snp=0;snp<nant;snp++) {    
	  for (int i=0;i<2048;i++) {
	    idxo = (k+gulp*NSAMPS_GULP)*2048+i;
	    idxg = k*nant*2048+snp*2048+i;
	    if (snapchoice !=-1) {
	      if (snp==snapchoice) out_data[idxo] = (unsigned short)(h_gulpdata[idxg]);
	    }
	    else {
	      if (snp==0) out_data[idxo] = (unsigned short)(h_gulpdata[idxg]);
	      else out_data[idxo] = out_data[idxo] + (unsigned short)(h_gulpdata[idxg]);
	    }
	    ts[k] += (h_gulpdata[idxg]);
	  }
	}
      }
int max_clipped = 0;
      if (snapchoice==-1) {
	for (int k=0;k<NSAMPS_GULP;k++) {
	  int ind_spec_clipped = 0; 
	  for (int i=0;i<2048;i++) {
	   if (out_data[(k+gulp*NSAMPS_GULP)*2048+i]>clipthresh_ss && ind_spec_clipped+bp_clipped < 820) {
	      out_data[(k+gulp*NSAMPS_GULP)*2048+i]=repval;
	      clipped++;
	      chans_clipped++;
	      spec_clipped++;
	      ind_spec_clipped++;
	    }
	}
		if (ind_spec_clipped > max_clipped) {
			max_clipped = ind_spec_clipped;
		}
	  
	  
	}
      }
   std::cout << "Channels clipped: " << bp_clipped+max_clipped<< std::endl;
std::cout << "bp clipped: " << bp_clipped << std::endl;
std::cout << "Spec clipped: " << max_clipped << std::endl;

      if (snapchoice==-1) {
	for (int k=0;k<NSAMPS_GULP/BT;k++) {

	  // time-series flagging
	  ts_sum = 0;
	  /*for (int i=k*BT;i<(k+1)*BT;i++)
	    ts_sum += ts[i];
	  if (ts_sum>thresh_ts) {
	    //cout << ts_sum << " " << thresh_ts << endl;
	    clipped_ts+=BT;
	    for (int j=k*BT;j<(k+1)*BT;j++) {
	      for (int i=0;i<2048;i++) 
		out_data[(j+gulp*NSAMPS_GULP)*2048+i] = repval;
	    }
	  }*/
	  // single-point flagging
	  for (int i=0;i<2048/BF;i++) {
	    ds_sum = 0;
	    for (int j=k*BT;j<(k+1)*BT;j++) {
	      for (int l=i*BF;l<(i+1)*BF;l++)
		ds_sum += out_data[(j+gulp*NSAMPS_GULP)*2048+l];
	    }
	    if (ds_sum>clipthresh) {
	      for (int j=k*BT;j<(k+1)*BT;j++) {
		//for (int l=i*BF;l<(i+1)*BF;l++)
		  //out_data[(j+gulp*NSAMPS_GULP)*2048+l] = repval;
		//clipped+=BT*BF;
	      }
	    }
	  }
	 if(k==0) {
                for(int i=0;i<2048;i++) {
                        out_data[(gulp*NSAMPS_GULP)*2048+i] = repval + 0.3;
                }
          } 
	}
      }
      
    
      // do logging
      //multilog(log, LOG_INFO, "main: processed %d bytes, clipped percentage %.10f, samples %lld, with zero-DM clipping %lld of %d\n", GULP, (float)(100.*clipped/(GULP/2.)), clipped, clipped_ts, NSAMPS_GULP);
//cout << "clipped percentage " << (float)(100.*clipped/(NSAMPS_GULP*2048)) << std::endl;
      fits_open_table(&fptr, fitsnam, READWRITE, &status);
      for (int i=0;i<2048*nant;i++) {
	//if (mmult[i]!=0) out_bp[i] = 1.*bp[i];
	if (bp[i]>=0) out_bp[i] = 1.*bp[i];  
	else out_bp[i] = -1.*bp[i];
      }
      out_pts[0] = (float)(100.*clipped_ts/(NSAMPS_GULP));
      out_psamp[0] = (float)(100.*clipped/(NSAMPS_GULP*2048));
      fits_write_col(fptr, TFLOAT, 1, rownum, 1, nant*2048, out_bp, &status);
      fits_write_col(fptr, TFLOAT, 2, rownum, 1, 1, out_pts, &status);
      fits_write_col(fptr, TFLOAT, 3, rownum, 1, 1, out_psamp, &status);
      if (status) cerr << "FITS error in write " << status << endl;

      rownum += 1;
      fits_update_key(fptr, TINT, "NAXIS2", &rownum, "", &status);
      fits_close_file(fptr, &status);
      nints++;
      
      /*if (intct>0) flog = fopen("/mnt/nfs/data/spectrometer_log.dat","a");
      for (int i=0;i<2048*nant;i++) {
	if (mmult[i]!=0) fprintf(flog,"%d\n",(int)(old_bp[i]));
	else fprintf(flog,"%d\n",-(int)(old_bp[i]));
	//fprintf(flog,"%d\n",bp[i]-old_bp[i]);
      }
      fprintf(flog,"%g\n",(float)(100.*clipped/(GULP/2.)));
      fprintf(flog,"%g\n",(float)(100.*clipped_ts/(NSAMPS_GULP)));
      fclose(flog);*/
      intct++;

      // do the start
      started_recording = 1;

      // close block for reading
      ipcio_close_block_read (hdu_in->data_block, bytes_read);
//cout << "Done with gulp" << std::endl;
      if (bytes_read < block_size) {
	observation_complete = 1;
	multilog(log, LOG_INFO, "main: finished, with bytes_read %llu < expected %llu\n", bytes_read, block_size);
	break;
      }
      
    }
    // DO THE WRITING
    written = ipcio_write (hdu_out->data_block, (char *) out_data, bytes_to_write);
    
    if (written < bytes_to_write)
      {
	multilog(log, LOG_INFO, "main: failed to write all data to datablock [output]\n");
	dsaX_dbgpu_cleanup (hdu_in, hdu_out, log);
	return EXIT_FAILURE;
      }
    bytes_written += written;
    cout << "Wrote " << bytes_written << " bytes"  << endl;
    //specnum += 100000;
    //multilog(log, LOG_INFO, "main: write %llu bytes, %llu total, specnum %llu\n", written, bytes_written, specnum*16-1000000);
    
    //sprintf(cmd,"echo p %lu | nc -4u -w1 10.10.1.11 11223",specnum);
    //system(cmd);
    
  }
    
  dsaX_dbgpu_cleanup (hdu_in, hdu_out, log);
  free(out_data);
  
}

void dsaX_dbgpu_cleanup (dada_hdu_t * in, dada_hdu_t * out, multilog_t * log)
{

  //dada_cuda_dbunregister (in);
  
  if (dada_hdu_unlock_read (in) < 0)
    {
      multilog(log, LOG_ERR, "could not unlock read on hdu_in\n");
    }
  dada_hdu_destroy (in);

  if (dada_hdu_unlock_write (out) < 0)
    {
      multilog(log, LOG_ERR, "could not unlock write on hdu_out\n");
    }
  dada_hdu_destroy (out);
}
