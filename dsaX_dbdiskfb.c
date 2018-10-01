/* Code to read from a raw data buffer and write to disk */

#include <time.h>
#include <sys/socket.h>
#include <math.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/syscall.h>
#include <sys/mman.h>
#include <sched.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <netdb.h>
#include <cmath>

#include "sock.h"
#include "tmutil.h"
#include "dada_client.h"
#include "dada_def.h"
#include "dada_hdu.h"
#include "multilog.h"
#include "ipcio.h"
#include "ipcbuf.h"
#include "dada_affinity.h"
#include "ascii_header.h"
#include "header.h"
#include "sigproc.h"

FILE *input, *output;
int swapout;

void send_string(char *string) /* includefile */
{
  int len;
  len=strlen(string);
  fwrite(&len, sizeof(int), 1, output);
  fwrite(string, sizeof(char), len, output);
  /*fprintf(stderr,"%s\n",string);*/
}

void send_float(char *name,float floating_point) /* includefile */
{
  send_string(name);
  fwrite(&floating_point,sizeof(float),1,output);
  /*fprintf(stderr,"%f\n",floating_point);*/
}

void send_double (char *name, double double_precision) /* includefile */
{
  send_string(name);
  fwrite(&double_precision,sizeof(double),1,output);
  /*fprintf(stderr,"%f\n",double_precision);*/
}

void send_int(char *name, int integer) /* includefile */
{
  send_string(name);
  fwrite(&integer,sizeof(int),1,output);
  /*fprintf(stderr,"%d\n",integer);*/
}

void send_long(char *name, long integer) /* includefile */
{
  send_string(name);
  fwrite(&integer,sizeof(long),1,output);
  /*fprintf(stderr,"%ld\n",integer);*/
}



void dsaX_dbgpu_cleanup (dada_hdu_t * in, multilog_t * log);
int dada_bind_thread_to_core (int core);

int dada_bind_thread_to_core(int core)
{

  cpu_set_t set;
  pid_t tpid;

  CPU_ZERO(&set);
    CPU_SET(core, &set);
      tpid = syscall(SYS_gettid);

  if (sched_setaffinity(tpid, sizeof(cpu_set_t), &set) < 0) {
    printf("failed to set cpu affinity: %s", strerror(errno));
    return -1;
  }
  
  CPU_ZERO(&set);
  if ( sched_getaffinity(tpid, sizeof(cpu_set_t), &set) < 0 ) {
    printf("failed to get cpu affinity: %s", strerror(errno));
    return -1;
  }

  return 0;
}

void dsaX_dbgpu_cleanup (dada_hdu_t * in, multilog_t * log)
{
  
  if (dada_hdu_unlock_read (in) < 0)
    {
      multilog(log, LOG_ERR, "could not unlock read on hdu_in\n");
    }
  dada_hdu_destroy (in);
  
}

void usage()
{
  fprintf (stdout,
	   "dsaX_dbdisk [options]\n"
	   " -c core   bind process to CPU core\n"
	   " -b blocksize for input [default 2415919104]\n"
	   " -h print usage\n");
}


int main (int argc, char *argv[]) {
  /* DADA Header plus Data Unit */
  dada_hdu_t* hdu_in = 0;
  
  /* DADA Logger */
  multilog_t* log = 0;
  
  // input data block HDU key
  //key_t in_key = 0x0000dcda; // for dsa5
  key_t in_key = 0x0000dada;
  // command line arguments
  uint64_t blocksize = 805306368;
  //uint64_t bout = 8192*2048; // output block size - assume input is a multiple.
  uint64_t bout = 6291456;
  int core = -1;
  int arg=0;
  int machine_id = 0;
  int telescope_id = 0;
  double refdm = 25.;
  int nchans = 2048;
  double fch1 = 1280.;
  double foff = 0.122;
  double tstart = 50000;
  double tsamp = 5e-9;
  int nifs = 1;
  char  source_name[20] = "zenith";

  while ((arg=getopt(argc,argv,"c:b:h")) != -1)
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
	      printf ("ERROR: -c flag requires argument\n");
	      return EXIT_FAILURE;
	    }
	case 'b':
	  blocksize = (uint64_t)(atoi(optarg));
	  break;
	case 'h':
	  usage();
	  return EXIT_SUCCESS;
	}
    }

  // DADA stuff
  
  log = multilog_open ("dsaX_dbdisk", 0);
  multilog_add (log, stderr);

  // open connection to the in/read DB
  hdu_in  = dada_hdu_create (log);
  dada_hdu_set_key (hdu_in, in_key);
  if (dada_hdu_connect (hdu_in) < 0) {
    printf ("dsaX_correlator_copy: could not connect to input buffer\n");
    return EXIT_FAILURE;
  }
  if (dada_hdu_lock_read (hdu_in) < 0) {
    printf ("dsaX_correlator_copy: could not lock to input buffer\n");
    return EXIT_FAILURE;
  }
  // Bind to cpu core
  if (core >= 0)
    {
      printf("binding to core %d\n", core);
      if (dada_bind_thread_to_core(core) < 0)
	printf("dsaX_correlator_copy: failed to bind to core %d\n", core);
    }

  bool observation_complete=0;

  // stuff for writing data
  char * cpbuf = (char *)malloc(sizeof(char)*blocksize);
  char * outbuf = (char *)malloc(sizeof(char)*bout);
  int ngulps = (int)(blocksize/bout);
  int gulp = 0, wseq = 0;;
  char *in_data;
  uint64_t written=0, written2=0;
  uint64_t block_id, bytes_read=0;
  FILE *fout;
  char fnam[100];
  

  multilog(log, LOG_INFO, "main: have ngulps %d, blocksize %llu, bout %llu\n",ngulps,blocksize,bout);
  
  // more DADA stuff - deal with headers
  
  uint64_t header_size = 0;
   multilog(log, LOG_INFO, "0\n");
   //multilog(log,LOG_INFO, "%s\n",hdu_in->header_block);
  // read the header from the input HDU
  char * header_in = ipcbuf_get_next_read (hdu_in->header_block, &header_size);
   multilog(log, LOG_INFO, "1\n");
  if (!header_in)
    {
      multilog(log ,LOG_ERR, "main: could not read next header\n");
      dsaX_dbgpu_cleanup (hdu_in, log);
      return EXIT_FAILURE;
    }
  multilog(log, LOG_INFO, "2\n");
  // mark the input header as cleared
  if (ipcbuf_mark_cleared (hdu_in->header_block) < 0)
    {
      multilog (log, LOG_ERR, "could not mark header block cleared [input]\n");
      dsaX_dbgpu_cleanup (hdu_in, log);
      return EXIT_FAILURE;
    }
multilog(log, LOG_INFO, "3");
  // main reading loop

  multilog(log, LOG_INFO, "main: starting read\n");

  while (!observation_complete) {

    // read a DADA block
    in_data = ipcio_open_block_read (hdu_in->data_block, &bytes_read, &block_id);
    // copy
    memcpy(cpbuf, in_data, blocksize);
    multilog(log, LOG_INFO, "main: starting new write (seq %d)\n",wseq);

    // open file for writing
    sprintf(fnam,"/home/user/fl_%d.filterbank",wseq);
    multilog(log, LOG_INFO, "10\n");
    fout = fopen(fnam,"wb");
output = fout;
                send_string("HEADER_START");
                send_string("source_name");
                send_string(source_name);
                send_int("machine_id",machine_id);
                send_int("telescope_id",telescope_id);
                if (nchans > 1) {
                        send_int("data_type",1);
                } else {
                        send_int("data_type",2);
                        send_double("refdm",psrdm);
                }
                send_double("fch1",fch1);
                send_double("foff",foff);
                send_int("nchans",nchans);
                send_int("nbits",nbits);
                send_double("tstart",tstart);
                send_double("tsamp",tsamp);
                send_int("nifs",nifs);
                send_string("HEADER_END");

    multilog(log, LOG_INFO, "11\n");
    for (gulp=0;gulp<ngulps;gulp++) {
    multilog(log, LOG_INFO, "12\n");
      // copy to outbuf
      memcpy(outbuf, cpbuf+gulp*bout, bout);
    multilog(log, LOG_INFO, "13\n");
      // write
      fwrite(outbuf, 1, bout, fout);
    multilog(log, LOG_INFO, "14\n");
    }
    fclose(fout);
    multilog(log, LOG_INFO, "15\n");
    wseq++;
    multilog(log, LOG_INFO, "main: finished new write to file %s\n",fnam);
    
    // for exiting
    if (bytes_read < blocksize) {
      observation_complete = 1;
      multilog(log, LOG_INFO, "main: finished, with bytes_read %llu < expected %llu\n", bytes_read, blocksize);
    }

    // close block for reading
    ipcio_close_block_read (hdu_in->data_block, bytes_read);

  }
  
  free(cpbuf);
  free(outbuf);
  dsaX_dbgpu_cleanup (hdu_in, log);
  
}
  
