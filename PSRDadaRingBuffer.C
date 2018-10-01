/***************************************************************************
 *
 *   Copyright (C) 2012 by Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *   Modified by VR Jan 2017 for dsaX
 *
 ***************************************************************************/

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <errno.h>
#include <string.h>

#include "PSRDadaRingBuffer.h"

#include "ascii_header.h"
#include "tmutil.h"

PSRDadaRingBuffer::PSRDadaRingBuffer (key_t dada_id) : DataSource ()
{
  dada_error = false;
  printed_first_line = true;

  dada_key = dada_id;

  // create the HDU
  hdu = dada_hdu_create (0);

  // set the shared memory key
  dada_hdu_set_key(hdu, dada_id);
}

PSRDadaRingBuffer::~PSRDadaRingBuffer ()
{
  if (hdu)
  {
		dada_hdu_unlock_write (hdu);
    dada_hdu_disconnect(hdu);
	}
	hdu = 0;
}


bool PSRDadaRingBuffer::connect()
{
  if (dada_hdu_connect (hdu) < 0)
  {
    dada_error = true;
    return false;
  }

  if (dada_hdu_lock_read (hdu) < 0)
  {
    dada_error = true;
    return false;
  }
	return true;
}

bool PSRDadaRingBuffer::disconnect()
{
  if (dada_hdu_unlock_read (hdu) < 0)
  {
    dada_error = true;
    return false;
  }

  if (dada_hdu_disconnect (hdu) < 0)
  {
    dada_error = true;
    return false;
  }  
	return true;
}

bool PSRDadaRingBuffer::read_header()
{
  uint64_t header_size = 0;
  char * header;

  // Wait for the next valid header block
  header = ipcbuf_get_next_read (hdu->header_block, &header_size);
  if (!header)
  {
    cerr << "PSRDadaRingBuffer::read_header could not get next header" << endl;
    dada_error = true;
    return false;
  }

#ifdef DEBUG
  cerr << "==================================================" << endl;
  cerr << header << endl;
  cerr << "==================================================" << endl;
#endif

  // get the required header params
  if (ascii_header_get (header, "NCHAN", "%d", &nchan) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract NCHAN from header" << endl;
    dada_error = true;
  }

	float bandwidth = 0;
  if (ascii_header_get (header, "BANDWIDTH", "%f", &bandwidth) < 0)
  {
    if (ascii_header_get (header, "BW", "%f", &bandwidth) < 0)
    {
      cerr << "WARNING: PSRDadaRingBuffer::read_header could not extract BANDWIDTH from header" << endl;
      dada_error = true;
    }
  }

	float cfreq = 0;
  if (ascii_header_get (header, "CFREQ", "%f", &cfreq) < 0)
  {
    if (ascii_header_get (header, "FREQ", "%f", &cfreq) < 0)
    {
      cerr << "WARNING: PSRDadaRingBuffer::read_header could not extract CFREQ from header" << endl;
      dada_error = true;
    }
  }

	if (!dada_error)
	{
		float start_freq = cfreq - (bandwidth / 2);
		float chan_width = bandwidth / nchan;
		float first_chan_cfreq = start_freq + (chan_width / 2);

		// set the cfreq of first chan and the delta freq between chans
		f0 = first_chan_cfreq;
		df = chan_width;
	}

  if (ascii_header_get (header, "NPOL", "%d", &npol) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract NPOL from header" << endl;
    dada_error = true; 
  }

  if (ascii_header_get (header, "NBIT", "%d", &nbit) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract NBIT from header" << endl;
    dada_error = true; 
  }

  if (ascii_header_get (header, "NANT", "%d", &nant) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract NANT from header" << endl;
    dada_error = true; 
  }

  

  if (ascii_header_get (header, "TSAMP", "%f", &tsamp) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract TSAMP from header" << endl;
    dada_error = true;
  }

  beam = 0;

  char utc_start_str[64] = "";
  if (ascii_header_get (header, "UTC_START", "%s", utc_start_str) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract UTC_START from header" << endl;
    dada_error = true;
  }
  else
    utc_start = str2utctime (utc_start_str);

  stride = nchan * (nbit / 8) * nant * npol;
  spectra_rate = 1000000 / (double) tsamp;

  // convert tsamp from usecs (DADA DEFAULT) to seconds
  tsamp /= 1000000;

  if (dada_error)
    return false;

  // if we sucessfully parsed the header, mark buffer cleared
  ipcbuf_mark_cleared (hdu->header_block);
  return true;
}

// 8.192 us sampling, so 65536 samples gets close to 1 s. 
size_t PSRDadaRingBuffer::get_data(size_t nsamps, unsigned char* data)
{

  uint64_t bytes_to_read = nsamps * nchan * (nbit / 8) * nant * npol;
  int64_t  bytes_read = 0;
  
  if (ipcbuf_eod((ipcbuf_t*)hdu->data_block))
    return 0;

  cerr << "PSRDadaRingBuffer::get_data: ipcio_read for " << bytes_to_read << " bytes" << endl;

  bytes_read = ipcio_read (hdu->data_block, data, bytes_to_read);

  if (!printed_first_line)
  {
    printed_first_line = true;
    unsigned char * ptr = (unsigned char *) data;
    for (unsigned i=0; i<1024; i++)
    {
      fprintf(stderr, "data[%d]=%d\n",i, ptr[i]);
    }
  }

  if (bytes_read < 0)
  {
    cerr << "PSRDadaRingBuffer::get_data: ipcio_read error: " <<  strerror(errno) << endl;
    dada_error = true;
    return -1;
  }
  
	// actually return the number of samples read, rounding down
	size_t nsamps_read = bytes_read / (nchan * (nbit / 8));

	if (nsamps_read != nsamps)
    if (!ipcbuf_eod((ipcbuf_t*)hdu->data_block))
		  cerr << "PSRDadaRingBuffer::get_data: returing fewer nsamps than requested!" << endl;

  return nsamps_read;
}
