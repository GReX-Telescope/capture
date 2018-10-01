#ifndef __DSAX_DEF_H
#define __DSAX_DEF_H

#include "dada_def.h"

// default port to connect to udpdb command interface
#define DSAX_DEFAULT_UDPDB_PORT   4010
#define DSAX_DEFAULT_PWC_LOGPORT  40123
#define DSAX_DEFAULT_BLOCK_KEY 0x0000dada
#define DSAX_DEFAULT_BLOCK_KEY2 0x0000dadd

#define UDP_HEADER   8              // size of header/sequence number
#define UDP_DATA     8192           // obs bytes per packet
#define UDP_PAYLOAD  8200           // header + datasize
#define UDP_IFACE    "10.10.2.1"    // default interface

#endif 

