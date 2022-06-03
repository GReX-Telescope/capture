BINDIR=/usr/local/dsaX/bin

CC=gcc
CFLAGS = -g -pthread -I/usr/local/psrdada/src
CDEPS=dsaX_correlator_def.h dsaX_correlator_udpdb_thread.h
SDEPS=dsaX_spectrometer_def.h dsaX_spectrometer_udpdb_thread.h
COBJ = dsaX_correlator_udpdb_thread.o
SOBJ = dsaX_spectrometer_udpdb_thread.o
LIBS = -L/usr/local/lib -lpsrdada

dsaX_correlator_udpdb_thread.o: dsaX_correlator_udpdb_thread.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

dsaX_correlator_udpdb_thread: $(COBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

dsaX_spectrometer_udpdb_thread.o: dsaX_spectrometer_udpdb_thread.c $(SDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

dsaX_spectrometer_udpdb_thread: $(SOBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

dsaX_spectrometer_reorder: dsaX_spectrometer_reorder.cu
	/usr/local/cuda-8.0/bin/nvcc -o $@ $^ -ccbin=g++ -I/usr/local/thrust/thrust -I/usr/local/psrdada/src -I/usr/local/include -L/usr/local/lib -lcfitsio -lpsrdada -g -O3

dsaX_spectrometer: dsaX_spectrometer.cu
	/usr/local/cuda-8.0/bin/nvcc -o $@ $^ -ccbin=g++ -I/home/user/sigproc-4.3/ -I/usr/local/thrust/thrust -I/usr/local/psrdada/src -I/usr/local/include -L/usr/local/lib -L/home/user/sigproc-4.3/ -lsigproc_linux -lpsrdada -lm -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran --gpu-architecture=sm_61 -O3

dsaX_correlator_trigger: dsaX_correlator_trigger.c
	gcc -o $@ $^ -I/usr/local/psrdada/src -I/usr/local/include -L/usr/local/lib -lpsrdada -lm -pthread -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran

dsaX_dbdisk: dsaX_dbdisk.c
	g++ -o $@ $^ -I/home/user/psrdada/src -I/usr/local/include -L/usr/local/lib -lpsrdada -lm -pthread -O3 -g -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran

dsaX_dbdisk_fb: dsaX_dbdisk_fb.c
	g++ -o $@ $^ -I/home/user/psrdada/src -I/usr/local/include -I/home/user/sigproc-4.3 -L/usr/local/lib -lpsrdada -lm -pthread -O3 -g -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran

.PHONY: clean

clean:
	rm -f *.o

install:
	cp dsaX_correlator_udpdb_thread dsaX_spectrometer_udpdb_thread dsaX_correlator_trigger dsaX_spectrometer dsaX_spectrometer_reorder dsaX_dbdisk $(BINDIR)
