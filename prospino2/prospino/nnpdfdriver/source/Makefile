
all: testcode.cc NNPDFDriver.cc testcode.f
	g++ -O3 -c NNPDFDriver.cc
#~ 	g++ -O3 -c testcode.cc -I /home/home4/institut_1b/dmeuser/master_code/LHAPDF_old/LHAPDF-6.1.6/../local/include
#~ 	g++ -O3 testcode.o NNPDFDriver.o -o testcode -I /home/home4/institut_1b/dmeuser/master_code/LHAPDF_old/LHAPDF-6.1.6/../local/include -L /home/home4/institut_1b/dmeuser/master_code/LHAPDF_old/LHAPDF-6.1.6/../local/lib -lLHAPDF
	gfortran -O3 -c testcode.f -o testcode.o
#~ 	gfortran -O3 testcode.o NNPDFDriver.o -o testcode_fortran -lstdc++
clean:
	rm -rf *.o *.fo testcode testcode_fortran *~
