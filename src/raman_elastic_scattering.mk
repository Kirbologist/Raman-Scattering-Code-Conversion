
UNSUPPORTED_FLAG=-I../lib/eigen-3.4.0/unsupported
EIGEN_FLAG=-I../lib/eigen-3.4.0
BOOST_FLAG=I../lib/boost_1_77_0

CC=g++
CPPFLAGS=-std=c++17 -Wall -pedantic -g -msse2 -O2 $(UNSUPPORTED_FLAG) $(EIGEN_FLAG) $(BOOST_FLAG)
OBJS=raman_elastic_scattering.o math.o raman_aux.o vsh.o sph.o rvh.o slv.o pst.o
PRODUCT=raman_elastic_scattering

.PHONY=build clean info

build : $(OBJS)
	$(CC) -o $(PRODUCT) $(OBJS)

clean :
	@rm -f $(OBJS)
	@rm -f $(PRODUCT)

info :

%.o : $.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS)
