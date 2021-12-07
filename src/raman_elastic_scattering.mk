#EIGEN_FLAG = -I/lib/include/eigen3/Eigen
EIGEN_FLAG =
GMP_FLAG =
MPFR_FLAG =

CC = g++
CPPFLAGS = -std=gnu++17 $(EIGEN_FLAG) $(GMP_FLAG) $(MPFR_FLAG)
OBJS = raman_elastic_scattering.o math.o low_level.o
PRODUCT = raman_elastic_scattering

.PHONY = build clean info

build : $(OBJS)
	$(CC) -o $(PRODUCT) $(OBJS)

clean :
	@rm -f $(OBJS)
	@rm -f $(PRODUCT)

info :

%.o : $.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS)
