UNSUPPORTED_FLAG=-Ilib/eigen-3.4.0/unsupported
EIGEN_FLAG=-Ilib/eigen-3.4.0
BOOST_FLAG=-Ilib/boost_1_77_0

CC=g++
CPPFLAGS=-std=c++17 -msse2 -O2 $(UNSUPPORTED_FLAG) $(EIGEN_FLAG) $(BOOST_FLAG)
OBJS=raman_elastic_scattering.o math.o smarties_aux.o vsh.o sph.o rvh.o slv.o pst.o
PRODUCT=raman_elastic_scattering

.PHONY=build clean info

build : $(OBJS)
	$(CC) -o $(PRODUCT) $(OBJS)

clean :
	@rm -f $(OBJS)
	@rm -f $(PRODUCT)

info :
	@echo "UNSUPPORTED EIGEN FLAG:" $(UNSUPPORTED_FLAG)
	@echo "EIGEN FLAG:" $(EIGEN_FLAG)
	@echo "BOOST FLAG:" $(BOOST_FLAG)
	@echo "COMPILER:" $(CC)
	@echo "COMPILER FLAGS:" $(CPPFLAGS)
	@echo "OBJECT FILES:" $(OBJS)
	@echo "PRODUCT:" $(PRODUCT)

$(PRODUCT).o : $(PRODUCT).cpp
	$(CC) -c -o $@ $< $(CPPFLAGS)

%.o : src/%.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS)
