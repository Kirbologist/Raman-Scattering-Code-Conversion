### Flags and directories
UNSUPPORTED_FLAG=-Ilib/eigen-3.4.0/unsupported
EIGEN_FLAG=-Ilib/eigen-3.4.0
BOOST_FLAG=-Ilib/boost_1_77_0
MP_LIBS=-lmpfr -lgmp

THREADS=1
PRECISION=113 # default value of 113 is number of significand bits in quadruple-precision

OUTDIR=output
CC=g++
CPPFLAGS=-std=c++17 -Wall -msse2 -O2 -ftree-parallelize-loops=$(THREADS) $(UNSUPPORTED_FLAG) $(EIGEN_FLAG) $(BOOST_FLAG)
MP_FLAGS=-DPRECISION=$(PRECISION)
LINK_FLAGS=-fopenmp
PRODUCT_NAME=raman_elastic_scattering
PRODUCT=$(OUTDIR)/$(PRODUCT_NAME)
DEF_OBJS=$(OUTDIR)/math.o $(OUTDIR)/defs_default.o $(OUTDIR)/defs_quad.o
MP_OBJS=$(OUTDIR)/defs_mp1.o $(OUTDIR)/defs_mp2.o $(OUTDIR)/defs_mp3.o $(OUTDIR)/defs_mp4.o


### Recipes
.PHONY=build clean info

build : $(OUTDIR) $(DEF_OBJS) $(OUTDIR)/main.o
	$(CC) -o $(PRODUCT) $(DEF_OBJS) $(OUTDIR)/main.o $(LINK_FLAGS)

mp : $(OUTDIR) $(DEF_OBJS) $(MP_OBJS) $(OUTDIR)/main_mp.o
	$(CC) -o $(PRODUCT) $(DEF_OBJS) $(MP_OBJS) $(OUTDIR)/main_mp.o $(LINK_FLAGS) $(MP_LIBS)

clean :
	@rm -f $(DEF_OBJS)
	@rm -f $(MP_OBJS)
	@rm -f $(OUTDIR)/main.o
	@rm -f $(OUTDIR)/main_mp.o
	@rm -f $(PRODUCT)

clean-mp :
	@rm -f $(MP_OBJS)
	@rm -f $(OUTDIR)/main_mp.o

info :
	@echo "UNSUPPORTED EIGEN FLAG:" $(UNSUPPORTED_FLAG)
	@echo "EIGEN FLAG:" $(EIGEN_FLAG)
	@echo "BOOST FLAG:" $(BOOST_FLAG)
	@echo "COMPILER:" $(CC)
	@echo "COMPILER FLAGS:" $(CPPFLAGS)
	@echo "COMPILER MP FLAGS:" $(MP_FLAGS)
	@echo "LIBRARIES:" $(MP_LIBS)
	@echo "PRODUCT:" $(PRODUCT)
	@echo "DEFAULT OBJS:" $(DEF_OBJS)
	@echo "MP OBJS:" $(MP_OBJS)


### Prerequisites
$(OUTDIR) :
	@mkdir $@

$(OUTDIR)/main.o : main.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS)

$(OUTDIR)/main_mp.o : main_mp.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS) $(MP_FLAGS) $(MP_LIBS)

$(OUTDIR)/%.o : src/%.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS)

$(OUTDIR)/%.o : src_mp/%.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS) $(MP_FLAGS) $(MP_LIBS)
