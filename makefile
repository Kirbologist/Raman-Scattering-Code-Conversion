### Flags and directories
UNSUPPORTED_FLAG=-Ilib/eigen-3.4.0/unsupported
EIGEN_FLAG=-Ilib/eigen-3.4.0
BOOST_FLAG=-Ilib/boost_1_77_0
MP_FLAGS=-lmpfr -lgmp

OUTDIR=output
CC=g++
CPP_FLAGS=-std=c++17 -Wall -g -msse2 -O2 $(UNSUPPORTED_FLAG) $(EIGEN_FLAG) $(BOOST_FLAG)
PRODUCT_NAME=raman_elastic_scattering
PRODUCT=$(OUTDIR)/$(PRODUCT_NAME)
DEF_OBJS=$(OUTDIR)/math.o $(OUTDIR)/defs_default.o
MP_OBJS=$(OUTDIR)/mp_part1.o $(OUTDIR)/mp_part2.o $(OUTDIR)/mp_part3.o $(OUTDIR)/mp_part4.o


### Recipes
.PHONY=build clean info

build : $(OUTDIR) $(DEF_OBJS) $(OUTDIR)/main.o
	$(CC) -o $(PRODUCT) $(DEF_OBJS) $(OUTDIR)/main.o

mp : $(OUTDIR) $(DEF_OBJS) $(MP_OBJS) $(OUTDIR)/main_mp.o
	$(CC) -o $(PRODUCT) $(DEF_OBJS) $(MP_OBJS) $(OUTDIR)/main_mp.o $(MP_FLAGS)

clean :
	@rm -f $(DEF_OBJS)
	@rm -f $(MP_OBJS)
	@rm -f $(OUTDIR)/main.o
	@rm -f $(OUTDIR)/main_mp.o
	@rm -f $(PRODUCT)

info :
	@echo "UNSUPPORTED EIGEN FLAG:" $(UNSUPPORTED_FLAG)
	@echo "EIGEN FLAG:" $(EIGEN_FLAG)
	@echo "BOOST FLAG:" $(BOOST_FLAG)
	@echo "COMPILER:" $(CC)
	@echo "COMPILER FLAGS:" $(CPP_FLAGS)
	@echo "PRODUCT:" $(PRODUCT)
	@echo "DEFAULT OBJS:" $(DEF_OBJS)
	@echo "MP_OBJS:" $(MP_OBJS)


### Prerequisites
$(OUTDIR) :
	@mkdir $@

$(OUTDIR)/main.o : main.cpp
	$(CC) -c -o $@ $< $(CPP_FLAGS)

$(OUTDIR)/main_mp.o : main_mp.cpp
	$(CC) -c -o $@ $< $(CPP_FLAGS) $(MP_FLAGS)

$(OUTDIR)/%.o : src/%.cpp
	$(CC) -c -o $@ $< $(CPP_FLAGS) $(MP_FLAGS)
