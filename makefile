### Flags and directories
UNSUPPORTED_FLAG=-Ilib/eigen-3.4.0/unsupported
EIGEN_FLAG=-Ilib/eigen-3.4.0
BOOST_FLAG=-Ilib/boost_1_77_0
MP_LIBS=-lmpfr -lgmp

OUTDIR=output
CC=g++
OBJ_FLAGS=-std=c++17 -Wall -msse2 -O2 -ftree-parallelize-loops=4 $(UNSUPPORTED_FLAG) $(EIGEN_FLAG) $(BOOST_FLAG)
LINK_FLAGS=-fopenmp
PRODUCT_NAME=raman_elastic_scattering
PRODUCT=$(OUTDIR)/$(PRODUCT_NAME)
DEF_OBJS=$(OUTDIR)/math.o $(OUTDIR)/defs_default.o $(OUTDIR)/defs_quad.o
MP_OBJS=$(OUTDIR)/mp_part1.o $(OUTDIR)/mp_part2.o $(OUTDIR)/mp_part3.o $(OUTDIR)/mp_part4.o


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
	@echo "COMPILER FLAGS:" $(OBJ_FLAGS)
	@echo "LIBRARIES:" $(MP_LIBS)
	@echo "PRODUCT:" $(PRODUCT)
	@echo "DEFAULT OBJS:" $(DEF_OBJS)
	@echo "MP_OBJS:" $(MP_OBJS)


### Prerequisites
$(OUTDIR) :
	@mkdir $@

$(OUTDIR)/main.o : main.cpp
	$(CC) -c -o $@ $< $(OBJ_FLAGS)

$(OUTDIR)/main_mp.o : main_mp.cpp
	$(CC) -c -o $@ $< $(OBJ_FLAGS) $(MP_LIBS)

$(OUTDIR)/%.o : src/%.cpp
	$(CC) -c -o $@ $< $(OBJ_FLAGS)

$(OUTDIR)/%.o : src_mp/%.cpp
	$(CC) -c -o $@ $< $(OBJ_FLAGS) $(MP_LIBS)
