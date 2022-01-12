### Flags and directories
UNSUPPORTED_FLAG=-Ilib/eigen-3.4.0/unsupported
EIGEN_FLAG=-Ilib/eigen-3.4.0
BOOST_FLAG=-Ilib/boost_1_77_0
MP_FLAGS=-lmpfr -lgmp

OUTDIR=output
CC=g++
CPP_FLAGS=-std=c++17 -Wall -pedantic -g -msse2 -O2 $(UNSUPPORTED_FLAG) $(EIGEN_FLAG) $(BOOST_FLAG)
PRODUCT_NAME=raman_elastic_scattering
PRODUCT=$(OUTDIR)/$(PRODUCT_NAME)
OBJS=$(OUTDIR)/$(PRODUCT_NAME).o
OBJS+=$(OUTDIR)/math.o
OBJS+=$(OUTDIR)/smarties_aux.o
OBJS+=$(OUTDIR)/vsh.o
OBJS+=$(OUTDIR)/sph.o
OBJS+=$(OUTDIR)/rvh.o
OBJS+=$(OUTDIR)/slv.o
OBJS+=$(OUTDIR)/pst.o


### Recipes
.PHONY=build clean info

build : $(OUTDIR) $(OBJS) $(OUTDIR)/main.o
	$(CC) -o $(PRODUCT) $(OBJS) $(OUTDIR)/main.o

mp : $(OUTDIR) $(OBJS) $(OUTDIR)/main_mp.o
	$(CC) -o $(PRODUCT) $(OBJS) $(OUTDIR)/main_mp.o $(MP_FLAGS)

clean :
	@rm -f $(OBJS)
	@rm -f $(OUTDIR)/main.o
	@rm -f $(OUTDIR)/main_mp.o
	@rm -f $(PRODUCT)

info :
	@echo "UNSUPPORTED EIGEN FLAG:" $(UNSUPPORTED_FLAG)
	@echo "EIGEN FLAG:" $(EIGEN_FLAG)
	@echo "BOOST FLAG:" $(BOOST_FLAG)
	@echo "COMPILER:" $(CC)
	@echo "COMPILER FLAGS:" $(CPP_FLAGS)
	@echo "OBJECT FILES:" $(OBJS)
	@echo "PRODUCT:" $(PRODUCT)


### Prerequisites
$(OUTDIR) :
	@mkdir $@

$(OUTDIR)/%.o : src/%.cpp
	$(CC) -c -o $@ $< $(CPP_FLAGS)

$(OUTDIR)/main.o : main.cpp
	$(CC) -c -o $@ $< $(CPP_FLAGS)

$(OUTDIR)/main_mp.o : main_mp.cpp
	$(CC) -c -o $@ $< $(CPP_FLAGS) $(MP_FLAGS)
