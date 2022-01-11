### Flags and directories
UNSUPPORTED_FLAG=-Ilib/eigen-3.4.0/unsupported
EIGEN_FLAG=-Ilib/eigen-3.4.0
BOOST_FLAG=-Ilib/boost_1_77_0

OUTDIR=output
CC=g++
CPPFLAGS=-std=c++17 -msse2 -O2 $(UNSUPPORTED_FLAG) $(EIGEN_FLAG) $(BOOST_FLAG)
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

build : $(OUTDIR) $(OBJS)
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


### Prerequisites
$(OUTDIR) :
	@mkdir $@

$(OUTDIR)/$(PRODUCT).o : $(PRODUCT_NAME).cpp
	$(CC) -c -o $@ $< $(CPPFLAGS)

$(OUTDIR)/%.o : src/%.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS)
