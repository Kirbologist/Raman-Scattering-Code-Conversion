### Flags and directories
UNSUPPORTED_FLAG=-Ilib/eigen-3.4.0/unsupported
EIGEN_FLAG=-Ilib/eigen-3.4.0
BOOST_FLAG=-Ilib/boost_1_77_0
MP_LIBS=-lmpfr -lgmp

PRECISION=113 # default value of 113 is number of significand bits in quadruple-precision

BUILDDIR=build
CC=g++
CPPFLAGS=-std=c++17 -Wall -msse2 -O2 -fopenmp
CPPFLAGS+=$(UNSUPPORTED_FLAG) $(EIGEN_FLAG) $(BOOST_FLAG)
MP_FLAGS=-DPRECISION=$(PRECISION)
PRODUCT=raman_elastic_scattering
UTILS_PRODUCT=store_GL_quadrature
DEF_OBJS=$(BUILDDIR)/misc.o $(BUILDDIR)/raman_elastic_scattering.o
DEF_OBJS+=$(BUILDDIR)/defs_single.o $(BUILDDIR)/defs_double.o $(BUILDDIR)/defs_quad.o
MP_OBJS=$(BUILDDIR)/defs_custom.o

### Recipes
.PHONY=regular clean info

all : mp utils

regular : $(BUILDDIR) $(DEF_OBJS) $(BUILDDIR)/main.o
	$(CC) -o $(PRODUCT) $(DEF_OBJS) $(BUILDDIR)/main.o $(CPPFLAGS)

mp : $(BUILDDIR) $(DEF_OBJS) $(MP_OBJS) $(BUILDDIR)/main_mp.o
	$(CC) -o $(PRODUCT) $(DEF_OBJS) $(MP_OBJS) $(BUILDDIR)/main_mp.o $(CPPFLAGS) $(MP_LIBS)

utils: $(BUILDDIR) $(DEF_OBJS) $(BUILDDIR)/utils.o
	$(CC) -o $(UTILS_PRODUCT) $(DEF_OBJS) $(BUILDDIR)/utils.o $(CPPFLAGS)

clean :
	@rm -f $(DEF_OBJS)
	@rm -f $(MP_OBJS)
	@rm -f $(BUILDDIR)/main.o
	@rm -f $(BUILDDIR)/main_mp.o
	@rm -f $(BUILDDIR)/utils.o
	@rm -f $(UTILS_PRODUCT)
	@rm -f $(PRODUCT)

clean-mp :
	@rm -f $(MP_OBJS)
	@rm -f $(BUILDDIR)/main_mp.o

info :
	@echo "UNSUPPORTED EIGEN FLAG:" $(UNSUPPORTED_FLAG)
	@echo "EIGEN FLAG:" $(EIGEN_FLAG)
	@echo "BOOST FLAG:" $(BOOST_FLAG)
	@echo "COMPILER:" $(CC)
	@echo "COMPILER REGULAR FLAGS:" $(CPPFLAGS)
	@echo "COMPILER MP FLAGS:" $(CPPFLAGS) $(MP_FLAGS)
	@echo "LIBRARIES:" $(MP_LIBS)
	@echo "PRODUCTS:" $(PRODUCT) $(UTILS_PRODUCT)
	@echo "REGULAR OBJS:" $(DEF_OBJS)
	@echo "MP OBJS:" $(DEF_OBJS) $(MP_OBJS)


### Prerequisites
$(BUILDDIR) :
	@mkdir $@

$(BUILDDIR)/defs_custom.o : src/defs_custom.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS) $(MP_FLAGS) $(MP_LIBS)

$(BUILDDIR)/main_mp.o : src/main_mp.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS) $(MP_FLAGS) $(MP_LIBS)

$(BUILDDIR)/%.o : src/%.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS)
