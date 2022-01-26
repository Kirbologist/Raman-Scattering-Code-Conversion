### Flags and directories
UNSUPPORTED_FLAG=-Ilib/eigen-3.4.0/unsupported
EIGEN_FLAG=-Ilib/eigen-3.4.0
BOOST_FLAG=-Ilib/boost_1_77_0
STATIC_FLAG=-include /home/siwanli/glibc_version_header/version_headers/x64/force_link_glibc_2.17.h
MP_LIBS=-l:libmpfr.a -l:libgmp.a

PRECISION=113 # default value of 113 is number of significand bits in quadruple-precision

BUILDDIR=build
CC=g++
CPPFLAGS=-std=c++17 -Wall -msse2 -O2 -fopenmp -static-libgcc -static-libstdc++
CPPFLAGS+=$(STATIC_FLAG) $(UNSUPPORTED_FLAG) $(EIGEN_FLAG) $(BOOST_FLAG)
MP_FLAGS=-DPRECISION=$(PRECISION)
LINK_FLAGS=-fopenmp
PRODUCT_NAME=raman_elastic_scattering
PRODUCT=$(PRODUCT_NAME)
DEF_OBJS=$(BUILDDIR)/math.o $(BUILDDIR)/raman_elastic_scattering.o
DEF_OBJS+=$(BUILDDIR)/defs_single.o $(BUILDDIR)/defs_double.o $(BUILDDIR)/defs_quad.o
MP_OBJS=$(BUILDDIR)/defs_custom.o

### Recipes
.PHONY=regular clean info

regular : $(BUILDDIR) $(DEF_OBJS) $(BUILDDIR)/main.o
	$(CC) -o $(PRODUCT) $(DEF_OBJS) $(BUILDDIR)/main.o $(CPPFLAGS)

mp : $(BUILDDIR) $(DEF_OBJS) $(MP_OBJS) $(BUILDDIR)/main_mp.o
	$(CC) -o $(PRODUCT) $(DEF_OBJS) $(MP_OBJS) $(BUILDDIR)/main_mp.o $(CPPFLAGS) $(MP_LIBS)

clean :
	@rm -f $(DEF_OBJS)
	@rm -f $(MP_OBJS)
	@rm -f $(BUILDDIR)/main.o
	@rm -f $(BUILDDIR)/main_mp.o
	@rm -f $(PRODUCT)

clean-mp :
	@rm -f $(MP_OBJS)
	@rm -f $(BUILDDIR)/main_mp.o

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
$(BUILDDIR) :
	@mkdir $@

$(BUILDDIR)/main.o : main.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS)

$(BUILDDIR)/main_mp.o : main_mp.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS) $(MP_FLAGS) $(MP_LIBS)

$(BUILDDIR)/%.o : src/%.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS)

$(BUILDDIR)/%.o : src_mp/%.cpp
	$(CC) -c -o $@ $< $(CPPFLAGS) $(MP_FLAGS) $(MP_LIBS)
