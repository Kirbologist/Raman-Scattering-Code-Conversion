#ifndef CORE_MP_HPP
#define CORE_MP_HPP

#include <boost/multiprecision/mpfr.hpp>

using namespace boost::multiprecision;

typedef number<mpfr_float_backend<PRECISION, mpfr_allocation_type::allocate_stack>> raman_float;

#endif
