#ifndef CORE_MP_HPP
#define CORE_MP_HPP

#include <boost/multiprecision/mpfr.hpp>

using namespace boost::multiprecision;

const int precision = 34;

typedef number<mpfr_float_backend<precision, mpfr_allocation_type::allocate_stack>> raman_float;

#endif
