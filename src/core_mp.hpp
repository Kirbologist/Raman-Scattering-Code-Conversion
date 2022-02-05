/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This code contains all common or basic header files and typedefs used for arbitrary precision functionality.
*/

#ifndef CORE_MP_HPP
#define CORE_MP_HPP

#include <boost/multiprecision/mpfr.hpp>

using namespace boost::multiprecision;

typedef number<mpfr_float_backend<PRECISION, mpfr_allocation_type::allocate_stack>> RamanFloat;

#endif
