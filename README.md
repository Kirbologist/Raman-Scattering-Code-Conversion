# Raman-Scattering-Code-Conversion
C++ code written as part of the UQ SMP project "Efficient code for calculation of Raman scattering by spheroids"

Project supervisor is A/Prof Taras Plakhotnik.

## Introduction

The original MATLAB code was written to calculate Raman scattering by spheroids with arbitrary precision. The original code (lent to me by T. Plakhotnik) uses the SMARTIES v1.01 MATLAB package, with some modifications made to support multiprecision using the Advanpix Multiprecision Toolbox. SMARTIES is an implementation of the T-matrix/Extended Boundary-Condition Method for light-scattering by spheroids. The goal of the project is to rewrite the code that was used for the calculation of Raman scattering in a more efficient compiled programming language that can support numerical computation. Hopefully this will accelerate the computation speeds by orders of magnitude and allow calculations for larger spheroids using higher precision.

## Progress

  - [x] aux* functions (2/2)
  - [ ] vsh* functions (8/9)
  - [x] sph* functions (11/11)
  - [ ] rvh* functions (0/5)
  - [ ] slv* functions (0/2)
  - [ ] main function (0/1)
  - [ ] arbitrary-precision support
  - [ ] parallel computing support

## Notable differences from smarties

- Expansion coefficients and other arrays that store using p-indices now include values for when n=0, m=0. This makes P (the length of the p-vectors) equal to (N+1)^2 and makes the code more convenient in a index-by-zero language without compromising any calculations. As a side effect, many other functions that only calculate values that depend only on n are affected as well. Functions that are affected by this include:
  - vshPinmTaunm
  - vshGetIncidentCoeffs
  - vshEgenThetaAllPhi
  - vshGetZnAll
  - sphGetBesselProductsPrimes
  - sphGetModifiedBesselProducts
  - sphCalculatePQ
- Functions that take 'scheme/mode' argument of type string now take them as enum types instead. This causes these functions to theoretically behave slightly differently for non-standard values of the 'scheme/mode' argument, although such behaviour shouldn't be able to manifest at compile-time. Functions that are affected by this include:
  - auxPrepareIntegrals
  - vshMakeIncidentParams
  - vshEgenThetaAllPhi
  - vshGetZnAll
- sphGetBesselProductsPrimes has its signature changed so that it now the int N is an argument, where N is the the number of Eigen Arrays in prods - 2, or the original int N_max that was used to produce prods. This is because of the limitations of C++ arrays. This issue would be mitigated if the array of Eigen Arrays can easily be implemented using Eigen's tensor module instead.
- Currently, auxPrepareIntegrals doesn't read from any pre-calculated values when preparing integrals.

## Dependencies

- Eigen 3.4.0, a free, open-source, efficient and comprehensive linear algebra library made for C++. Website: https://eigen.tuxfamily.org/index.php
- Boost (C++ Libraries) 1.77.0, a free, open-source set of libraries with applications in a wide variety of areas. Currently only used for some special mathematical functions. Website: https://www.boost.org
- The GNU Multiple Precision Floating-Point Reliable Library (GNU MPFR), a C library that provides support for arbitrary-precision floating-point computation. Website: https://www.mpfr.org

## Notes

- Although most of the code is based on the SMARTIES v1.01 MATLAB package, only the functions necessary for calculating Raman scattering will be converted.
- The coding style is loosely based on Google's C++ style guide, with variable, struct and function names matching closely with the ones given in the original SMARTIES code, which themselves closely match the names given in "Scattering, Absorption, and Emission of Light by Small Particles" by Mishchenko et al.
- The code was written to be compiled with G++ and on hardware that correctly implements the IEEE 754 standard. If IEEE 754 is not supported, compiling or executing the code may cause a divide-by-zero hardware exception.

## Copyright Disclaimer

SMARTIES was written by Walter Somerville, Baptiste Auguié, and Eric Le Ru (copyright 2015); neither me, nor T. Plakhotnik wrote any of the original code. The package is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License. You may find the original SMARTIES package here:
https://www.wgtn.ac.nz/scps/research/research-groups/raman-lab/numerical-tools/smarties

The code in TensorMatrixCast.h was written by DavidAce on the following stackoverflow thread: https://stackoverflow.com/questions/48795789/eigen-unsupported-tensor-to-eigen-matrix (accessed 17 Dec 2021)
