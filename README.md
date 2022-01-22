# Raman-Scattering-Code-Conversion
C++ code written as part of the UQ SMP project "Efficient code for calculation of Raman scattering by spheroids"

Author is Siwan Li

Project supervisor is A/Prof Taras Plakhotnik.

## Introduction

The original MATLAB code was written to calculate Raman scattering by spheroids with arbitrary precision. The original code (shared with me by T. Plakhotnik) uses the SMARTIES v1.01 MATLAB package, with some modifications made to support arbitrary precision using the Advanpix Multiprecision Toolbox. SMARTIES is an implementation of the T-matrix/Extended Boundary-Condition Method for light-scattering by spheroids. The goal of the project is to rewrite the code that was used for the calculation of Raman scattering in a more efficient compiled programming language that can support numerical computation. Hopefully this will accelerate the computation speeds by orders of magnitude and allow calculations for larger spheroids using higher precision.

## Instructions

Open the terminal to the Raman-Scattering-Code-Conversion directory. For calculating using C++'s inbuilt floating point types, simply run the following command to compile the program:
```
make [OPTIONS]
```

`[OPTIONS]` can contain any of the following parameters:
- `mp`: the program will have the option to do calculations using customised arbitrary-precision floating points. By default, 113 bits of precision are used. GMP and MPFR must be manually installed for this option to work. (To install them, download both libraries using the links in the 'Dependencies' section below, extract the files and follow the instructions given in their respective 'INSTALL' files.)
- `SUBTHREADS=<value>`: within each top-level thread, the program will perform the inner calculations using `<value>` many sub-threads. By default, if this option is not specified, the number of sub-threads used is 1.
- `PRECISION=<value>`: if used with `mp`, the program will customise the custom arbitrary type to use `<value>` many bits. The precision allowed is only limited by the amount of memory available on your computer.

Note that to change the number of threads or bits of precision used, the program must be recompiled. In either case, once the binaries have been built (this may take a few minutes), enter the following command to run the program:
```
output/raman_elastic_scattering
```
By default, the program also writes the output to `output/results.txt`.

You can change the calculation parameters by editing the `config.txt` file. Parameters must be entered in the format `Parameter:<value>`, where `<value>` is either a number or a fully lower-case word, no spaces. If no value is given, the program uses the default value.

Here's what each 'run' parameter under does:
- `No. of CPUs` is the number of CPUs (threads) to be used to run the calculations. Specifically, it's the radii of the particles to calculate for is allocated among `<value>` many CPUs. Within each of these threads, more threads can be made to calculate the inner-most for loops of the code in parallel (this is controlled using the `SUBTHREADS=<value>` option above). If `<value>` is 0, the program is run with the maximum number of CPUs. By default, this parameter is 1.
- `Minimum diameter` is the minimum diameter in nanometres of the largest axis of the spheroids to calculate the scattering of. By default, this parameter is 1000.
- `Maximum diameter` is the maximum diameter in nanometres of the largest axis of the spheroids to calculate the scattering of. By default, this parameter is 2000.
- `No. of radii` is the number of radii between `DIA_MIN`/2 and `DIA_MAX`/2 inclusive to calculate with. The lengths of these radii are regularly spaced. By default, this parameter is 100.
- `No. of thetas` is the number of spheroid orientations to calculate, with angles between 0 degrees and 90 degrees inclusive. The angles of these orientations are regularly spaced. By default, this parameter is 19.
- `Calculation type` is the type of precision used for floating-point calculations. By default, this parameter is `double`. The value of this parameter must be of the form `XXXX` to run all calculations with type `XXXX` or `XXXX-YYYY` to calculate T-matrices with type `XXXX` and integrals with type `YYYY`, where `XXXX` and `YYYY` can be one of the following calculation types:
  - `single` for single-precision floating points (using C++'s float types).
  - `double` for double-precision floating points.
  - `quad` for quadruple-precision floating points (using C++'s long double types).
  - `custom` for performing all calculations with customised floating points of arbitrary precision. (Only works if code was compiled using `make mp`)
- `Print output to file` is a `yes`/`no` parameter. If `yes`, the program writes the output to `output/results.txt`. Otherwise, it doesn't write to any file.

## Progress

  - [x] aux* functions (2/2)
  - [x] vsh* functions (9/9)
  - [x] sph* functions (11/11)
  - [x] rvh* functions (5/3)
  - [x] slv* functions (2/2)
  - [x] pst* functions (2/2)
  - [x] Non-SMARTIES functions (3/3)
  - [x] Minimum viable product
  - [x] Arbitrary-precision support
  - [x] Automated parallel computing

## Notable differences from SMARTIES

- Expansion coefficients and other arrays that store using p-indices now include values for when n=0, m=0. This makes P (the length of the p-vectors) equal to (N+1)^2 and makes the code more convenient in a index-by-zero language without compromising any calculations. As a side effect, many other functions that calculate values for various values n are affected as well. Functions that are affected by this include the following. (Note that in almost all cases, the extra values are 0, so they act as padding.)
  - vshPinmTaunm
  - vshGetIncidentCoeffs
  - vshEgenThetaAllPhi
  - vshGetZnAll
  - sphGetBesselProductsPrimes
  - sphGetModifiedBesselProducts
  - sphCalculatePQ
  - rvhGetFieldCoefficients
  - pstScatteringMatrixOA
- Functions that take an argument of type string (for specifying method of calculation) now take them as enum types instead. This causes these functions to theoretically behave slightly differently for unexpected values of the argument, although such behaviour shouldn't be able to manifest at compile-time. Functions that are affected by this include:
  - auxPrepareIntegrals
  - vshMakeIncidentParams
  - vshEgenThetaAllPhi
  - vshGetZnAll
- Some structs that originally contained [1 x X] matrix members may now contain [X x 1] matrix members instead for slightly better syntax. Structs that have this include (but are not limited to):
  - stAbcdnm
- Currently, auxPrepareIntegrals doesn't read from any pre-calculated values when preparing integrals.
- sphCalculatePQ doesn't try to access stParams.output, since stParams is expected to be the same struct type as the one given in the specification, in which case output is not a member of stParams. Such a member does exist in stOptions however, so future implementations may take stOptions as an argument type. (For calculating Raman scattering, this option is true by default.)
- The slvGetOptionsFromStruct function cannot be called on its own and is instead implemented into the stOptions constructors.
- The pstMakeStructForField function currently puts stIncPar into returning struct stRes by using std::move(); this means that the stIncPar will be made empty after pstMakeStructForField is used. This is done, since in the final program, the input stIncPar doesn't need to be kept.
- sphEstimateDelta doesn't get called in slvForT in the final program for this, so this function isn't implemented.
- rvhGetSymmetricMat is not debugged since it's never used in the final raman_elastic_scattering program.
- rvhTruncateMatrices, rvhGetSymmetricMat, rvhGetFieldCoefficients, rvhGetAverageCrossSections can only take stTR vectors as arguments and not stPR vectors. (Overloading these functions with versions that can take stPR vectors is relatively simple however.)

## Dependencies

- A compiler with C++17 support (tested with gcc 9.3.0-17 and clang 10.0.0-4 targeting x86_64-pc-linux-gnu)
- Eigen 3.4.0, a free, open-source, efficient and comprehensive linear algebra library made for C++. Website: https://eigen.tuxfamily.org/index.php
- Boost (C++ Libraries) 1.77.0, a free, open-source set of libraries with applications in a wide variety of areas. Used for some template maths functions and its MPFR class wrapper. Website: https://www.boost.org
- The GNU Multi Precision Arithmetic Library (GMP), a C library that provides support for arbitrary precision arithmetic. While it has arbitrary-precision floating-point types, MPFR is preferred. This library is a prerequisite for MPFR. Website: https://gmplib.org/
- The GNU Multiple Precision Floating-Point Reliable Library (GNU MPFR), a C library that provides support for arbitrary-precision floating-point computation. Website: https://www.mpfr.org

## Notes

- Although most of the code is based on the SMARTIES v1.01 MATLAB package, only the functions necessary for calculating Raman scattering will be converted.
- The coding style is loosely based on Google's C++ style guide, with variable, struct and function names matching closely with the ones given in the original SMARTIES code.
- The code was written to be compiled with GCC and on hardware that correctly implements the IEEE 754 standard. If IEEE 754 is not supported, compiling or executing the code may cause a divide-by-zero hardware exception.

## Copyright Disclaimer

The following libraries are open-source and were written by their respective developers under their respective licenses. Neither S. Li (myself), nor T. Plakhotnik wrote any of the following code.

SMARTIES was written by Walter Somerville, Baptiste Auguié, and Eric Le Ru (copyright 2015). The package is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License. You may find the original SMARTIES package here:
https://www.wgtn.ac.nz/scps/research/research-groups/raman-lab/numerical-tools/smarties

The code in `tensor_matrix_cast.hpp` was written by DavidAce on the following stackoverflow thread: https://stackoverflow.com/questions/48795789/eigen-unsupported-tensor-to-eigen-matrix (accessed 17 Dec 2021)

Parts of the Eigen 3.4.0 library have been included under the MPL2 license and parts of the Boost 1.77.0 library have been included under the Boost license. Eigen and Boost were made by their respective developers listed on their websites.
