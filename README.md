# Raman-Scattering-Code-Conversion
C++ code written as part of the UQ SMP project "Efficient code for calculation of Raman scattering by spheroids"

Author is Siwan Li

Project supervisor is A/Prof Taras Plakhotnik.

## Introduction

The original MATLAB code was written to numerically calculate Raman scattering by spheroids with arbitrary precision. The original code (shared with me by T. Plakhotnik) uses the SMARTIES v1.01 MATLAB package, with some modifications made to support arbitrary precision using the Advanpix Multiprecision Toolbox. SMARTIES is a code suite used to calculate the optical properties of spheroidal particles. The goal of the project is to rewrite the code that was used for the calculation of Raman scattering in a more efficient compiled programming language that can support numerical computation. Hopefully this will accelerate the computation speeds by orders of magnitude and allow calculations for larger spheroids using higher precision.

Raman scattering is the inelastic scattering of light, i.e. a form of scattering where incident photons exchanges energy with the incident particle before propagating off in a different direction. Typically, this exchange of energy means that the particle gains or loses heat energy, and the incident photons change in frequency (usually heat energy is gained and frequency is lowered). Here, monochromatic light is scattering off an oblate spheroidal particle that behaves like a raindrop in the air. The program then calculates the amount of light scattered into some new frequency. Although the code is originally intended to model raindrop-like particles, minimal work should be needed to use the code for calculating the Raman scattering off of any spheroid (oblate or prolate) in any medium. Either change around the parameter settings shown below, or edit the `LoadParams` function in `raman_elastic_scattering.hpp`.

## Instructions

### Setup

Open the terminal to the Raman-Scattering-Code-Conversion directory. For calculating using C++'s inbuilt floating point types, simply run the following command to compile the `raman_elastic_scattering` binaries:
```bash
make [OPTIONS]
```

`[OPTIONS]` can contain any of the following parameters:
- `mp`: `raman_elastic_scattering` will have the option to do calculations using customised arbitrary-precision floating points. By default, 113 bits of precision are used. GMP and MPFR must be manually installed for this option to work. (To install them, download both libraries using the links in the 'Dependencies' section below, extract the files and follow the instructions given in their respective 'INSTALL' files.)
- `PRECISION=<value>`: if used with `mp`, `raman_elastic_scattering` will customise the custom arbitrary type to use `<value>` many bits. The precision allowed is only limited by the amount of memory available on your computer.

Note that to change the number of threads or bits of precision used, `raman_elastic_scattering` must be recompiled.

### Running `raman_elastic_scattering`

In either case, once the binaries have been built (this may take a few minutes), enter the following command to run `raman_elastic_scattering`:
```bash
./raman_elastic_scattering [OPTIONS]
```

`[OPTIONS]` can contain any of the following parameters:
- `--input=</path/to/file>` means that `raman_elastic_scattering` reads the parameters from the file specified by `<path/to/file>`. More information about these parameters can be found in the 'Configuration' subsection. By default, the input file is `config.txt`.
- `--output-dir=</path/to/directory` means that `raman_elastic_scattering` writes the program output to the directory specified by `<path/to/drectory>`, in parallel with printing them to the terminal. Each top-level thread prints to its own output file. By default, this directory is `output`.

### Configuration

You can change the calculation parameters by editing the `config.txt` file, or the file specified under the flag `--input=</path/to/file>`. Parameters must be entered in the format `Parameter:<value>`, where `<value>` is either a number or a fully lower-case word, no spaces. If no value is given, `raman_elastic_scattering` uses the default value. Otherwise if invalid values are given, it will fall back to the default values in the best case, or interpret the value in unexpected ways in the worst case. Examples of allowable formats for numerical values are given in its own subsection.

Here's what each 'run' parameter does:
- `No. of CPUs` is a non-negative integer representing the number of CPUs (threads) to be used to run the calculations. Some algorithms like matrix products can take advantage of the extra threads to parallelise calculations. If `<value>` is 0, `raman_elastic_scattering` will run with the maximum number of CPUs. By default, this parameter is 1.
- `No. of CPUs to partition particle radii for` is a non-negative integer. The different radii of the particles is partitioned among `<value>` threads, so each thread performs all calculations for particles of its allocated radii to its own output file. By default, this parameter is 1.
- `Calculation type` is the type used for floating-point calculations. By default, this parameter is `double`. The value of this parameter must be of the form `XXXX` to run all calculations with type `XXXX` or `XXXX-YYYY` to calculate T-matrices with type `XXXX` and integrals with type `YYYY`, where `XXXX` and `YYYY` can be one of the following calculation types:
  - `single` for single-precision floating points (using C++'s float types).
  - `double` for double-precision floating points.
  - `quad` for quadruple-precision floating points (using C++'s long double types).
  - `custom` for performing all calculations with customised floating points of arbitrary precision. (Only works if code was compiled using `make mp`).
- `Print output to file` is a `yes`/`no` parameter. If `yes`, `raman_elastic_scattering` writes the output to `output/results.txt`. Otherwise, it doesn't write to any file. By default, this parameter is `yes`.

The following parameters are for determining what sizes, shapes and orientations of spheroidal particles to calculate the Raman scattering of. These particles have parameters `a`, `c`, `h`, `theta` and `phi`. `a` is the radius of the particle along the x and y semi-axes. `c` is the radius of the particle along the z semi-axis. `h` is the ratio `a/c`. `theta` is the angle of tilt of the z semi-axis with respect to the surrounding medium (zero tilt means that the z semi-axis of the particle is aligned with the z axis of the surrounding medium). `phi` is the angle of rotation along the z semi-axis, the particles are rotationally symmetric along the z semi-axis, changing this parameter should not change the results of the calculations, outside of small rounding errors.
- `Minimum diameter` is a floating-point parameter representing the minimum diameter in nanometres of the largest semi-axis. By default, this parameter is 1000.
- `Maximum diameter` is a floating-point parameter representing the maximum diameter in nanometres of the largest semi-axis. By default, this parameter is 2000.
- `Minimum h` is a floating-point value representing the minimum value of `h`. By default, this value is 1/3.
- `Maximum h` is a floating-point value representign the maximum value of `h`. By default, this value is 1/3.
- `No. of particle radii` is a positive integer representing the number of particle radii between `Minimum diameter`/2 and `Maximum diameter`/2 inclusive to calculate with. The lengths of these radii are regularly spaced. By default, this parameter is 100.
- `No. of particle thetas` is a positive integer representing the number of different `theta` to calculate with, with `theta` between 0 radians and 𝝅/2 radians inclusive. The values of `theta` are regularly spaced. By default, this parameter is 19.
- `Particle phi` is a floating-point value representing the value of `phi` in radians. By default, this parameter is 0.
- `No. of h ratios` is a positive integer representing the number of `h` values between `Minimum h` and `Maximum h` inclusive to calculate with. These values are linearly spaced. By default, this parameter is 1.

The following parameters are for determining the spherical coordinates of the sample of points inside the particles to perform field calculations with. The physics definition of spherical coordinates are used, so `r` is the radial distance from the centre of the particle, `phi` is the azimuthal/longitudinal angle and `theta` is the polar/colatitudinal angle.
- `No. of r-coordinates` is a positive integer representing the number of radial coordinates to use within the particle when performing field calculations/integrations. By default, this parameter is 100.
- `No. of phi-coordinates` is a positive integer representing the number of azimuthal angle coordinates to use within the particle when performing field calculations/integrations. By default, this parameter is 320.
- `No. of theta-coordinates` is a positive integer representing the number of polar angle/colatitude coordinates to use within the particle when performing field calculations/integrations. By default, this parameter is 320.

The following parameters are used for calculating the T-matrices.
- `Nb_theta` is a positive integer, which is the number of angles used in Gaussian quadratures for the evaluation of P and Q matrix integrals.
- `Nb_theta_pst` is a positive integer, which is the number of angles used for surface field averaging in postprocessing.

The following 'calculation' parameters are for describing the incident light and the mediums that the light passes through.
- `epsilon1` is a floating-point value representing the dielectric constant of the surrounding medium.
- `epsilon2` is a floating-point value representing the dielectric constant of the scattering medium (the particle).
- `Raman epsilon2` is a floating-point value representing the dielectric constant of the scattering medium (the particle) at the Raman wavelength (i.e. `Raman lambda`).
- `lambda` is a floating-point value representing the wavelength of the incident light.
- `Raman lambda` is a floating-point value representing the wavelength of the scattered light, i.e. the Raman wavelength.

### Examples of valid numerical values

Integers must be typed using the syntax of C++'s integer literals. While negative integers are readable by the program, any integer parameter must be non-negative, so there's no reason to use them. The following examples give the same value:
```cpp
42 // decimal value
052 // leading 0 makes it an octal value
0x2a // leading 0x makes it a hexadecimal value
0X2A // 0X is also makes it hexadecimal, the case of hexadecimal digits are ignored
0b101010 // leading 0b makes it a binary value
```
Any value given will be converted into a `signed int` type.

Floating-point values must be typed using the syntax of C++'s floating-point literals or integer literals. Alternatively, fractions of two floating-point values can be used by typing two floating-point values with a `/` delimiter between them. Here are some examples:
```cpp
42
42. // 42
3.14
4.e2 // 400
1E-5 // 0.00001
0x1ffp2 // 130816, leading 0x means to read it as a sequence of hexadecimal values
// 'p' indicates a hex exponent part and is always necessary for hex floats.
0X10.1P0 // 16.0625
1/3 // 0.333333...
2E10/2E12 //0.25
```
Any value given will be converted to the types given in the `Calculation type` parameter in the `config.txt` file, or some other specified input file.

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

  - [x] Make the program more user-friendly
  - [ ] Make the code more readable
  - [ ] Run more rigorous tests
  - [x] Optimise code
  - [ ] Create a release

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
- auxPrepareIntegrals reads from binary `.dat` files instead of `.mat` files. It can only read files made using `storeGLquadrature`.
- sphCalculatePQ doesn't try to access stParams.output, since stParams is expected to be the same struct type as the one given in the specification, in which case 'output' is not a member of stParams. Such a member does exist in stOptions however, so future implementations may take stOptions as an argument type. (For calculating Raman scattering, this option is true by default.)
- The slvGetOptionsFromStruct function cannot be called on its own and is instead implemented into the stOptions constructors.
- The pstMakeStructForField function currently puts stIncPar into returning struct stRes by using std::move(); this means that the stIncPar will be made empty after pstMakeStructForField is used. This is done, since in the final program, the input stIncPar doesn't need to be kept.
- sphEstimateDelta doesn't get called in slvForT in the original Raman scattering MATLAB program, so this function isn't implemented.
- rvhGetSymmetricMat is not debugged since it's never used in the final `raman_elastic_scattering` program.
- rvhTruncateMatrices, rvhGetSymmetricMat, rvhGetFieldCoefficients, rvhGetAverageCrossSections can only take stTR vectors as arguments and not stPR vectors. (Overloading these functions with versions that can take stPR vectors is relatively simple however.)

## Dependencies

- A compiler with C++17 support (tested with gcc 9.3.0-17 and clang 10.0.0-4 targeting a x86_64 Ubuntu-based linux PC)
- Eigen 3.4.0, a free, open-source, efficient and comprehensive linear algebra library made for C++. Website: https://eigen.tuxfamily.org/index.php
- Boost (C++ Libraries) 1.77.0, a free, open-source set of libraries with applications in a wide variety of areas. Used for some template maths functions and its MPFR class wrapper. Website: https://www.boost.org
- The GNU Multi Precision Arithmetic Library (GMP), a C library that provides support for arbitrary precision arithmetic. While it has arbitrary-precision floating-point types, MPFR is preferred. This library is a prerequisite for MPFR. Website: https://gmplib.org/
- The GNU Multiple Precision Floating-Point Reliable Library (GNU MPFR), a C library that provides support for arbitrary-precision floating-point computation. Website: https://www.mpfr.org

## Notes

- Although most of the code is based on the SMARTIES v1.01 MATLAB package, only the functions necessary for calculating Raman scattering will be converted.
- The coding style is loosely based on Google's C++ style guide, with variable, struct and function names matching closely with the ones given in the original SMARTIES code.
- The code was written to be compiled with GCC and on hardware that correctly implements the IEEE 754 standard.

## Copyright Disclaimer

The following libraries are open-source and were written by their respective developers under their respective licenses. Neither S. Li (myself), nor T. Plakhotnik wrote any of the following code.

SMARTIES was written by Walter Somerville, Baptiste Auguié, and Eric Le Ru (copyright 2015). The package is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License. You may find the original SMARTIES package here:
https://www.wgtn.ac.nz/scps/research/research-groups/raman-lab/numerical-tools/smarties

The casting code in misc.hpp is based on code that was written by DavidAce on the following stackoverflow thread: https://stackoverflow.com/questions/48795789/eigen-unsupported-tensor-to-eigen-matrix (accessed 17 Dec 2021)

Parts of the Eigen 3.4.0 library have been included under the MPL2 license and parts of the Boost 1.77.0 library have been included under the Boost license. Eigen and Boost were made by their respective developers listed on their websites.
