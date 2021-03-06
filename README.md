# Raman-Scattering-Code-Conversion
C++ code written as part of the UQ SMP project "Efficient code for calculation of Raman scattering by spheroids"

Author is Siwan Li

Project supervisor is A/Prof Taras Plakhotnik.

## Introduction

The original MATLAB code was written to numerically calculate Raman scattering by spheroids with arbitrary precision. The original code (shared with me by T. Plakhotnik) uses the SMARTIES v1.01 MATLAB package, with some modifications made to support arbitrary precision using the Advanpix Multiprecision Toolbox. SMARTIES is a code suite used to calculate the optical properties of spheroidal particles. The goal of the project is to rewrite the code that was used for the calculation of Raman scattering in a more efficient compiled programming language that can support numerical computation. Hopefully this will accelerate the computation speeds by orders of magnitude and allow calculations for larger spheroids using higher precision.

Raman scattering is the inelastic scattering of light, i.e. a form of scattering where incident photons exchanges energy with the incident particle before propagating off in a different direction. Typically, this exchange of energy means that the particle gains or loses heat energy, and the incident photons change in frequency (usually heat energy is gained and frequency is lowered). Here, monochromatic light is scattering off an oblate spheroidal particle that behaves like a raindrop in the air. The program then calculates the amount of light scattered into some new frequency. Although the code is originally intended to model raindrop-like particles, minimal work should be needed to use the code for calculating the Raman scattering off of any spheroid (oblate or prolate) in any medium. Either change around the parameter settings shown below, or edit the `LoadParams` function in `raman_elastic_scattering.hpp`.

## Instructions

### Setup

For running the code using just single, double or quadruple-precision floating points, getting the basic, portable version from the list of releases is enough. The binaries can directly be executed as is.

If your computer doesn't use x86_64 architecture, or if you want to use floating points of arbitrary precision, you must get the source code yourself, either by cloning this repo or getting the advanced version from the list of releases, and then build the code. The code uses the GMP and MPFR libraries for doing multi-precision. Arbitrary-precision functionality must to be installed manually, since these two libraries have to be configured and compiled to the machine to fully work, and changing the amount of precision requires recompiling the Raman code. To install them, download both libraries using the links in [the dependencies section](#dependencies), extract the files and follow the instructions given in their respective 'INSTALL' files.

Open the terminal to the Raman-Scattering-Code-Conversion directory. Then simply run the following command to compile the `raman_elastic_scattering` binaries:
```bash
make [OPTIONS]
```

`[OPTIONS]` can contain any of the following parameters:
- `PRECISION=<value>`: `raman_elastic_scattering` will customise the custom arbitrary type to use `<value>` many bits. The precision allowed is only limited by the amount of memory available on your computer. If this option isn't specified, the binary will compile with a precision of 113 bits by default.

To change the number of threads or bits of precision used, `raman_elastic_scattering` must be recompiled. This as simple as running

```
make clean-mp
```

to remove the existing multiprecision binaries, then running `make` again, specifying the amount of precision with the `PRECISION` flag.

### Running `raman_elastic_scattering`

In either case, once the binaries have been built (this may take a few minutes), execute `raman_elastic_scattering` via the command:
```bash
raman_elastic_scattering [OPTIONS]
```

`[OPTIONS]` can contain any of the following parameters:
- `--input=</path/to/file>` means that `raman_elastic_scattering` reads the parameters from the file specified by `<path/to/file>`. More information about these parameters can be found in the [Configuration subsection](#configuration). By default, the input file is `config.txt`.
- `--output-dir=</path/to/directory` means that `raman_elastic_scattering` writes the program output to the directory specified by `<path/to/drectory>`, in parallel with printing them to the terminal. Each top-level thread prints to its own output file. By default, this directory is `output`.

### Configuration

You can change the calculation parameters by editing the `config.txt` file, or the file specified under the flag `--input=</path/to/file>`. Parameters must be entered in the format `Parameter:<value>`, where `<value>` is either a number or a fully lower-case word, no spaces. If no value is given, `raman_elastic_scattering` uses the default value. Otherwise if invalid values are given, it will fall back to the default values in the best case, or interpret the value in unexpected ways in the worst case. Examples of allowable formats for numerical values are given in [its own subsection](#examples-of-valid-numerical-values).

Here's what each 'run' parameter does:
- `No. of CPUs` is a non-negative integer representing the number of CPUs (threads) to be used to run the calculations. Some algorithms like matrix products can take advantage of the extra threads to parallelise calculations. If `<value>` is 0, `raman_elastic_scattering` will run with the maximum number of CPUs. By default, this parameter is 1.
- `No. of CPUs to partition particle radii for` is a non-negative integer. The different radii of the particles is partitioned among ``<value>`` threads, so each thread performs all calculations for particles of its allocated radii to its own output file. Note that each thread that was allocated a set of radii can fork into more threads to compute for the same set of radii. By default, this parameter is 1.
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
- `No. of particle thetas` is a positive integer representing the number of different `theta` to calculate with, with `theta` between 0 radians and ????/2 radians inclusive. The values of `theta` are regularly spaced. By default, this parameter is 19.
- `Particle phi` is a floating-point value representing the value of `phi` in radians. By default, this parameter is 0.
- `No. of h ratios` is a positive integer representing the number of `h` values between `Minimum h` and `Maximum h` inclusive to calculate with. These values are linearly spaced. By default, this parameter is 1.

The following parameters are for determining the spherical coordinates of the sample of points inside the particles to perform field calculations with. The physics definition of spherical coordinates are used, so `r` is the radial distance from the centre of the particle, `phi` is the azimuthal/longitudinal angle and `theta` is the polar/colatitudinal angle.
- `No. of r-coordinates` is a positive integer representing the number of radial coordinates to use within the particle when performing field calculations/integrations. By default, this parameter is 100.
- `No. of phi-coordinates` is a positive integer representing the number of azimuthal angle coordinates to use within the particle when performing field calculations/integrations. By default, this parameter is 320.
- `No. of theta-coordinates` is a positive integer representing the number of polar angle/colatitude coordinates to use within the particle when performing field calculations/integrations. By default, this parameter is 320.

The following parameters are used for calculating the T-matrices.
- `Nb_theta` is a positive integer, which is the number of angles used in Gaussian quadratures for the evaluation of P and Q matrix integrals.
- `Nb_theta_pst` is a positive integer, which is the number of angles used for surface field averaging in postprocessing.
- `Delta` is the number of extra multipoles to use for P and Q matrices. If `Delta = -1`, then the code tries to estimate it itself by calling `sphEstimateDelta`.

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

### Running `StoreGLquadrature` script

Run the command `make utils` to compile the binary, then execute the file `store_GL_quadrature`. Running this will calculate and store Gauss-Legendre quadrature Eigen::Arrays for some values of `N_theta` as raw binary files in the `data/` directory. These will be calculated for the one or two calculation types given in `config.txt`. This allows these quadratures to just be looked up in the future, instead of being calculated from scratch. Unfortunately, since the data is stored as raw binary, this makes it rather unsafe to use. Therefore, it's highly recommended to re-calculate this data when using the code on a different machine. There are also problems reading quadratures that were calculated for the custom type (i.e. `RamanFloat`), as such, storing quadratures of that type is disabled. More details are given in the source code in `smarties_aux.hpp`.

## Reading the code

`lib` contains some third-party libraries; all original code is in the `src` folder. The original functions from the SMARTIES MATLAB package are listed in files by their prefix. These SMARTIES files include the `.hpp` files `smarties_aux`, `vsh`, `sph`, `rvh`, `slv` and `pst`. However, `sphEstimateDelta` is located in `rvh.hpp`, as otherwise, `sph.hpp` and `rvh.hpp` would depend on each other. Common typedefs, structs and functions are stored in the `core.hpp`, `core_mp.hpp`, `misc.hpp` and `misc.cpp` files, so it's worth reading those ones first. Template functions are all instantiated within `defs_*.cpp` files. The functions specifically for Raman scattering is in the `raman_elastic_scattering` files, and the main functions are located in the files `main.cpp`, `main_mp.cpp` and `utils.cpp`.

Much of the code makes reference to equations from [`Mishchenko 2002`][5], [`JQSRT 2013`][3] or [`JQSRT 2015`][4], which can be found in [the references section](#references). A smaller user-guide for a more concise explanation of the theories and terminologies used can be found as the file ["SMARTIES User Guide.pdf"][2] with the original SMARTIES code. Most functions also return structs. The definition of these structs are given at the top of the file (sometimes of a different file), so more details for each struct can usually be found there. Functions or structs may also return or contain mathematical functions; this is done so by calculating the result of the function for many parameter values and tabulating these results in an Eigen::Array or Eigen::Tensor. The exact number of dimensions used depends on how many parameters there are to the function, but usually there are just 1 or 2 parameters for Arrays and sometimes 3 for Tensors.

Some functions (especially the ones relating to the individual terms and coefficients of the multipole expansion of the electric fields) have parameters `n` and `m`, where `n` goes from 1 to `N` for some specified value of `N`, and `|m| <= n`. The two parameters are usually compressed down to one parameter, using what's called the 'p-index', where `p = n*(n+1) + m`. This also reduces the amount of dimensions needed to represent the function in an Eigen::Array or Eigen::Tensor. Unlike the original MATLAB code, and in this C++ code however, the case where `(n,m) = (0,0)` is included as padding. This case gives a p-index of 0 and makes it so that `P` (the number of possible p-indices) is equal to `(N+1)^2` instead of `N*(N+2)`. This case is included in C++ so that in the arrays/tensors representing these functions, the indices `0, 1, 2, ..., N` actually correspond to the value `n = 0, 1, 2, ..., N` and the indices `0, 1, 2, ..., P-1` actually corresponds to the indices the value `p = 0, 1, 2, ..., P-1`. In the original MATLAB code, this case wasn't included, since MATLAB uses index-from-1 instead of index-from-0.

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
  - [x] Make the code more readable
  - [x] Run more rigorous tests
  - [x] Optimise code
  - [x] Create a release

## Future directions of the project

If this project were to continue, here are a few ideas for ways to improve or speed up the code, or do something else related:

- In the source code `sphMakeGeometry`, half of the terms that are calculated in some parts are always zero. This wasn't optimised in the orginal MATLAB code, since the way MATLAB ran the script made it so that it would take roughly equal time to calculate anyway. In C++, optimising this part of the code may or may not give some slight speed improvements (it's possible that these zero cases are optimised out during compile time anyway).
- Implement multi-processing and/or GPU acceleration at a low-level and throughout the entire code. Since much of the code involves doing computations with a lot of data, having them done simultaneously via multi-processing and/or GPU acceleration may significantly improve the computation speeds. NVIDIA CUDA already has a C++ library as well and Eigen has some support for CUDA.
- Analytically reduce/simplify the maths required to get the desired results. Spheroids have many rotation and reflection symmetries that could still be exploited to help simplify the maths. This would especially be helpful for the field calculation parts, where the triple integral could be reduced down to a double integral. Out of all methods, this could give the greatest improvement for the least amount of effort.
- Convert the whole SMARTIES package to a C++ or FORTRAN library. Of course, this would mean that a lot of restructuring would need to be done to ensure cohesiveness of code, and many functions would have to be re-analysed to make them work well with C++/FORTRAN.

## Notable differences from SMARTIES

- Expansion coefficients and other arrays that store using p-indices now include values for when n=0, m=0. This makes P (the length of the p-vectors) equal to (N+1)^2 and makes the code more convenient in a index-by-zero language without compromising any calculations. As a side effect, many other functions that calculate values for various values n are affected as well. It is not recommended to use these values, as they can often be undefined or not calculated. Functions that are affected by this include the following. (Note that in almost all cases, the extra values are 0, so they act as padding.)
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
- Some structs that originally contained [1 x X] matrix members may now contain [X x 1] matrix members instead for better syntax and brevity of code. Structs that have this include (but are not limited to):
  - stIncEabnm
  - stEforPhi
  - stAbcdnm
- auxPrepareIntegrals reads from binary `.dat` files instead of `.mat` files. It can only read files made using `storeGLquadrature`.
- The slvGetOptionsFromStruct function cannot be called on its own and is instead implemented into the stOptions constructors.
- The pstMakeStructForField function currently puts stIncPar into returning struct stRes by using std::move(); this means that the stIncPar will be made empty after pstMakeStructForField is used. This is done, since in the final program, the input stIncPar doesn't need to be kept.
- rvhGetSymmetricMat, rvhGetFieldCoefficients, rvhGetAverageCrossSections can only take stTR vectors as arguments and not stPR vectors. (Overloading these functions with versions that can take stPR vectors is relatively simple however.)
- Many of the original functions had extra functionality for when the code is ran on Octave instead of MATLAB. Since this is only C++, such code (including the `isOctave` function itself) is not implemented.
- sphEstimateDelta may give different values of Delta between the MATLAB version and the C++ version. This is because sphEstimateNB is very sensitive to rounding errors, and exactly how and when floating-points are rounded differ between different implementations even in the same language. Although this shouldn't matter in the expected use cases and in this program, where Delta is just added on to some maximum multipole order `N` to get an even larger multipole order to compute the P, Q, T and R matrices with, before those same matrices get truncated down to a maximum multipole of order `N` anyway. If it absolutely matters, manually specifying a value of Delta is easy to do anyway.

## Dependencies

- The basic, portable executable requires glibc version 2.17 or above to work.
- Building from source code requires a compiler with proper C++17 support (tested with gcc 9.3.0-17 and clang 10.0.0-4 targeting a x86_64 Ubuntu-based gnu-linux PC)
- [Eigen 3.4.0][https://eigen.tuxfamily.org/index.php], a free, open-source, efficient and comprehensive linear algebra library made for C++. The necessary files from this library are already included with the source code.
- [Boost (C++ Libraries) 1.77.0][https://www.boost.org], a free, open-source set of libraries with applications in a wide variety of areas. Used for some template maths functions and its MPFR class wrapper. The necessary files from this library are already included with the source code.
- [The GNU Multi Precision Arithmetic Library (GMP)][https://gmplib.org/], a C library that provides support for arbitrary precision arithmetic. While it has arbitrary-precision floating-point types, MPFR is preferred. This library is a prerequisite for MPFR.
- [The GNU Multiple Precision Floating-Point Reliable Library (GNU MPFR)][https://www.mpfr.org], a C library that provides support for arbitrary-precision floating-point computation.

## Notes

- Although most of the code is based on the SMARTIES v1.01 MATLAB package, only the functions necessary for calculating Raman scattering will be converted.
- The coding style is loosely based on Google's C++ style guide, with variable, struct and function names matching closely with the ones given in the original SMARTIES code.
- The code was written to be compiled with GCC and on hardware that correctly implements the IEEE 754 standard.

## Copyright Disclaimers and Acknowledgements

I acknowledge the support of the School of Maths and Physics of the University of Queensland in organising, funding and providing equipment for this research project.

The following libraries are open-source and were written by their respective developers under their respective licenses. Neither S. Li (myself), nor T. Plakhotnik wrote any of the following code.

SMARTIES was written by Walter Somerville, Baptiste Augui??, and Eric Le Ru (copyright 2015). The package is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License. You may find the original SMARTIES package [here][1]

SMARTIES itself contains converted/borrowed Fortran code originally written by Mishchenko, Travis and Mackowski ([accessible here][6]) and MATLAB code originally written by von Winckel ([accessible here][7]). Von Winckel's code is covered by the BSD license.

The casting code in `misc.hpp` is based on code that was written by DavidAce on stackoverflow [in this thread][8] (accessed 17 Dec 2021). The binary I/O code in `smarties_aux.hpp` is based on code that was written by andrea also on stackoverflow [in this thread][9] (accessed 24 Jan 2022).

Parts of the Eigen 3.4.0 library have been included under the MPL2 license and parts of the Boost 1.77.0 library have been included under the Boost license. Eigen and Boost were made by their respective developers listed on their websites.

## References

1. SMARTIES [Internet]. Wellington: School of Chemical and Physical Sciences, Victoria University of Wellington; May 2016. Available from: https://www.wgtn.ac.nz/scps/research/research-groups/raman-lab/numerical-tools/smarties
2. Somerville WRC, Augui?? B, Le Ru EC. 2016 May. SMARTIES: User-friendly codes for fast and accurate calculations of light scattering by spheroids. Journal of Quantitative Spectroscopy and Radiative Transfer 174: 39-55
3. Somerville WRC, Augui?? B, Le Ru EC. 2013 July. A new numerically stable implementation of the T-matrix method for electromagnetic scattering by spheroidal particles. Journal of Quantitative Spectroscopy and Radiative Transfer 123: 153-168
4. Somerville WRC, Augui?? B, Le Ru EC. 2015 July. Accurate and convergent T-matrix calculations of light scattering by spheroids. Journal of Quantitative Spectroscopy and Radiative Transfer 160: 29-35
5. Mishchenko MI, Travis LD, Lacis AA. 2002. Scattering, absorption, and emission of light by small particles. Cambridge: Cambridge University Press. 445 p.
6.  Mishchenko MI, Travis LD, Mackowski DW. T-Matrix Codes for Computing Electromagnetic Scattering by Nonspherical and Aggregated Particles [Internet]. New York City: NASA Goddard Institute for Space Studies; 2020 July. Available from: https://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html
7. Von Winckel G. Legendre-Gauss Quadrature Weights and Nodes [Internet]. MATLAB Central File Exchange. 2004 May 11. Available from: https://au.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes
8. Eigen unsupported Tensor to Eigen matrix [Internet]. Stack Overflow; 2021 Dec 17. Available from: https://stackoverflow.com/questions/48795789/eigen-unsupported-tensor-to-eigen-matrix
9. How to write/read an Eigen matrix from binary file [Internet]. Stack Overflow; 2022 Jan 24. Available from: https://stackoverflow.com/questions/25389480/how-to-write-read-an-eigen-matrix-from-binary-file

## External Links
[1]: https://www.wgtn.ac.nz/scps/research/research-groups/raman-lab/numerical-tools/smarties "SMARTIES page"

[2]: https://doi.org/10.1016/j.jqsrt.2016.01.005 "SMARTIES user guide"

[3]: https://www.sciencedirect.com/science/article/pii/S0022407313000423 "JQSRT 2013"

[4]: https://www.sciencedirect.com/science/article/pii/S0022407315001120 "JQSRT 2015"

[5]: https://scholar.google.com/scholar_lookup?title=Scattering%2C%20absorption%20and%20emission%20of%20light%20by%20small%20particles&author=M.I.%20Mishchenko&publication_year=2002 "Mishchenko 2002"

[6]: https://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html "Mishchenko Fortran"

[7]: https://au.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes "lgwt.m"

[8]: https://stackoverflow.com/questions/48795789/eigen-unsupported-tensor-to-eigen-matrix "Eigen tensor/matrix casting"

[9]: https://stackoverflow.com/questions/25389480/how-to-write-read-an-eigen-matrix-from-binary-file "Eigen matrix I/O"
