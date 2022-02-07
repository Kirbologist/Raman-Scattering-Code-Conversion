/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This code contains all 'aux' SMARTIES functions that are used in Raman scattering calculations,
as well as all the 'utils' functions related to storing data about
Gauss-Legendre quadratures in a lookup-table-like fashion.
These 'aux' functions are auxiliary functions (low-level).
Note that the file is named 'smarties_aux' instead of 'aux', since 'aux' is already a reserved
file name in some operating systems.
*/

#ifndef SMARTIES_AUX_HPP
#define SMARTIES_AUX_HPP

#include "core.hpp"
#include "misc.hpp"
#include <filesystem>
#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;

namespace Smarties {

  enum sInt {GAUSS, RECTANGLE, PTS};

  /* Struct containing information for a Gaussian Legendre quadrature of order N */
  template <class Real>
  struct stGLQuad {
    ArrayXr<Real> x;  // Legendre-Gauss nodes, [(N + 1) X 1]
    ArrayXr<Real> w; // Legendre-Gauss weights, [(N + 1) X 1]
  };

  /* Struct containing geometric information for spheroids */
  template <class Real>
  struct stRtfunc {
    int Nb_theta; // number of theta
    ArrayXr<Real> theta; // theta angles where things are computed
    ArrayXr<Real> w_theta; // weights from quadrature
    ArrayXr<Real> r; // r(theta) defining the geometry
    ArrayXr<Real> dr_dt; // derivative of r(theta)
    Real a; // semi-axis length for axes along x, y
    Real c; // semi-axis length along z
    Real h; // the aspect ratio between max(a, c) and min(a, c)
    Real r0; // equivalent-volume-sphere radius
    sInt type; // describes how the points were computed
  };

  /*
  Calculates nodes and weights for Legendre Gaussian quadrature of order N.
  The original MATLAB script was written by Greg von Winckel and distributed on MatlabCentral file exchange
  website (script name lgwt.m) and is covered by the BSD license:

  Copyright (c) 2009, Greg von Winckel
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

  This script is for computing definite integrals using Legendre-Gauss
  Quadrature. Computes the Legendre-Gauss nodes and weights on an interval
  [a,b] with truncation order N

  Suppose you have a continuous function f(x) which is defined on [a,b]
  which you can evaluate at any x in [a,b]. Simply evaluate it at all of
  the values contained in the x vector to obtain a vector f. Then compute
  the definite integral using sum(f.*w);

  Written by Greg von Winckel - 02/25/2004
  */
  template <class Real>
  unique_ptr<stGLQuad<Real>> auxInitLegendreQuad(int N, Real a = -1.0, Real b = 1.0) {
    assert(N >= 1);
    N--;
    int N1 = N + 1;
    int N2 = N + 2;

    ArrayXr<Real> xu = ArrayXr<Real>::LinSpaced(N1, -1, 1);

    // Initial guess
    ArrayXr<Real> y = cos((2*ArrayXr<Real>::LinSpaced(N1, 0, N) + 1)*mp_pi<Real>()/(2*N + 2)) +
        (0.27/N1)*sin(mp_pi<Real>()*xu*N/N2);

    ArrayXXr<Real> L(N1, N2); // Legendre-Gauss Vandermonde Matrix
    L.col(0).setOnes();
    ArrayXr<Real> Lp(N1); // Derivative of LGVM

    // Compute the zeros of the N + 1 Legendre polynomial
    // using the recursion relation and the Newton-Raphson method
    ArrayXr<Real> y0(N1);
    y0.fill(2);
    int n_iter = 0;
    while ((n_iter < 15) && mp_eps<Real>() <
        abs(acos(y) - acos(y0.template cast<complex<Real>>())).maxCoeff()) {
      n_iter++;
      L.col(1) = y;

      for (int i = 1; i < N1; i++)
        L.col(i + 1) = ((2*i + 1)*y * L.col(i) - i*L.col(i - 1))/(i + 1);

      Lp = N2*(L.col(N) - y*L.col(N1))/(1 - y.pow(2));

      y0 = y;
      y = y0 - L.col(N1)/Lp;
    }

    auto output = make_unique<stGLQuad<Real>>();

    // Linear map from [-1, 1] to [a, b]
    output->x = (a*(1 - y) + b*(1 + y))/2;
    // Compute the weights
    output->w = ArrayXr<Real>::Constant(N1, b - a) /
        ((1 - y.pow(2))*Lp.pow(2)) * pow(static_cast<Real>(N2) / N1, 2);
    return output;
  }

  /*
  Writes the raw binary data of an Eigen::MatrixBase object to a file, with a header specifying the number
  of rows and columns of the object. Raw binary data is used, as writing and reading
  the matrix entries as decimal values loses some precision,
  even if it's formatted to write with 'full precision'.
  Inputs:
    filename - path of the data file to write to
    matrix - matrix thats data is getting written
  */
  template<class Matrix>
  void WriteBinary(const string filename, const Matrix& matrix) {
    ofstream out(filename, ios::out | ios::binary | ios::trunc);
    typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
    out.write((char*) (&rows), sizeof(typename Matrix::Index));
    out.write((char*) (&cols), sizeof(typename Matrix::Index));
    out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar));
    out.close();
  }

  /*
  Reads the raw binary data of an Eigen::MatrixBase object from a file, with a header specifying the number
  of rows and columns of the Eigen::MatrixBase. The read data gets written into an Eigen::MatrixBase object.
  Raw binary data is used, as writing and reading the matrix entries as decimal values loses some precision,
  even if it's formatted to write with 'full precision'.
  Inputs:
    filename - path of the data file to read from
    matrix - matrix to overwrite the data of
  */
  template<class Matrix>
  void ReadBinary(const string filename, Matrix& matrix) {
    ifstream in(filename, ios::in | ios::binary);
    typename Matrix::Index rows=0, cols=0;
    in.read((char*) (&rows),sizeof(typename Matrix::Index));
    in.read((char*) (&cols),sizeof(typename Matrix::Index));
    matrix.resize(rows, cols);
    in.read( (char *) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar));
    in.close();
  }

  /*
  Calculates and stores points and weights for integral quadrature.
  This calculates for values of Nt not already stored in the 'data/<calc_type>/quad_table.dat' file.
  If the file doesn't exist, it creates it then it populates it with the values Nt.
  Inputs:
    Nt - vector containing the nodes to be tabulated
    keep_old - if true, add to existing data, otherwise overwrite any existing data
  Dependencies:
    auxInitLegendreQuad
  */
  template <class Real>
  void UpdateGLquadrature(ArrayXi N_theta, bool keep_old = true) {
    string calc_type = GetTypeName<Real>();
    string data_path = "data/" + calc_type;
    string quad_table_file_name = data_path + "/quad_table.dat";
    std::filesystem::create_directory(data_path);
    ifstream quad_table_file(quad_table_file_name, ios::in | ios::binary);

    ArrayXi new_N_theta(N_theta.size());
    if (keep_old && quad_table_file.good()) {
      ArrayXi old_N_theta;
      ReadBinary(quad_table_file_name, old_N_theta);
      // we only want to add in entries that aren't already there.
      auto it = set_difference(N_theta.data(), N_theta.data() + N_theta.size(),
          old_N_theta.data(), old_N_theta.data() + old_N_theta.size(), new_N_theta.data());
      new_N_theta.conservativeResize(distance(new_N_theta.data(), it));
    } else
      new_N_theta = N_theta;

    for (int i = 0; i < new_N_theta.size(); i++) {
      unique_ptr<stGLQuad<Real>> GL_quad = auxInitLegendreQuad<Real>(new_N_theta(i));
      ArrayXXr<Real> values(GL_quad->x.size(), 2);
      values << acos(GL_quad->x), GL_quad->w;
      string values_file_name = data_path + "/quad_table_values_" + to_string(new_N_theta(i)) + ".dat";
      WriteBinary(values_file_name, values);
    }
    quad_table_file.close();

    WriteBinary(quad_table_file_name, N_theta);
  }

  /*
  Calculates and stores points and weights for intgral quadrature.
  WARNING: the data gets stored as raw binary, the format of which may or may not differ from system to system.
  Thus, if the program is being used on a new system, it's highly recommended that you re-run this.
  This calculates Nt in steps of 5 from 50 to 505, then in steps of 100 from 600 to 2000,
  as well as 5 above each of those values (605, 705, etc.).
  */
  template <class Real>
  void StoreGLquadrature() {
    ArrayXi N_theta1 = Seq2Array(50, 505, 5);
    ArrayXXi N_theta2((2000 - 600) / 100 + 1, 2);
    N_theta2.col(0) = Seq2Array(600, 2000, 100);
    N_theta2.col(1) = N_theta2.col(0) + 5;
    ArrayXi N_theta3 = N_theta2.reshaped();
    ArrayXi N_theta(N_theta1.size() + N_theta3.size());
    N_theta << N_theta1, N_theta3;

    UpdateGLquadrature<Real>(N_theta, false);
  }

  /*
  Calculates points and weights for integral quadrature, using the given method.
  When using Gaussian quadrature, precalculated values stored in files in the 'data/' folder
  will be used if available. Those values include N_int in the range 50 to 500 by steps of 5,
  and from 600 to 2000 by steps of 100 with pairs 600, 605, 700, 705, ..., 2000, 2005 to facilitate convergence tests
  Inputs:
    N_int - the number of integration points required
    type - the integration scheme, either `GAUSS` for Gaussian quadrature, or `rectangle` for Simpson
  Output:
    Returns unique pointer to a stRTfunc struct
  Dependencies:
    auxInitLegendreQuad
  */
  template <class Real>
  unique_ptr<stRtfunc<Real>> auxPrepareIntegrals(int N_int, sInt type) {
    auto output = make_unique<stRtfunc<Real>>();
    output->Nb_theta = N_int;

    switch (type) {
      case GAUSS: { // For Gauss-Legendre quadrature
        if (N_int < 50) {
          unique_ptr<stGLQuad<Real>> GL_quad = auxInitLegendreQuad<Real>(N_int);
          output->theta = acos(GL_quad->x); // [T X 1]
          output->w_theta = GL_quad->w; // [T X 1]
        } else {
          string calc_type = GetTypeName<Real>();
          string data_path = "data/" + calc_type;
          string values_file_name = data_path + "/quad_table_values_" + to_string(N_int) + ".dat";
          ifstream values_file(values_file_name, ios::in | ios::binary);
          if (values_file.good()) { // if tabulated values for the given N_int exist
            ArrayXXr<Real> values;
            ReadBinary(values_file_name, values);
            output->theta = values.col(0); // [T X 1]
            output->w_theta = values.col(1); // [T X 1]
          } else {
            unique_ptr<stGLQuad<Real>> GL_quad = auxInitLegendreQuad<Real>(N_int);
            output->theta = acos(GL_quad->x); // [T X 1]
            output->w_theta = GL_quad->w; // [T X 1]
          }
        }
        break;
      }

      case RECTANGLE: {
        // For Simpson integration
        output->theta = ArrayXr<Real>::LinSpaced(N_int, 0, mp_pi<Real>());
        Real d_theta = mp_pi<Real>()/(N_int - 1);
        output->w_theta = d_theta*sin(output->theta);
        break;
      }

      default:
        cout << "Integration type not recognized" << endl;
    }

    return output;
  }
}

#endif
