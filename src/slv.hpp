/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This code contains all 'slv' SMARTIES functions that are used in Raman scattering calculations,
i.e. high level functions for solving a specific class of problems.
Note that although the original MATLAB code used the function slvGetOptionsFromStruct,
no such function is defined here, and the functionality is instead defined using
the default constructor of stOptions.
*/

#ifndef SLV_HPP
#define SLV_HPP

#include "core.hpp"
#include "rvh.hpp"
#include "vsh.hpp"
#include "sph.hpp"

using namespace Eigen;
using namespace std;

namespace Smarties {

  /* Struct containing a T-matrix with orientation-averaged properties */
  template <class Real>
  struct stTmatrix {
    unique_ptr<stCrossSection<Real>> st_C_oa;
    vector<unique_ptr<stTR<Real>>> st_TR_list;
  };

  /*
  Calculates the T-matrix with orientation-averaged properties.
  Inputs:
    params - a unique pointer to a stParams struct.
             Here, only the `a`, `c`, `k1`, `s`, `N` and `Nb_theta` member variables are required.
             More details of the struct can be found in core.hpp
    options - a unique pointer to a stOptions struct.
  Output:
    Unique pointer to a new stTmatrix struct, containing member variables
    st_C_oa and st_TR_list with size M = abs_m_vec.size().
  Dependencies:
    rvhGetAverageCrossSections, rvhGetSymmetricMat, rvhGetTRfromPQ,
    rvhTruncateMatrices, sphCalculatePQ, sphEstimateNB, sphMakeGeometry, sphEstimateDelta
  */
  template <class Real>
  unique_ptr<stTmatrix<Real>> slvForT(const unique_ptr<stParams<Real>>& params,
      const unique_ptr<stOptions>& options,
      unique_ptr<stRtfunc<Real>> st_geometry = unique_ptr<stRtfunc<Real>>()) {
    Real c = params->c;
    Real a = params->a;

    ArrayXr<Real> s = params->s;
    ArrayXr<Real> k1 = params->k1;

    auto st_k1_s = make_unique<stParams<Real>>();
    st_k1_s->k1 = k1;
    st_k1_s->s = s;

    int N = params->N;
    int Nb_theta = params->Nb_theta;
    // If `options` member `abs_m_vec` is non-zero in size, i.e. is defined, then use member `abs_m_vec`.
    // Otherwise construct a new Eigen::Array.
    ArrayXi abs_m_vec = options->abs_m_vec.size() ? options->abs_m_vec : ArrayXi::LinSpaced(N + 1, 0, N);

    st_k1_s->output = options->output;

    // Make structure describing spheroidal geometry and quadrature points for numerical integrations
    if (!st_geometry)
      st_geometry = sphMakeGeometry(Nb_theta, a, c);

    if (options->delta < 0) { // then need to estimate delta
      stDelta<Real> st_delta = sphEstimateDelta(st_geometry, st_k1_s);
      if (st_delta.delta == -1) {
        throw runtime_error(string("ERROR: Delta could not be found. Results are likely to be non-converged. ") +
            string("Try choosing Delta manually instead."));
      }
      options->delta = st_delta.delta;
      cout << "Delta estimated to Delta = " << st_delta.delta <<
          " with relative error in T_{11}^{22,m=1} of " << st_delta.err << endl;
    }

    int NQ = N + options->delta; // NQ>=N: Maximum multipole order for computing P and Q matrices
    int NB = options->NB;

    // Estimating NB, the number of multipoles to compute the Bessel functions (NB >= NQ)
    if (NB <= 0)
      NB = sphEstimateNB(NQ, st_geometry, st_k1_s);
    if (NB < NQ)
      NB = NQ; // NB must be at least NQ

    // Calculate P and Q matrices
    vector<unique_ptr<stPQ<Real>>> st_PQ_list = sphCalculatePQ(NQ, abs_m_vec, st_geometry, st_k1_s, NB);
    // Calculate T (and possibly R) matrices
    vector<unique_ptr<stTR<Real>>> st_TR_list = rvhGetTRfromPQ(st_PQ_list, options->get_R);

    // If needed, discard higher order multipoles (which are affected by the finite size of P and Q)
    if (NQ > N)
      st_TR_list = rvhTruncateMatrices(st_TR_list, N);
    // T and R matrices now go up to N multipoles

    // If required, symmetrize the T-matrix
    if (options->get_symmetric_T)
      st_TR_list = rvhGetSymmetricMat(st_TR_list);

    // Calculate the (Ext, Abs, Sca) cross-sections for orientation-averaged excitation
    // st_TR_list converted into a two-dimensional std::vector for compatibility with rvhGetAverageCrossSections
    vector<vector<unique_ptr<stTR<Real>>>> st_TR_array;
    st_TR_array.push_back(move(st_TR_list));
    unique_ptr<stCrossSection<Real>> st_coa = rvhGetAverageCrossSections(params->k1, st_TR_array);
    st_TR_list = move(st_TR_array[0]);

    auto output = make_unique<stTmatrix<Real>>();
    output->st_C_oa = move(st_coa);
    output->st_TR_list = move(st_TR_list);

    return output;
  }
}

#endif
