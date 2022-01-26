#ifndef SMARTIES_AUX_HPP
#define SMARTIES_AUX_HPP

#include "core.hpp"
#include "math.hpp"
#include <filesystem>
#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;

namespace Smarties {

  enum sInt {GAUSS, RECTANGLE, PTS};

  template <class Real>
  struct stGLQuad {
    ArrayXr<Real> x;
    ArrayXr<Real> w;
  };

  template <class Real>
  struct stRtfunc {
    int Nb_theta;
    ArrayXr<Real> theta;
    ArrayXr<Real> w_theta;
    ArrayXr<Real> r;
    ArrayXr<Real> dr_dt;
    Real a;
    Real c;
    Real h;
    Real r0;
    sInt type;
  };

  template <class Real>
  unique_ptr<stGLQuad<Real>> auxInitLegendreQuad(int N, Real a = -1.0, Real b = 1.0) {
    assert(N >= 1);
    N--;
    int N1 = N + 1;
    int N2 = N + 2;

    ArrayXr<Real> xu = ArrayXr<Real>::LinSpaced(N1, -1, 1);
    ArrayXr<Real> y = cos((2*ArrayXr<Real>::LinSpaced(N1, 0, N) + 1)*mp_pi<Real>()/(2*N + 2)) +
        (0.27/N1)*sin(mp_pi<Real>()*xu*N/N2);
    ArrayXr<Real> y0(N1);
    y0.fill(2);
    ArrayXXr<Real> L(N1, N2);
    L.col(0).setOnes();
    ArrayXr<Real> Lp(N1);

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

    output->x = (a*(1 - y) + b*(1 + y))/2;
    output->w = ArrayXr<Real>::Constant(N1, b - a) /
        ((1 - y.pow(2))*Lp.pow(2)) * pow(static_cast<Real>(N2) / N1, 2);
    return output;
  }

  template<class Matrix>
  void writeBinary(const string filename, const Matrix& matrix) {
    ofstream out(filename, ios::out | ios::binary | ios::trunc);
    typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
    out.write((char*) (&rows), sizeof(typename Matrix::Index));
    out.write((char*) (&cols), sizeof(typename Matrix::Index));
    out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar));
    out.close();
  }

  template<class Matrix>
  void readBinary(const string filename, Matrix& matrix) {
    ifstream in(filename, ios::in | ios::binary);
    typename Matrix::Index rows=0, cols=0;
    in.read((char*) (&rows),sizeof(typename Matrix::Index));
    in.read((char*) (&cols),sizeof(typename Matrix::Index));
    matrix.resize(rows, cols);
    in.read( (char *) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar));
    in.close();
  }

  template <class Real>
  string getTypeName() {
    string calc_type = typeid(static_cast<Real>(0)).name();
    if (calc_type.find("N5boost14multiprecision6numberINS0_8backends18mpfr_float_backend") != string::npos) {
      size_t prec_begin = calc_type.find("ILj") + 3;
      size_t prec_end = calc_type.find("ELNS");
      string precision = calc_type.substr(prec_begin, prec_end - prec_begin);
      calc_type = "custom_" + precision;
    } else if (calc_type == "f")
      calc_type = "single";
    else if (calc_type == "d")
      calc_type = "double";
    else if (calc_type == "e")
      calc_type = "quad";
    return calc_type;
  }

  template <class Real>
  void updateGLquadrature(ArrayXi N_theta, bool keep_old = true) {
    string calc_type = getTypeName<Real>();
    string data_path = "data/" + calc_type;
    string quad_table_file_name = data_path + "/quad_table.dat";
    std::filesystem::create_directory(data_path);
    ifstream quad_table_file(quad_table_file_name, ios::in | ios::binary);

    ArrayXi new_N_theta(N_theta.size());
    if (keep_old && quad_table_file.good()) {
      ArrayXi old_N_theta;
      readBinary(quad_table_file_name, old_N_theta);
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
      writeBinary(values_file_name, values);
    }
    quad_table_file.close();

    writeBinary(quad_table_file_name, N_theta);
  }

  template <class Real>
  void storeGLquadrature() {
    ArrayXi N_theta1 = seq2Array(50, 505, 5);
    ArrayXXi N_theta2((2000 - 600) / 100 + 1, 2);
    N_theta2.col(0) = seq2Array(600, 2000, 100);
    N_theta2.col(1) = N_theta2.col(0) + 5;
    ArrayXi N_theta3 = N_theta2.reshaped();
    ArrayXi N_theta(N_theta1.size() + N_theta3.size());
    N_theta << N_theta1, N_theta3;

    updateGLquadrature<Real>(N_theta, false);
  }

  template <class Real>
  unique_ptr<stRtfunc<Real>> auxPrepareIntegrals(int N_int, sInt type) {
    auto output = make_unique<stRtfunc<Real>>();
    output->Nb_theta = N_int;

    switch (type) {
      case GAUSS: {
        if (N_int < 50) {
          unique_ptr<stGLQuad<Real>> GL_quad = auxInitLegendreQuad<Real>(N_int);
          output->theta = acos(GL_quad->x);
          output->w_theta = GL_quad->w;
        } else {
          string calc_type = getTypeName<Real>();
          string data_path = "data/" + calc_type;
          string values_file_name = data_path + "/quad_table_values_" + to_string(N_int) + ".dat";
          ifstream values_file(values_file_name, ios::in | ios::binary);
          if (values_file.good()) {
            ArrayXXr<Real> values;
            readBinary(values_file_name, values);
            output->theta = values.col(0);
            output->w_theta = values.col(1);
          } else {
            unique_ptr<stGLQuad<Real>> GL_quad = auxInitLegendreQuad<Real>(N_int);
            output->theta = acos(GL_quad->x);
            output->w_theta = GL_quad->w;
          }
        }
        break;
      }

      case RECTANGLE: {
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
