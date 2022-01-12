#include "core.h"

using namespace std;
using namespace Eigen;
using namespace Smarties;

template <class Real>
unique_ptr<stParams<Real>> loadParam(string type = "");

template <class Real>
void RamanElasticScattering(int argc, char** argv);
