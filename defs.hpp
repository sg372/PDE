
#ifndef DEFS_HPP_
#define DEFS_HPP_

#include "Eigen/Dense"


using namespace Eigen;

#define complex std::complex<double>

//Definition of the matrix type using eigen
typedef Matrix<unsigned, Dynamic, Dynamic> UnsignedMatrix;

typedef Matrix<double, Dynamic, Dynamic> DoubleMatrix;

typedef Matrix<complex, Dynamic, Dynamic> ComplexMatrix;

typedef Matrix<double, Dynamic, 1> DoubleVector;

#endif   /* DEFS_HPP_ */

