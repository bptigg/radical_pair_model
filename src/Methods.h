#pragma once

#include "QuantumObject.h"

#define _USE_MATH_DEFINES
#include <math.h>


#include <array>
#define MATRIX3x3 std::array<std::array<double,3>,3>

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> commutator(MATRIX& A, MATRIX& B);

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeHamiltonian(std::vector<int32_t> dims, int ind, std::array<double, 3> parvec );
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeHamiltonian(std::vector<int32_t> dims, int ind_1, int ind_2, std::array<std::array<double, 3>, 3> parmat);

MATRIX3x3 PointDipoleDipoleCoupling(double r);
MATRIX3x3 PointDipoleDipoleCoupling(std::array<double, 3> r);