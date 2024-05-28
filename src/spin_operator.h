#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

#include <complex>
#include <iostream>
#include <numeric>

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> jmat(double j, char* args);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> _jplus(double j);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> _jz(double j);

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> spin_jx(double j);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> spin_jy(double j);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> spin_jz(double j);

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> identity(std::vector<int> dims);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> identity(int dims);

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> get_hermition(Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>& matrix);