#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor>

#include <complex>
#include <iostream>

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> jmat(double j, char* args);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> _jplus(double j);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> _jz(double j);

Eigen::Tensor<double, Eigen::RowMajor> identity(std::vector<int> dims);
Eigen::SparseMatrix<double, Eigen::RowMajor> identity(int dims);

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> get_hermition(Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>& matrix);