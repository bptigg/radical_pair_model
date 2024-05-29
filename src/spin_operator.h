#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

#include <complex>
#include <iostream>
#include <numeric>

#define MATRIX Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>

struct specs
{
	int ind;
	char opstr;
};

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> jmat(double j, char* args);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> _jplus(double j);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> _jz(double j);

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> spin_jx(double j);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> spin_jy(double j);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> spin_jz(double j);

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> identity(std::vector<int> dims);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> identity(int dims);

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> get_hermition(Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>& matrix);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> OperatorFunctionMapping(char args, int dims);

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeSpinOperator(std::vector<int32_t> dims, std::vector<specs> spec);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeZeroOperator(std::vector<int32_t> dims);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeTensor(std::vector<Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>>& MatrixVec, std::vector<int32_t> dims);
