#pragma once

#include "QuantumObject.h"
#include "ThreadPool.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <array>
#include <Eigen/Dense>
#include <thread>

#define MATRIX3x3 std::array<std::array<double,3>,3>

#define MatrixCol Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor>

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> commutator(Matrix& A, Matrix& B);

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeHamiltonian(std::vector<int32_t> dims, int ind, std::array<double, 3> parvec);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeHamiltonian(std::vector<int32_t> dims, int ind_1, int ind_2, std::array<std::array<double, 3>, 3> parmat);
Eigen::VectorXcd FlattenMatrix(Matrix mat);
std::vector<std::complex<double>> FlattenMatrixVec(Matrix mat);
Matrix ReformMatrix(std::vector<std::complex<double>> vec);

Eigen::VectorXcd DotProduct(Matrix& mat, const Eigen::VectorXcd& vec);
//std::vector<std::complex<double>> DotProductVec(Eigen::MatrixXcd& mat, const std::vector<std::complex<double>>& vec);
std::vector<std::complex<double>> DotProductVec(Matrix& mat, const std::vector<std::complex<double>>& vec = {}, bool term = false);

void StartThreadPool(int max_threads);
void KillThreadPool();


MATRIX3x3 PointDipoleDipoleCoupling(double r);
MATRIX3x3 PointDipoleDipoleCoupling(std::array<double, 3> r);

Eigen::Vector4d SingletState();

double simpson_integration(std::vector<double> x_list, std::vector<double> y_list = {});
std::vector<std::array<double,3>> FibonacciSphere(int n);

Matrix BlockInverse(Matrix mat, int dim, int inner_block_size);

template<typename t, int s1, int s2>
Eigen::Vector<t, s1* s2> TensorProduct(Eigen::Vector<t, s1> v1, Eigen::Vector<t, s2> v2)
{
	Eigen::Vector<t, s1* s2> ReturnVec;

	for (int i = 0; i < s1; i++)
	{
		int offset = i * s2;
		for (int e = 0; e < s2; e++)
		{
			ReturnVec[offset + e] = v1[i] * v2[e];
		}
	}

	return ReturnVec;
}




