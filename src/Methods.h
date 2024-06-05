#pragma once

#include "QuantumObject.h"
#include "ThreadPool.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <array>
#include <Eigen/Dense>
#include <thread>

#define MATRIX3x3 std::array<std::array<double,3>,3>

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> commutator(Matrix& A, Matrix& B);

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeHamiltonian(std::vector<int32_t> dims, int ind, std::array<double, 3> parvec);
Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeHamiltonian(std::vector<int32_t> dims, int ind_1, int ind_2, std::array<std::array<double, 3>, 3> parmat);
Eigen::VectorXcd FlattenMatrix(Matrix mat);
std::vector<std::complex<double>> FlattenMatrixVec(Matrix mat);
Matrix ReformMatrix(std::vector<std::complex<double>> vec);

Eigen::VectorXcd DotProduct(Matrix& mat, const Eigen::VectorXcd& vec);
//std::vector<std::complex<double>> DotProductVec(Eigen::MatrixXcd& mat, const std::vector<std::complex<double>>& vec);
std::vector<std::complex<double>> DotProductVec(Matrix& mat, const std::vector<std::complex<double>>& vec);

void StartThreadPool(int max_threads);
void KillThreadPool();


MATRIX3x3 PointDipoleDipoleCoupling(double r);
MATRIX3x3 PointDipoleDipoleCoupling(std::array<double, 3> r);

Eigen::Vector4d SingletState();

double simpson_integration(std::vector<double> x_list, std::vector<double> y_list = {});

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

template<typename t>
void Merge(std::vector<std::pair<int, t>>& arr, int left, int mid, int right)
{
	int l_size = mid - left + 1;
	int r_size = right - mid;

	std::vector<std::pair<int, t>> l(l_size);
	std::vector<std::pair<int, t>> r(r_size);

	for (int i = 0; i < l_size; i++)
	{
		l[i] = arr[left + i];
	}
	for (int i = 0; i < r_size; i++)
	{
		r[i] = arr[mid + 1 + i];
	}

	int LeftIndex = 0, RightIndex = 0;
	int CurrentIndex = left;

	while ((LeftIndex < l_size) && (RightIndex < r_size))
	{
		if (l[LeftIndex].first <= r[RightIndex].first)
		{
			arr[CurrentIndex] = l[LeftIndex];
			LeftIndex++;
		}
		else
		{
			arr[CurrentIndex] = r[RightIndex];
			RightIndex++;
		}
		CurrentIndex++;
	}

	while (LeftIndex < l_size) { arr[CurrentIndex++] = l[LeftIndex++]; }
	while (RightIndex < r_size) { arr[CurrentIndex++] = r[RightIndex++]; }
	
}

template<typename t>
void MergeSort(std::vector<std::pair<int, t>>& arr, int left, int right)
{
	if (left >= right)
	{
		return;
	}
	
	int mid = (left + right) / 2;
	MergeSort(arr, left, mid);
	MergeSort(arr, mid + 1, left);
	Merge(arr, left, mid, right);
}

void sort(std::vector<std::pair<int, std::complex<double>>>& arr);


