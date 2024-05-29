#pragma once

#include "spin_operator.h"

enum class type
{
	bra = 0,
	ket = 1,
	oper = 2,
	supoper = 3,
};

class QuantumObject
{
private:
	type m_type;
	Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> m_data;
	std::pair<int32_t, int32_t> m_dims; 

public:
	QuantumObject(Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>& data, type t);
	QuantumObject(Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>& data, type t, std::pair<int32_t, int32_t> dims);

	Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> get_data();
	std::pair<int32_t, int32_t> get_dims();
	type get_type();
};

