#include "QuantumObject.h"

QuantumObject::QuantumObject(Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>& data, type t)
	:m_type(t)
{
	int32_t row, col = 0;
	row = data.rows();
	col = data.cols();

	m_dims = { row,col };
	QuantumObject(m_data, m_type, m_dims);
}

QuantumObject::QuantumObject(Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>& data, type t, std::pair<int32_t, int32_t> dims)
	:m_type(t), m_dims(dims)
{
	m_data = Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>(m_dims.first, m_dims.second);
	m_data = data;
}

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> QuantumObject::get_data()
{
	return m_data;
}

std::pair<int32_t, int32_t> QuantumObject::get_dims()
{
	return m_dims;
}

type QuantumObject::get_type()
{
	return m_type;
}
