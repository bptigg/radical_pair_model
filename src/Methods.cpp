#include "Methods.h"

#include "constants.h"

static std::array<char, 3> axis = { 'x', 'y', 'z' };

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> commutator(MATRIX& A, MATRIX& B)
{
	if (A.rows() != B.rows() and B.cols() != A.cols())
	{
		std::cout << "MATRICIES NOT THE SAME DIMENSION" << std::endl;
		return Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>();
	}

	MATRIX ReturnMat(A.rows(), A.cols());
	ReturnMat = (A * B) - (B * A);
	return ReturnMat;
}

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeHamiltonian(std::vector<int32_t> dims, int ind, std::array<double, 3> parvec)
{
	std::vector<MATRIX> components = {};

	int dim = 1;
	int size = dims.size();
	for (int i = 0; i < size; i++)
	{
		dim = dim * dims[i];
	}

	MATRIX ReturnMatrix(dim, dim);

	std::vector<std::pair<float, char>> comp = {};

	for (int i = 0; i < 3; i++)
	{
		if (parvec[i] != 0)
		{
			comp.push_back({ parvec[i], axis[i] });
		}
	}
	
	for (auto[v, ax] : comp)
	{
		specs spec = { ind, ax };
		MATRIX c = v * MakeSpinOperator(dims, { spec });
		components.push_back(c);
	}

	if (components.size() != 0)
	{
		for (auto mat : components)
		{
			ReturnMatrix = ReturnMatrix + mat;
		}
	}
	else
	{
		ReturnMatrix = MakeZeroOperator(dims);
	}

	return ReturnMatrix;
}


Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeHamiltonian(std::vector<int32_t> dims, int ind_1, int ind_2, std::array<std::array<double, 3>, 3> parmat)
{
	std::vector<MATRIX> components = {};
	int dim = 1;
	int size = dims.size();
	for (int i = 0; i < size; i++)
	{
		dim = dim * dims[i];
	}

	MATRIX ReturnMatrix(dim, dim);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (parmat[i][j] != 0)
			{
				MATRIX c = parmat[i][j] * MakeSpinOperator(dims, { {ind_1, axis[i]}, {ind_2, axis[j]} });
				components.push_back(c);
			}
		}
	}

	if (components.size() != 0)
	{
		for (auto mat : components)
		{
			ReturnMatrix = ReturnMatrix + mat;
		}
	}
	else
	{
		ReturnMatrix = MakeZeroOperator(dims);
	}

	return ReturnMatrix;
}

MATRIX3x3 PointDipoleDipoleCoupling(double r)
{
	double C = -1 * ((GSL_CONST_MKSA_VACUUM_PERMEABILITY * M_1_PI) / (4e-30)) * std::pow((G_FACTOR_ELECTRON * GSL_CONST_MKSA_BOHR_MAGNETON), 2) * (1 / (1e6 * GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR * 2 * M_PI));
	double d = C / std::pow(r, 3);

	MATRIX3x3 A;
	A[0][0] = -d;
	A[1][1] = -d;
	A[2][2] = -2 * d;
	
	return A;
}

MATRIX3x3 PointDipoleDipoleCoupling(std::array<double, 3> r)
{
	double C = -1 * ((GSL_CONST_MKSA_VACUUM_PERMEABILITY * M_1_PI) / (4e-30)) * std::pow((G_FACTOR_ELECTRON * GSL_CONST_MKSA_BOHR_MAGNETON), 2) * (1 / (1e6 * GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR * 2 * M_PI));
	double r_norm = 0;

	for (int i = 0; i < 3; i++)
	{
		r_norm = r_norm + std::pow(r[i], 2);
	}

	r_norm = std::sqrt(r_norm);
	double d = C / std::pow(r_norm, 3);
	std::array<double, 3>e;
	for (int i = 0; i < 3; i++)
	{
		e[i] = r[i] / r_norm;
	}

	MATRIX3x3 ReturnMatrix;

	for (int i = 0; i < 3; i++)
	{
		for (int i2 = 0; i2 < 3; i2++)
		{
			ReturnMatrix[i][i2] = 3*d*(e[i] * e[i2]);
			if (i == i2)
			{
				ReturnMatrix[i][i2] = ReturnMatrix[i][i2] - d;
			}
		}
	}
	
	return ReturnMatrix;
}
