#include "Methods.h"
#include <mutex>
#include "constants.h"

static std::array<char, 3> axis = { 'x', 'y', 'z' };
static int num_jobs = 0;
static std::vector<std::pair<int, std::vector<std::complex<double>>>> work;
static std::mutex work_lock;
static std::vector<std::pair<int, std::complex<double>>> DotReturnVec;
static ThreadPool* pool = nullptr;

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> commutator(Matrix& A, Matrix& B)
{
	if (A.rows() != B.rows() and B.cols() != A.cols())
	{
		std::cout << "MATRICIES NOT THE SAME DIMENSION" << std::endl;
		return Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>();
	}

	Matrix ReturnMat(A.rows(), A.cols());
	ReturnMat = (A * B) - (B * A);
	return ReturnMat;
}

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> MakeHamiltonian(std::vector<int32_t> dims, int ind, std::array<double, 3> parvec)
{
	std::vector<Matrix> components = {};

	int dim = 1;
	int size = dims.size();
	for (int i = 0; i < size; i++)
	{
		dim = dim * dims[i];
	}

	Matrix ReturnMatrix(dim, dim);

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
		Matrix c = v * MakeSpinOperator(dims, { spec });
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
	std::vector<Matrix> components = {};
	int dim = 1;
	int size = dims.size();
	for (int i = 0; i < size; i++)
	{
		dim = dim * dims[i];
	}

	Matrix ReturnMatrix(dim, dim);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (parmat[i][j] != 0)
			{
				Matrix c = parmat[i][j] * MakeSpinOperator(dims, { {ind_1, axis[i]}, {ind_2, axis[j]} });
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

Eigen::VectorXcd FlattenMatrix(Matrix mat)
{

	int row = mat.rows();
	int col = mat.cols();
	Eigen::VectorXcd ReturnVec(row * col);
	for (int i = 0; i < row; i++)
	{
		int offset = i * col;
		for (int e = 0; e < col; e++)
		{
			ReturnVec[offset + e] = mat.coeff(i, e);
		}
	}
	return ReturnVec;
}

std::vector<std::complex<double>> FlattenMatrixVec(Matrix mat)
{
	int row = mat.rows();
	int col = mat.cols();
	std::vector<std::complex<double>> ReturnVec;
	ReturnVec.reserve(row * col);
	for (int i = 0; i < row; i++)
	{
		for (int e = 0; e < col; e++)
		{
			ReturnVec.push_back(mat.coeff(i, e));
		}
	}
	return ReturnVec;
}

Matrix ReformMatrix(std::vector<std::complex<double>> vec)
{
	int size = vec.size();
	int dim = (int)std::sqrt(size);
	Matrix rho(dim, dim);
	
	typedef Eigen::Triplet<std::complex<double>, int32_t> T;
	std::vector<T> entries;
	
	for (int i = 0; i < dim; i++)
	{
		int offset = i * dim;
		for (int e = 0; e < dim; e++)
		{
			if (vec[e + offset] != (std::complex<double>)0)
			{
				entries.push_back(T(i, e, vec[e + offset]));
			}
		}
	}

	rho.setFromTriplets(entries.begin(), entries.end());
	return rho;
}

Eigen::VectorXcd DotProduct(Matrix& mat, const Eigen::VectorXcd& vec)
{
	Eigen::VectorXcd ReturnVec(vec.rows());

	for (int i = 0; i < mat.rows(); i++)
	{
		std::complex<double> sum = 0;
		for (int e = 0; e < mat.cols(); e++)
		{
			sum = sum + mat.coeff(i, e) * vec[e];
		}
		ReturnVec[i] = sum;
	}
	return ReturnVec;
}



//std::vector<std::complex<double>> DotProductVec(Eigen::MatrixXcd& mat, const std::vector<std::complex<double>>& vec)
std::vector<std::complex<double>> DotProductVec(Matrix& mat, const std::vector<std::complex<double>>& vec)
{
	unsigned int size = vec.size();

	if (pool == nullptr)
	{	
		unsigned int nthreads = std::thread::hardware_concurrency();
		StartThreadPool(nthreads);
		work.clear();
		for (int i = 0; i < size; i++)
		{
			std::vector<std::complex<double>> temp;
			auto temp2 = mat.row(i);
			for (int i = 0; i < size; i++)
			{
				temp.push_back(temp2.coeff(i));
			}
			work.push_back({ i,temp });
		}
	}

	auto dot = [&](int i)
		{
			std::complex<double> sum = 0;
			for (int e = 0; e < work[i].second.size(); e++)
			{
				auto val1 = work[i].second[e], val2 = vec[e];
				const uint32_t bits = *(reinterpret_cast<uint32_t*>(&val1));
				if ((bits + bits) == 0)
				{
					sum = sum + std::complex<double>(0.0,0.0);
				}
				else
				{
					sum = sum + work[i].second[e] * val2;
				}
			}
			work_lock.lock();
			DotReturnVec.push_back({ i, sum });
			work_lock.unlock();
		};

	DotReturnVec.clear();

	for (int i = 0; i < work.size(); i++)
	{
		pool->QueueJob(dot, i);
	}

	pool->start();
	while (pool->Busy()) {};
	pool->Stop();

	typedef std::pair<int, int> TempVecType;


	auto SortPair = [&](TempVecType a, TempVecType b)
		{
			return b.first > a.first;
		};

	std::vector<TempVecType> TempVec;
	TempVec.reserve(DotReturnVec.size() * sizeof(TempVecType));
	for (int i = 0; i < DotReturnVec.size(); i++)
	{
		TempVec.push_back({ DotReturnVec[i].first, i });
	}

	std::sort(TempVec.begin(), TempVec.end(), SortPair);

	//MergeSort(DotReturnVec, 0, DotReturnVec.size() - 1);
	//sort(DotReturnVec); //put quicksort in tommorow
	std::vector<std::complex<double>> ReturnVec;

	for (int i = 0; i < DotReturnVec.size(); i++)
	{
		ReturnVec.push_back(DotReturnVec[TempVec[i].second].second);
	}

	return ReturnVec;
}

void StartThreadPool(int max_threads)
{
	pool = new ThreadPool(max_threads);
}

void KillThreadPool()
{
	pool->Stop();
	delete pool;
	pool = nullptr;
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

Eigen::Vector4d SingletState()
{
	Eigen::Vector4d singlet;

	Eigen::Vector2d up = { 1,0 };
	Eigen::Vector2d down = { 0,1 };

	auto sum = [](int a, int b) { return a + std::pow(b, 2); };

	Eigen::Vector4d SingletState = TensorProduct<double, 2, 2>(up, down) - TensorProduct<double, 2, 2>(down, up);
	double NormFactor = 1.0 / std::sqrt((double)std::reduce(SingletState.begin(), SingletState.end(), 0, sum));
	SingletState = NormFactor * SingletState;

	return SingletState;
}

double simpson_integration(std::vector<double> x_list, std::vector<double> y_list)
{
	double area = 0;
	for (int i = 0; i < x_list.size()-1; i++)
	{
		double diff = x_list[i + 1] - x_list[i];
		double ab = y_list[i] + y_list[i + 1];

		area = area + (ab * 0.5) * diff;
	}
	return area;
}

void sort(std::vector<std::pair<int, std::complex<double>>>& arr)
{
	auto sortpair = [&](std::pair<int, std::complex<double>> a, std::pair<int, std::complex<double>>b)
		{
			return b.first > a.first;
		};

	std::vector<std::pair<int, std::complex<double>>>::iterator b = arr.begin();
	std::vector<std::pair<int, std::complex<double>>>::iterator e = arr.end();

	std::sort(b, e, sortpair);
}
