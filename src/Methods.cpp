#include "Methods.h"
#include <mutex>
#include "constants.h"

static std::array<char, 3> axis = { 'x', 'y', 'z' };
static int num_jobs = 0;
//static std::vector<std::pair<int, std::vector<std::complex<double>>>> work;
static std::vector<std::pair<int, Eigen::SparseVector<std::complex<double>>>> work;
static std::vector<Eigen::SparseVector<std::complex<double>>> work_col;
static std::mutex work_lock;
static std::mutex cols_lock;
static std::vector<std::pair<int, Eigen::SparseVector<std::complex<double>>>> cols_vec;
static std::vector<std::pair<int, std::complex<double>>> DotReturnVec;
static ThreadPool* pool = nullptr;
int BaseSize = 32;

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
std::vector<std::complex<double>> DotProductVec(Matrix& mat, const std::vector<std::complex<double>>& vec, bool term)
{
	unsigned int size = vec.size();

	if (pool == nullptr)
	{	
		unsigned int nthreads = std::thread::hardware_concurrency();
		size = mat.rows();
		StartThreadPool(nthreads);
		work.clear();
		for (int i = 0; i < size; i++)
		{
			//std::vector<std::complex<double>> temp;
			Eigen::SparseVector<std::complex<double>> temp2(size);
			temp2 = mat.row(i);
			//for (int i = 0; i < size; i++)
			//{
			//	temp.push_back(temp2.coeff(i));
			//}
			work.push_back({ i,temp2 });
		}
	}

	if (term)
	{
		return { {std::complex<double>(0,0)} };
	}

	auto dot = [&](int i)
		{
			std::complex<double> sum = 0;

			int ind_size = work[i].second.data().size();
			for (int e = 0; e < ind_size; e++)
			{
				int index = work[i].second.data().index(e);
				
				auto val1 = work[i].second.coeff(index), val2 = vec[index];
				sum = sum + (val1 * val2);
			}

			//for (int e = 0; e < work[i].second.size(); e++)
			//{
			//	auto val1 = work[i].second.coeff(e), val2 = vec[e];
			//	const uint32_t bits = *(reinterpret_cast<uint32_t*>(&val1));
			//	if ((bits + bits) == 0)
			//	{
			//		sum = sum + std::complex<double>(0.0,0.0);
			//	}
			//	else
			//	{
			//		sum = sum + val1 * val2;
			//	}
			//
			//}

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

	typedef std::pair<int, std::complex<double>> TempVecType;


	auto SortPair = [&](TempVecType a, TempVecType b)
		{
			return b.first > a.first;
		};

	//std::vector<TempVecType> TempVec;
	//TempVec.reserve(DotReturnVec.size() * sizeof(TempVecType));
	//for (int i = 0; i < DotReturnVec.size(); i++)
	//{
	//	TempVec.push_back({ DotReturnVec[i].first, i });
	//}

	//std::sort(TempVec.begin(), TempVec.end(), SortPair);
	std::sort(DotReturnVec.begin(), DotReturnVec.end(), SortPair);
#if _DEBUG
	for (int i = 0; i < DotReturnVec.size(); i++)
	{
		if (DotReturnVec[i].first != i)
		{
			std::cin.get();
		}
	}
#endif

	std::vector<std::complex<double>> ReturnVec;

	for (int i = 0; i < DotReturnVec.size(); i++)
	{
		//ReturnVec.push_back(DotReturnVec[TempVec[i].second].second);
		ReturnVec.push_back(DotReturnVec[i].second);
	}

	return ReturnVec;
}

void StartThreadPool(int max_threads)
{
	pool = new ThreadPool(max_threads);
}

void KillThreadPool()
{
	if (pool == nullptr)
	{
		return;
	}
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

std::vector<std::array<double, 3>> FibonacciSphere(int n)
{
	std::vector<std::array<double, 3>> points;
	double phi = M_PI * (3.0 - std::sqrt(5.0)); //golden angle in radians

	for (int i = 0; i < n; i++)
	{
		double y = 1.0 - ((double)i / double(n - 1)) * 2;
		double r = std::sqrt(1.0 - (y * y));

		double theta = phi * (double)i;

		double x = std::cos(theta) * r;
		double z = std::sin(theta) * r;

		points.push_back({x,y,z});
	}

	return points;
}

Matrix BlockInverse(Matrix mat, int dim, int inner_block_size) //finds the inverse of a matrix that's been built using a block layout 
{
	std::vector<Matrix> blocks;

	std::vector<std::pair<int, int>> dimensions;
	dimensions.push_back({ 1,1 });
	dimensions.push_back({ dim - 1, 1 });
	dimensions.push_back({ 1, dim - 1 });
	dimensions.push_back({ dim - 1, dim - 1 });

	//doing the biggest block first

	typedef Eigen::Triplet<std::complex<double>, int32_t> T;
	std::vector<T> entries;

	int RightLowerDim = (dim - 1) * inner_block_size;
	Matrix S(RightLowerDim, RightLowerDim);
	for (int i = inner_block_size; i < (dim * inner_block_size); i++)
	{
		Eigen::SparseVector<std::complex<double>> row = mat.row(i);
		int size = row.data().size();
		for (int e = 0; e < size; e++)
		{
			int index = row.data().index(e);

			if (index < inner_block_size)
			{
				continue;
			}

			std::complex<double> val = row.coeff(index);
			entries.push_back(T(i - inner_block_size, index - inner_block_size, val));
		}
	}
	S.setFromTriplets(entries.begin(), entries.end());

	entries.clear();
	Matrix B_00_inverse(inner_block_size, inner_block_size);
	for (int i = 0; i < inner_block_size; i++)
	{
		Eigen::SparseVector<std::complex<double>> row = mat.row(i);
		int size = row.data().size();
		for (int e = 0; e < size; e++)
		{
			int index = row.data().index(e);

			if (index >= inner_block_size)
			{
				continue;
			}

			std::complex<double> val = row.coeff(index);

			entries.push_back(T(i, index, val));
		}
	}
	B_00_inverse.setFromTriplets(entries.begin(), entries.end()); //B_00 is diagonal so the inverse is the reciprocal of the matrix elements 
	for (int i = 0; i < inner_block_size; i++)
	{
		B_00_inverse.coeffRef(i, i) = std::complex(1.0) / (std::complex<double>)B_00_inverse.coeff(i, i);
	}

	entries.clear();
	Matrix B_01(inner_block_size, dimensions[1].first * inner_block_size);
	for (int i = 0; i < inner_block_size; i++)
	{
		Eigen::SparseVector<std::complex<double>> row = mat.row(i);
		int size = row.data().size();
		for (int e = 0; e < size; e++)
		{
			int index = row.data().index(e);

			if (index < inner_block_size)
			{
				continue;
			}

			std::complex<double> val = row.coeff(index);

			entries.push_back(T(i, index - inner_block_size, val));
		}
	}
	B_01.setFromTriplets(entries.begin(), entries.end());

	entries.clear();
	Matrix B_10(dimensions[2].second * inner_block_size, inner_block_size);
	for (int i = inner_block_size; i < dim * inner_block_size; i++)
	{
		Eigen::SparseVector<std::complex<double>> row = mat.row(i);
		int size = row.data().size();
		for (int e = 0; e < size; e++)
		{
			int index = row.data().index(e);

			if (index >= inner_block_size)
			{
				continue;
			}

			std::complex<double> val = row.coeff(index);

			entries.push_back(T(i - inner_block_size, index, val));
		}
	}
	B_10.setFromTriplets(entries.begin(), entries.end());
	entries.clear();

	S = S - (B_10 * B_00_inverse * B_01);
	Matrix S_inverse(RightLowerDim, RightLowerDim);
	if (dim > 2) //not in 2x2 block form 
	{
		//S_inverse = BlockInverse(S, dimensions[3].first, inner_block_size);
	}
	else
	{
		for (int i = 0; i < inner_block_size; i++)
		{
			S_inverse.coeffRef(i, i) = std::complex(1.0) / (std::complex<double>)S.coeff(i, i);
		}
	}
	blocks = { B_00_inverse, B_01, B_10, S_inverse };

	std::vector<Matrix> corners;

	{
		Matrix q1 = B_00_inverse + (B_00_inverse * B_01 * S_inverse * B_10 * B_00_inverse);
		Matrix q2 = std::complex(-1.0) * B_00_inverse * B_01 * S_inverse;
		Matrix q3 = std::complex(-1.0) * S_inverse * B_10 * B_00_inverse;
		//Matrix q4 = S_inverse

		corners = { q1,q2,q3,S_inverse };
	}

	entries.clear();

	for (int i = 0; i < 4; i++)
	{
		int rows = corners[i].rows();
		for (int j = 0; j < rows; j++)
		{
			Eigen::SparseVector<std::complex<double>> row = corners[i].row(j);
			int size = row.data().size();
			for (int k = 0; k < size; k++)
			{
				int index = row.data().index(k);
				entries.push_back(T((std::floor((double)i / 2.0) * inner_block_size) + j, ((i) % 2 * inner_block_size) + index, row.coeff(index)));
			}
		}
	}

	Matrix InvertedMat(dim* inner_block_size, dim* inner_block_size);
	InvertedMat.setFromTriplets(entries.begin(), entries.end());

	return InvertedMat;
}

bool Diagonal(const Matrix* mat)
{
	int rows = mat->rows();
	for (int i = 0; i < rows; i++)
	{
		Eigen::SparseVector<std::complex<double>> row = mat->row(i);
		if (row.data().size() > 1)
		{
			return false;
		}
		
		if (row.data().index(0) != i)
		{
			return false;
		}
	}

	return true;
}

Matrix GetInverse(const Matrix* mat, bool diagonal)
{
	int rows = mat->rows();
	Matrix inverse(rows, rows);

	typedef Eigen::Triplet<std::complex<double>, int32_t> T;
	std::vector<T> entries;

	if (diagonal)
	{
		for (int i = 0; i < rows; i++)
		{
			entries.push_back(T(i, i, std::complex(1.0) / (std::complex<double>)mat->coeff(i, i)));
		}
		inverse.setFromTriplets(entries.begin(), entries.end());
		return inverse;
	}

	//Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>> solver;
	//solver.compute(*mat);
	auto ide = identity_col(rows);
	work_col.clear();

	auto SolveCol = [mat](int i) //update this code to use the LU decompastition rather than the doing the decomposition each time
		{
			Eigen::SparseLU<Matrix> solver;
			solver.compute(*mat);
			Eigen::SparseVector<std::complex<double>> xCol = solver.solve(work_col[i]);
			cols_lock.lock();
			cols_vec.push_back({ i,xCol });
			cols_lock.unlock();
		};

	ThreadPool* SolvePool = new ThreadPool(std::thread::hardware_concurrency());

	for (int i = 0; i < ide.cols(); i++)
	{
		work_col.push_back(ide.col(i));
		SolvePool->QueueJob(SolveCol, i);
		//cols.push_back(solver.solve(ide.col(i)));
	}

	SolvePool->start();
	while (SolvePool->Busy()) {};
	SolvePool->Stop();
	work_col.clear();
	
	delete SolvePool;

	typedef std::pair<int, Eigen::SparseVector<std::complex<double>>> TempVecType;


	auto SortPair = [&](TempVecType a, TempVecType b)
		{
			return b.first > a.first;
		};
	std::sort(cols_vec.begin(), cols_vec.end(), SortPair);

	for (int i = 0; i < cols_vec.size(); i++)
	{
		for (int e = 0; e < cols_vec[i].second.data().size(); e++)
		{
			int index = cols_vec[i].second.data().index(e);
			std::complex<double> val = cols_vec[i].second.coeff(index);
			if (std::abs(val) <= 1e-6)
			{
				continue;
			}
			entries.push_back(T(index, i, val));
		}
	}
	inverse.setFromTriplets(entries.begin(), entries.end());
	cols_vec.clear();
	return inverse;
}

Matrix BlockInverse(const Matrix& mat, int dim) //Mat - matrix to invert, dim - number of rows/columns 
{
	std::vector<Matrix> blocks;

	int BlockSize = 0;
	BlockSize = (int)std::floor((double)mat.rows() / 2.0);

	std::vector<std::pair<int, int>> dimensions;

	dimensions.push_back({ BlockSize, BlockSize });
	dimensions.push_back({ BlockSize, dim - BlockSize });
	dimensions.push_back({ dim - BlockSize, BlockSize });
	dimensions.push_back({ dim - BlockSize, dim - BlockSize });

	typedef Eigen::Triplet<std::complex<double>, int32_t> T;
	std::vector<T> entries;

	int RightLowerDim = dimensions[3].first;
	Matrix S(RightLowerDim, RightLowerDim);
	for (int i = BlockSize; i < dim; i++)
	{
		Eigen::SparseVector<std::complex<double>> row = mat.row(i);
		int size = row.data().size();
		for (int e = 0; e < size; e++)
		{
			int index = row.data().index(e);

			if (index < BlockSize)
			{
				continue;
			}

			std::complex<double> val = row.coeff(index);
			entries.push_back(T(i - BlockSize, index - BlockSize, val));
		}
	}
	S.setFromTriplets(entries.begin(), entries.end());

	entries.clear();
	Matrix B_00(BlockSize, BlockSize);
	for (int i = 0; i < BlockSize; i++)
	{
		Eigen::SparseVector<std::complex<double>> row = mat.row(i);
		int size = row.data().size();
		for (int e = 0; e < size; e++)
		{
			int index = row.data().index(e);

			if (index >= BlockSize)
			{
				continue;
			}

			std::complex<double> val = row.coeff(index);

			entries.push_back(T(i, index, val));
		}
	}
	B_00.setFromTriplets(entries.begin(), entries.end());
	
	bool diag = Diagonal(&B_00);
	Matrix B_00_inverse(BlockSize, BlockSize);
	if (!diag && B_00.rows() > BaseSize)
	{
		B_00_inverse = BlockInverse(B_00, B_00.rows());
	}
	else if (diag)
	{
		B_00_inverse = GetInverse(&B_00, true);
	}
	else
	{
		B_00_inverse = GetInverse(&B_00, false);
	}

	{
		Matrix empty(BlockSize, BlockSize);
		B_00 = B_00 * empty;
	}

	entries.clear();
	Matrix B_01(dimensions[1].first, dimensions[1].second);
	for (int i = 0; i < BlockSize; i++)
	{
		Eigen::SparseVector<std::complex<double>> row = mat.row(i);
		int size = row.data().size();
		for (int e = 0; e < size; e++)
		{
			int index = row.data().index(e);

			if (index < BlockSize)
			{
				continue;
			}

			std::complex<double> val = row.coeff(index);

			entries.push_back(T(i, index - BlockSize, val));
		}
	}
	B_01.setFromTriplets(entries.begin(), entries.end());

	entries.clear();
	Matrix B_10(dimensions[2].first, dimensions[2].second);
	for (int i = BlockSize; i < dim; i++)
	{
		Eigen::SparseVector<std::complex<double>> row = mat.row(i);
		int size = row.data().size();
		for (int e = 0; e < size; e++)
		{
			int index = row.data().index(e);
	
			if (index >= BlockSize)
			{
				continue;
			}
	
			std::complex<double> val = row.coeff(index);
	
			entries.push_back(T(i - BlockSize, index, val));
		}
	}
	B_10.setFromTriplets(entries.begin(), entries.end());
	entries.clear();

	S = S - (B_10 * B_00_inverse * B_01);
	diag = Diagonal(&S);
	Matrix S_inverse(dimensions[3].first, dimensions[3].second);
	if (!diag && B_00.rows() > BaseSize)
	{
		S_inverse = BlockInverse(S, S.rows());
	}
	else if (diag)
	{
		S_inverse = GetInverse(&S, true);
	}
	else
	{
		S_inverse = GetInverse(&S, false);
	}

	{
		Matrix empty(S.rows(), S.cols());
		S = S * empty;
	}

	std::vector<Matrix> corners;

	{
		Matrix q1 = B_00_inverse + (B_00_inverse * B_01 * S_inverse * B_10 * B_00_inverse);
		Matrix q2 = std::complex(-1.0) * B_00_inverse * B_01 * S_inverse;
		Matrix q3 = std::complex(-1.0) * S_inverse * B_10 * B_00_inverse;
		//Matrix q4 = S_inverse

		corners = { q1,q2,q3,S_inverse };
	}

	entries.clear();

	for (int i = 0; i < 4; i++)
	{
		int rows = corners[i].rows();
		for (int j = 0; j < rows; j++)
		{
			Eigen::SparseVector<std::complex<double>> row = corners[i].row(j);
			int size = row.data().size();
			for (int k = 0; k < size; k++)
			{
				int index = row.data().index(k);
				
				int rmod = 0, cmod = 0;
				switch (i)
				{
				case 1:
					rmod = 0;
					cmod = BlockSize;
					break;
				case 2:
					rmod = BlockSize;
					cmod = 0;
					break;
				case 3:
					rmod = BlockSize;
					cmod = BlockSize;
					break;
				}
				entries.push_back(T(rmod + j, cmod + index, row.coeff(index)));
			}
		}
	}

	Matrix InvertedMat(dim , dim);
	InvertedMat.setFromTriplets(entries.begin(), entries.end());
	return InvertedMat;
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
