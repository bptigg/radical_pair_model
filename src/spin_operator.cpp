#include "spin_operator.h"

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> jmat(double j, char* args)
{
	if ((2.0 * j) != (int)(2.0 * j) or (j < 0)) {
		std::cout << "j must be a half integer or a non negative integer" << std::endl;
		//return nullptr
	}

	if (args == nullptr) {
		char x = 'x', y = 'y', z = 'z';
		return jmat(j, &x), jmat(j, &y), jmat(j, &z);
	}

	Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> matrix;

	if (args[0] == '+') {
		matrix = _jplus(j);
	}
	else if (args[0] == '-') {
		auto TempMatrix = _jplus(j);
		matrix = get_hermition(TempMatrix);
	}
	else if (args[0] == 'x') {
		auto TempMatrix = _jplus(j);
		matrix = 0.5 * (TempMatrix + get_hermition(TempMatrix));
	}
	else if (args[0] == 'y') {
		auto TempMatrix = _jplus(j);
		matrix = -0.5 * j * (TempMatrix - get_hermition(TempMatrix));
	}
	else if (args[0] == 'z') {
		matrix = _jz(j);
	}
	else {
		std::cout << "Invalid type" << std::endl;
	}
	return matrix;
}

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> _jplus(double j)
{
	std::vector<std::complex<double>> m = {};
	std::vector<std::complex<double>> data = {};
	int diff = 2 * j + 1;
	for (int i = 0; i < diff; i++) {
		m.push_back((std::complex<double>)(j - (double)i));
		
		if (i == 0) {
			continue;
		}
		std::complex<double> val = std::sqrt((j * (j + 1))- (m[i] + std::complex<double>(1)) * m[i]);
		data.push_back(val);
	}

	int N = m.size();
	int32_t* ind = (int32_t*)malloc(N * sizeof(int32_t));
	for (int i = 1; i < N; i++)
	{
		ind[i - 1] = i;
	}

	int32_t* ptr = (int32_t*)malloc((N + 1) * sizeof(int32_t));
	for (int32_t i = 0; i < N - 1; i++)
	{
		ptr[i] = i;
	}
	ptr[N - 1]	= (int32_t)(N - 1);
	ptr[N]		= (int32_t)(N - 1);

	typedef Eigen::Triplet<std::complex<double>, int32_t> T;
	std::vector<T> entries;
	int32_t row = 0;
	int current_row = 0;
	std::vector<int32_t> rows;
	while (rows.size() < data.size())
	{
		if(current_row + 1 == ptr[row + 1] - ptr[row])
		{
			rows.push_back(row);
			current_row = current_row + 1;
		}
		else
		{
			row = row + 1;
			current_row = 0;
		}
	}

	for (int i = 0; i < data.size(); i++)
	{
		entries.push_back(T(rows[i], ind[i], data[i]));
	}
	Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> mat(N, N);
	mat.setFromTriplets(entries.begin(), entries.end());

	free(ind);
	free(ptr);

	return mat;

}

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> _jz(double j)
{
	int N = (int)((2 * j) + 1);

	std::vector<std::complex<double>> data;
	for (int i = 0; i < N; i++)
	{
		if (j - (double)i == 0)
		{
			continue;
		}
		else
		{
			data.push_back(std::complex<double>(j - (double)i, 0));
		}
	}

	std::vector<int32_t> ind;
	std::vector<int32_t> ptr;

	if (N % 2 == 0) {
		ind = std::vector<int32_t>(N);
		std::iota(std::begin(ind), std::end(ind), 0);
		
		ptr = std::vector<int32_t>(N + 1);
		std::iota(std::begin(ptr), std::end(ptr), 0);
		ptr[N] = N;
	}
	else
	{
		j = (int)j;
		for (int32_t i = 0; i < N; i++)
		{
			if ((double)i == j)
			{
				continue;
			}

			ind.push_back(i);
		}
		bool go_back = true;
		for (int32_t i = 0; i < N; i++)
		{
			ptr.push_back(i);
			if ((double)i == j && go_back)
			{
				i = i - 1;
				go_back = false;
			}
		}
		ptr[N] = N - 1;
	}

	typedef Eigen::Triplet<std::complex<double>, int32_t> T;
	std::vector<T> entries;
	int32_t row = 0;
	int current_row = 0;
	std::vector<int32_t> rows;
	while(rows.size() < data.size())
	{
		if (current_row + 1 <= ptr[row + 1] - ptr[row])
		{
			rows.push_back(row);
			current_row = current_row + 1;
		}
		else
		{
			row = row + 1;
			current_row = 0;
		}
	}
	for (int i = 0; i < data.size(); i++)
	{
		entries.push_back(T(rows[i], ind[i], data[i]));
	}
	Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> mat(N, N);
	mat.setFromTriplets(entries.begin(), entries.end());

	return mat;
}

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> spin_jx(double j)
{
	char x = 'x';
	return jmat(j, &x);
}

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> spin_jy(double j)
{
	char x = 'y';
	return jmat(j, &x);
}

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> spin_jz(double j)
{
	char x = 'z';
	return jmat(j, &x);
}

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> identity(std::vector<int> dims)
{
	//Eigen::Tensor<double, dims.size()> IdentityTensor(dims);

	std::vector<Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>> IdentityMatrix = {};
	for (auto i : dims) 
	{
		IdentityMatrix.push_back(identity(i));
	}

	auto multiply = [&](int a, int b) {return a * b; };

	int dimension = std::reduce(dims.begin(), dims.end(), 1, multiply);
	Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> Identiy(dimension, dimension);

	std::vector<Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>> StepUp;

	for (int i = 0; i < dims.size(); i++)
	{
		dimension = std::reduce(dims.begin(), dims.end()-dims.size() + i+1, 1, multiply);
		StepUp.push_back(Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>(dimension, dimension));
	}

	StepUp[0] = IdentityMatrix[0];

	for (int i = 0; i < dims.size()-1 ; i++) 
	{
		StepUp[i+1] = Eigen::KroneckerProductSparse(StepUp[i], IdentityMatrix[i + 1]);
	}
	Identiy = StepUp[StepUp.size()-1];

	return Identiy;
}

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> identity(int dims)
{
	Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> identity(dims, dims);
	for (int i = 0; i < dims; i++)
	{
		identity.coeffRef(i, i) = (std::complex<double>)1.0;
	}
	//std::cout << Eigen::MatrixXcd(identity) << std::endl;
	return identity;
}

Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> get_hermition(Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>& matrix)
{
	Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> ConjugateTransposeMatrix = matrix.transpose().conjugate();
	return ConjugateTransposeMatrix;
}

