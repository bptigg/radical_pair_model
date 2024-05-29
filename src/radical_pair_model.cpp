#include <Eigen/Core>
#include <random>
#include <iostream>

//#include "spin_operator.h"
#include "Methods.h"
#include "constants.h"

int main() {
	//std::default_random_engine generator;
	//std::poisson_distribution<int> distribution(4.1);
	//auto poisson = [&]() { return distribution(generator); };
	//
	//Eigen::RowVectorXi v = Eigen::RowVectorXi::NullaryExpr(10, poisson);
	////std::cout << v << "\n";
	//
	//auto mat = PointDipoleDipoleCoupling({ 8.51061, -14.251621, 6.5492562 });
	//
	//char x = 'x';
	//auto spin_x = jmat(2.5, &x);
	//std::cout << Eigen::MatrixXcd(spin_x) << std::endl;
	//
	//std::cout << "------------------------------------" << std::endl;
	//
	//x = 'z';
	//spin_x = jmat(1.5, &x);
	//std::cout << Eigen::MatrixXcd(spin_x) << std::endl;
	//
	//auto i = identity({ 2, 2 });
	//
	//std::cout << Eigen::MatrixXcd(i) << std::endl;
	//
	//std::vector<int32_t> dims = { 2,2,3,3 };
	//std::vector<specs> spec = { { 0, 'x' } , { 1, 'y' } };
	//
	//auto b = MakeSpinOperator(dims, spec);
	//
	//std::cout << Eigen::MatrixXcd(b) << std::endl;
	//
	//auto c = commutator(b, spin_x);
	//std::cout << Eigen::MatrixXcd(c) << std::endl;
	//
	//Eigen::SparseMatrix<double, Eigen::ColMajor> test(4, 2);
	//std::cout << test.rows() << " , " << test.cols() << std::endl;
	//
	//auto testmat = MakeHamiltonian(dims, 0, {0.0, 0.0, 8.79645943005142 });
	//std::cout << Eigen::MatrixXcd(testmat) << std::endl;

	auto DDSE = PointDipoleDipoleCoupling({ 8.51061, -14.251621, 6.5492562 });
	std::vector<int32_t>dims = { 2,2,3,3 };

	MATRIX3x3 ide;

	for (int i = 0; i < 3; i++)
	{
		for (int e = 0; e < 3; e++)
		{
			DDSE[i][e] = DDSE[i][e] * 2 * M_PI;
			if (i == e)
			{
				ide[i][e] = 1;
				continue;
			}
			ide[i][e] = 0;
		}
		
	}

	MATRIX intitial_projection_operator = 0.25 * MakeSpinOperator(dims, {}) - MakeHamiltonian(dims, 0, 1, ide);
	MATRIX identity = MakeSpinOperator(dims, {});
	MATRIX triplet_projection_operator = identity - intitial_projection_operator;

	MATRIX DipoleDipoleElectronHamiltonian = MakeHamiltonian(dims, 0, 1, DDSE);

	//std::cout << Eigen::MatrixXcd(intitial_projection_operator) << std::endl;
	//std::cout << Eigen::MatrixXcd(triplet_projection_operator) << std::endl;
	std::cout << Eigen::MatrixXcd(DipoleDipoleElectronHamiltonian) << std::endl;

	MATRIX H1 = MakeHamiltonian(dims, 0, 2, N5);
	MATRIX H2 = MakeHamiltonian(dims, 1, 3, N1);
	MATRIX HyperfineHamiltonian = H1 + H2;

	MATRIX H_0 = HyperfineHamiltonian + DipoleDipoleElectronHamiltonian;
	std::cout << Eigen::MatrixXcd(H_0) << std::endl;


	std::cin.get();
}

