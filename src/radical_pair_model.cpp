#include <Eigen/Core>
#include <random>
#include <iostream>

#include "spin_operator.h"

int main() {
	std::default_random_engine generator;
	std::poisson_distribution<int> distribution(4.1);
	auto poisson = [&]() { return distribution(generator); };

	Eigen::RowVectorXi v = Eigen::RowVectorXi::NullaryExpr(10, poisson);
	//std::cout << v << "\n";

	char x = 'x';
	auto spin_x = jmat(2.5, &x);
	std::cout << Eigen::MatrixXcd(spin_x) << std::endl;

	std::cout << "------------------------------------" << std::endl;

	x = 'z';
	spin_x = jmat(2.5, &x);
	std::cout << Eigen::MatrixXcd(spin_x) << std::endl;
}


