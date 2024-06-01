#include <Eigen/Core>
#include <random>
#include <iostream>
#include <tuple>
#include <fstream>
#include <ctime>
#include <chrono>
#include <string>

//#include "spin_operator.h"
#include "Methods.h"
#include "constants.h"

#include <Eigen/Dense>
#include <omp.h>
#include <boost\numeric\odeint.hpp>
#include <boost\numeric\odeint\external\eigen\eigen.hpp>
#include <boost\numeric\odeint\external\openmp\openmp.hpp>
#include <Eigen/Eigenvalues>

#pragma warning( disable : 4996)


typedef std::vector<std::complex<double>> state_type;
using namespace boost::numeric;
typedef odeint::runge_kutta4<state_type> stepper_type;


class QuantumMasterEquation
{
	Eigen::MatrixXcd m_supor;

public:
	QuantumMasterEquation(Eigen::MatrixXcd& SuperOperator)
		:m_supor(SuperOperator)
	{
	}

	void operator() (const state_type& rho, state_type& drhodt, const double t)
	{
		state_type vec = DotProductVec(m_supor, rho);
		for (int i = 0; i < vec.size(); i++)
		{
			vec[i] = vec[i] * std::complex<double>(-1,0);
		}
		drhodt = vec;
	}
};

struct observer
{
	std::vector<state_type>& m_states;
	std::vector<double>& m_times;

	observer(std::vector<state_type>& states, std::vector<double>& times)
		:m_states(states), m_times(times) {}

	void operator() (state_type& rho, double t)
	{
 		m_states.push_back(rho);
		m_times.push_back(t);

		if (m_times.size() % 10 == 0)
		{
			std::cout << "Current time step: " << m_times[m_times.size() - 1] << std::endl;
		}
	}

	std::tuple<std::vector<state_type>, std::vector<double>> get_data()
	{
		return std::make_tuple(m_states, m_times);
	}
};

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

	std::time_t t = std::time(0);
	std::tm now = *std::localtime(&t);
	char buf[20];
	strftime(buf, sizeof(buf), "%H_%M_%S__%d_%m_%y", &now);
	std::string time_string = buf;

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

	Matrix singlet_projection_operator = 0.25 * MakeSpinOperator(dims, {}) - MakeHamiltonian(dims, 0, 1, ide);
	Matrix identity = MakeSpinOperator(dims, {});
	Matrix triplet_projection_operator = identity - singlet_projection_operator;

	Matrix DipoleDipoleElectronHamiltonian = MakeHamiltonian(dims, 0, 1, DDSE);

	//std::cout << Eigen::MatrixXcd(intitial_projection_operator) << std::endl;
	//std::cout << Eigen::MatrixXcd(triplet_projection_operator) << std::endl;
	//std::cout << Eigen::MatrixXcd(DipoleDipoleElectronHamiltonian) << std::endl;

	Matrix H1 = MakeHamiltonian(dims, 0, 2, N5);
	Matrix H2 = MakeHamiltonian(dims, 1, 3, N1);
	Matrix HyperfineHamiltonian = H1 + H2;

	Matrix H_0 = HyperfineHamiltonian + DipoleDipoleElectronHamiltonian;
	Matrix density_matrxix = singlet_projection_operator / singlet_projection_operator.diagonal().sum();

	//will put into functions later

	std::array<std::array<double, 3>, 1> oris = { {1.0,0.0,0.0} };

	auto multiply = [](std::array<double, 3> arr, double val) {for (int i = 0; i < 3; i++) { arr[i] = arr[i] * val; return arr; }};

	double t_max = 12 / KF;

	for (auto ori : oris)
	{
		auto b0vec = multiply(ori, B0);
		Matrix ZeemanHamiltionan = MakeHamiltonian(dims, 0, b0vec) + MakeHamiltonian(dims, 1, b0vec);
		Matrix Heff = H_0 + ZeemanHamiltionan;
		
		std::vector<Matrix> vec1 = { Heff, identity };
		std::vector<Matrix> vec2 = { identity, Heff.transpose() };
		Matrix Leff = std::complex<double> (0,1) * (MakeTensor(vec1, dims) - MakeTensor(vec2, dims));
		auto Leff_data = Eigen::MatrixXcd(Leff);
		auto vec = FlattenMatrixVec(density_matrxix);
		//Eigen::VectorXcd vec_eff = DotProduct(Leff, vec);
		//vec_eff = std::complex<double>(0, -1) * vec_eff;
		//std::cout << vec_eff << std::endl;
		QuantumMasterEquation eq(Leff_data);

		std::vector<state_type> x_vec = {};
		std::vector<double> time = {};
		observer ob(x_vec, time);
		std::vector<std::pair<double, std::complex<double>>> traj = { {0.0, (std::complex<double>)3.0 } };

		//int chunk_size = vec.size() / omp_get_max_threads();
		//omp_set_num_threads()
		double abs_error = 1e-8, rel_error = 1e-6;
		auto steps = odeint::integrate_const(stepper_type(), eq, vec, 0.0, t_max, 0.001, ob);
		KillThreadPool();
		auto [data, times] = ob.get_data();
		for (int i = 0; i < time.size(); i++)
		{
			Matrix rho = ReformMatrix(data[i]);
			Matrix result = singlet_projection_operator * rho;

			traj.push_back({ time[i], result.diagonal().sum() });
		}

		std::vector<double> tlist;
		std::vector<double> ps;

		std::ofstream file;
		std::string filename = "radical_pair_" + time_string + ".txt";
		file.open("radical_pair.txt");

		for (int i = 0; i < traj.size(); i++)
		{
			tlist.push_back(traj[i].first);
			ps.push_back(traj[i].second.real());
			file << tlist[i] << " , " << ps[i] << "\n";
		}
		file.close();

		double yeild = 0;

		auto f = [](double ps, double t, double kr) {return ps * std::exp(-kr * t); };
		std::vector<double> y_list = {};
		for (int i = 0; i < ps.size(); i++)
		{
			y_list.push_back(f(ps[i], tlist[i], KF));
		}
		yeild = KF * simpson_integration(tlist, y_list);

	}
	

}

