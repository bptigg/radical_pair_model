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
typedef odeint::runge_kutta_dopri5<state_type> stepper_type;


typedef std::vector<std::tuple<double, std::complex<double>, std::complex<double>>> TrajectoryType;


class QuantumMasterEquation
{
	//Eigen::MatrixXcd m_supor;
	Matrix m_supor;

public:
	QuantumMasterEquation(Matrix& SuperOperator)
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

class output_lock_guard
{
private:
	bool m_kill = false;
	int m_write_limit = 0;

	std::mutex OutputLock;
public:
	const bool& kill = m_kill;
	const int& write_limit = m_write_limit;
public:
	void UpdateKill(bool input)
	{
		OutputLock.lock();
		m_kill = input;
		OutputLock.unlock();
	}

	void UpdateWriteLimit(int limit)
	{
		OutputLock.lock();
		if (limit < 0)
		{
			limit = 0;
		}
		m_write_limit = limit;
		OutputLock.unlock();

	}
};

struct observer
{
	std::vector<state_type>& m_states;
	std::vector<double>& m_times;

	output_lock_guard* m_olg;
	uint32_t m_write_limit;

	observer(std::vector<state_type>& states, std::vector<double>& times, output_lock_guard* olg)
		:m_states(states), m_times(times), m_olg(olg)
	{
		m_write_limit = 0;
	}

	void operator() (state_type& rho, double t)
	{
 		m_states.push_back(rho);
		m_times.push_back(t);

		m_write_limit++;

		if (m_write_limit % 10 == 0)
		{
			std::cout << "Current time step: " << m_times[m_times.size() - 1] << std::endl;
		}

		m_olg->UpdateWriteLimit(m_write_limit);

	}

	std::tuple<std::vector<state_type>, std::vector<double>> get_data()
	{
		return std::make_tuple(m_states, m_times);
	}
};

class Simulation
{
private:
	TrajectoryType* m_traj;
	output_lock_guard* m_olg;
	observer* m_obv;

	Matrix* Ps;
public:
	Simulation(TrajectoryType* traj, output_lock_guard* olg, observer* obv, Matrix* ps)
		:m_traj(traj), m_olg(olg), m_obv(obv), Ps(ps) {}
	
	void state_fraction()
	{
		uint32_t loop = 0;
		while (!m_olg->kill)
		{
			while (loop < m_olg->write_limit)
			{
				state_type result_all = m_obv->m_states[0];
				const size_t vec_size = result_all.size() / 2;

				state_type result_1(result_all.begin(), result_all.begin() + vec_size);
				state_type result_2(result_all.begin() + vec_size, result_all.end());

				double time = m_obv->m_times[0];

				m_obv->m_states.erase(m_obv->m_states.begin());
				m_obv->m_times.erase(m_obv->m_times.begin());

				Matrix rho_1 = ReformMatrix(result_1);
				Matrix result_mat_1 = *Ps * rho_1;

				Matrix rho_2 = ReformMatrix(result_2);
				Matrix result_mat_2 = *Ps * rho_2;

				m_traj->push_back(std::make_tuple(time, result_mat_1.diagonal().sum(), result_mat_2.diagonal().sum()));
				loop++;
			}
		}
	}
};

int main() 
{
	std::time_t t = std::time(0);
	std::tm now = *std::localtime(&t);
	char buf[20];
	strftime(buf, sizeof(buf), "%H_%M_%S__%d_%m_%y", &now);
	std::string time_string = buf;

	//auto DDSE = PointDipoleDipoleCoupling({ 8.51061, -14.251621, 6.5492562 });
	
	auto DDSE_c = DDSE_Wc;
	auto DDSE_d = DDSE_Wd;
	
	//std::vector<int32_t>dims = { 2,2,3,3,3};
	std::vector<int32_t>dims = { 2,2,2,2,2 };

	auto product = [](int a, int b) {return a * b; };
	uint32_t dim = std::reduce(dims.begin(), dims.end(), 1, product);

	MATRIX3x3 ide;

	for (int i = 0; i < 3; i++)
	{
		for (int e = 0; e < 3; e++)
		{
			DDSE_c[i][e] = DDSE_c[i][e] * 2 * M_PI;
			DDSE_d[i][e] = DDSE_d[i][e] * 2 * M_PI;

			if (i == e)
			{
				ide[i][e] = 1;
				continue;
			}
			ide[i][e] = 0;
		}

	}

	Matrix singlet_projection_operator = 0.25 * MakeSpinOperator(dims, {}) - MakeHamiltonian(dims, 0, 1, ide);
	Matrix identity_mat = MakeSpinOperator(dims, {});
	Matrix triplet_projection_operator = identity_mat - singlet_projection_operator;

	Matrix H_0_C(dim, dim);
	{
		Matrix DipoleDipoleElectronHamiltonian = MakeHamiltonian(dims, 0, 1, DDSE_c);

		Matrix H1 = MakeHamiltonian(dims, 0, 2, N5_Wc);
		Matrix H2 = MakeHamiltonian(dims, 1, 3, N1_Wc);
		Matrix HyperfineHamiltonian = H1 + H2;

		H_0_C = HyperfineHamiltonian + DipoleDipoleElectronHamiltonian;
	}

	Matrix H_0_D(dim, dim);
	{
		Matrix DipoleDipoleElectronHamiltonian = MakeHamiltonian(dims, 0, 1, DDSE_d);

		Matrix H1 = MakeHamiltonian(dims, 0, 2, N5_Wd);
		Matrix H2 = MakeHamiltonian(dims, 1, 4, N1_Wd);
		Matrix HyperfineHamiltonian = H1 + H2;

		H_0_D = HyperfineHamiltonian + DipoleDipoleElectronHamiltonian;
	}

	state_type vec = {};
	
	{
		Matrix rho_c = singlet_projection_operator / singlet_projection_operator.diagonal().sum();
		Matrix rho_d = MakeZeroOperator(dims);

		auto vec_c = FlattenMatrixVec(rho_c);
		auto vec_d = FlattenMatrixVec(rho_d);

		vec.insert(vec.end(), std::make_move_iterator(vec_c.begin()), std::make_move_iterator(vec_c.end()));
		vec.insert(vec.end(), std::make_move_iterator(vec_d.begin()), std::make_move_iterator(vec_d.end()));
	}

	//will put into functions later

	std::array<std::array<double, 3>, 1> oris = { {1.0,0.0,0.0} };

	auto multiply = [](std::array<double, 3> arr, double val) {for (int i = 0; i < 3; i++) { arr[i] = arr[i] * val; return arr; }};

	//double t_max = 1 / KF;
	double t_max = 1;

	output_lock_guard* olg = new output_lock_guard;

	for (auto ori : oris)
	{
		olg->UpdateKill(false);

		auto b0vec = multiply(ori, B0);
		Matrix ZeemanHamiltionan = MakeHamiltonian(dims, 0, b0vec) + MakeHamiltonian(dims, 1, b0vec);
		
		Matrix Heff_c = H_0_C + ZeemanHamiltionan;
		Matrix Heff_d = H_0_D + ZeemanHamiltionan;
		
		int dim2 = std::pow(dim, 2);
		Matrix Leff_c(dim2, dim2);
		{
			std::vector<Matrix> vec1_c = { Heff_c, identity_mat };
			std::vector<Matrix> vec2_c = { identity_mat, Heff_c.transpose() };
			std::vector<Matrix> vec3_c = { singlet_projection_operator, identity_mat };
			std::vector<Matrix> vec4_c = { identity_mat, singlet_projection_operator };
			Leff_c = (std::complex<double>(0, 1) * (MakeTensor(vec1_c, dims) - MakeTensor(vec2_c, dims))) + (0.5 * KR_S1 * (MakeTensor(vec3_c, dims) + MakeTensor(vec4_c, dims)));
		}
		//auto Leff_data = Eigen::MatrixXcd(Leff_c);
		Matrix Leff_d(dim2, dim2);
		Matrix identity_supor = identity(std::pow(dim, 2));
		{
			std::vector<Matrix> vec1_d = { Heff_d, identity_mat };
			std::vector<Matrix> vec2_d = { identity_mat, Heff_d.transpose() };
			Leff_d = (std::complex<double>(0, 1) * (MakeTensor(vec1_d, dims) - MakeTensor(vec2_d, dims))) + (KF * identity_supor);
		}

		int dim_size = 2 * std::pow(dim, 2);
		Matrix Leff_data(dim_size, dim_size);
		{
			
			Matrix q1 = Leff_c + KCD * identity_supor;
			Matrix q2 = -1 * KDC * identity_supor;
			Matrix q3 = -1 * KCD * identity_supor;
			Matrix q4 = Leff_d + KDC * identity_supor;

			std::vector<Matrix*> mat_corner = { &q1, &q2, &q3, &q4 };
			int Leff_size = dim_size * 0.5;
			
			typedef Eigen::Triplet<std::complex<double>, int32_t> T;
			std::vector<T> entries;

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < Leff_size; j++)
				{
					for (int k = 0; k < Leff_size; k++)
					{
						std::complex<double> value = mat_corner[i]->coeff(j, k);
						if (value == std::complex<double>(0.0, 0.0))
						{ 
							continue;
						}
						std::pair<int, int> pos = { (std::floor((double)i / 2.0) * Leff_size) + j, ((i) % 2 * Leff_size) + k };
						entries.push_back(T(pos.first, pos.second, value));
					}
				}
			}

			Leff_data.setFromTriplets(entries.begin(), entries.end());

		}

		QuantumMasterEquation eq(Leff_data);

		std::vector<state_type> x_vec = {};
		std::vector<double> time = {};
		observer ob(x_vec, time, olg);

		TrajectoryType traj = {std::make_tuple(0.0, std::complex<double>(1.0,0), std::complex<double>(0.0,0.0))};
		Simulation* sim = new Simulation(&traj, olg, &ob, &singlet_projection_operator);

		std::thread sim_thread(&Simulation::state_fraction, sim);

		double abs_error = 1e-8, rel_error = 1e-6;
		auto steps = odeint::integrate_const(stepper_type(), eq, vec, 0.0, t_max, 0.0001, ob);
		KillThreadPool();

		olg->UpdateKill(true);
		sim_thread.join();
		delete sim;

		auto [data, times] = ob.get_data();
		for (int i = 0; i < time.size(); i++)
		{
			state_type result_all = data[i];
			const size_t vec_size = result_all.size() / 2;

			state_type result_1(result_all.begin(), result_all.begin() + vec_size);
			state_type result_2(result_all.begin() + vec_size, result_all.end());

			Matrix rho_1 = ReformMatrix(result_1);
			Matrix result_mat_1 = singlet_projection_operator * rho_1;
			Matrix rho_2 = ReformMatrix(result_1);
			Matrix result_mat_2 = singlet_projection_operator * rho_2;

			traj.push_back(std::make_tuple(time[i], result_mat_1.diagonal().sum(), result_mat_2.diagonal().sum()));
		}

		std::vector<double> tlist;
		std::vector<std::vector<double>> ps;

		std::ofstream file;
		std::string filename = "radical_pair_" + time_string + ".txt";
		file.open(filename);

		for (int i = 0; i < traj.size(); i++)
		{
			auto [t, tr1, tr2] = traj[i];
			tlist.push_back(t);
			ps.push_back({ tr1.real(), tr2.real() });
			file << tlist[i] << " , " << ps[i][0] << " , " << ps[i][1] << "\n";
		}
		file.close();

		double yeild = 0;

		auto f = [](double ps, double t, double kr) {return ps * std::exp(-kr * t); };
		
		std::vector<double> y_list = {};
		for (int i = 0; i < ps.size(); i++)
		{
			y_list.push_back(f(ps[i][0], tlist[i], KF));
		}
		yeild = KF * simpson_integration(tlist, y_list);

	}

	delete olg;
	

}

