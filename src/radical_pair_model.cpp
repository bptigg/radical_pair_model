#include <Eigen/Core>
#include <random>
#include <iostream>
#include <tuple>
#include <fstream>
#include <ctime>
#include <chrono>
#include <string>
#include <filesystem>


//#include "spin_operator.h"
#include "Methods.h"
#include "constants.h"
#include "Structure.h"

#include <Eigen/Dense>
#include <omp.h>
#include <boost\numeric\odeint.hpp>
#include <boost\numeric\odeint\external\eigen\eigen.hpp>
#include <boost\numeric\odeint\external\openmp\openmp.hpp>
#include <Eigen/Eigenvalues>

//#include <Eigen\SuperLUSupport>

#pragma warning( disable : 4996)


typedef std::vector<std::complex<double>> state_type;
using namespace boost::numeric;
typedef odeint::runge_kutta_dopri5<state_type> stepper_type;


#define MODE 1


//typedef std::vector<std::tuple<double, std::complex<double>, std::complex<double>>> TrajectoryType;
typedef std::vector<std::pair<double, std::vector<std::complex<double>>>> TrajectoryType;

class QuantumMasterEquation
{
	//Eigen::MatrixXcd m_supor;
	Matrix m_supor;
	Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor>* m_cached_matrix;
	int m_dims;

public:
	QuantumMasterEquation(Matrix& SuperOperator)
		:m_supor(SuperOperator)
	{
		m_cached_matrix = nullptr;
		m_dims = 0;
	}

	QuantumMasterEquation(Matrix& SuperOperator, int dims)
		:m_supor(SuperOperator)
	{
		m_cached_matrix = nullptr;
		m_dims = dims;
	}

	~QuantumMasterEquation()
	{
		if (m_cached_matrix != nullptr)
		{
			delete m_cached_matrix;
		}
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

	double SparseSolver(const state_type& rho, const Matrix& projection_operator, system_setup& setup)
	{
		Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>> Solver;

		int n = 16;
		omp_set_num_threads(n);
		Eigen::setNbThreads(n);

		//Eigen::BiCGSTAB<Matrix> Solver2;
		//
		//int leff_size = m_dims * 0.5;
		//Matrix m_supor_preconditioner(m_dims, m_dims);
		//Matrix top_left = m_supor.topLeftCorner(leff_size, leff_size);
		//Matrix bottom_right = m_supor.bottomRightCorner(leff_size, leff_size);
		//
		//std::vector<Matrix*> mat_corner = { &top_left, &bottom_right };
		//
		//typedef Eigen::Triplet<std::complex<double>, int32_t> T;
		//std::vector<T> entries;
		//
		//for (int i = 0; i < 2; i++)
		//{
		//	for (int j = 0; j < leff_size; j++)
		//	{
		//		//int ind_size = mat_corner[i]->row(j).coeffs().size();
		//		//std::cout << ind_size << std::endl;
		//		for (int k = 0; k < leff_size; k++)
		//		{
		//			std::complex<double> value = mat_corner[i]->coeff(j, k);
		//			if (value == std::complex<double>(0.0, 0.0))
		//			{
		//				continue;
		//			}
		//			std::pair<int, int> pos = { (i * leff_size) + j, (i * leff_size) + k };
		//			entries.push_back(T(pos.first, pos.second, value));
		//		}
		//	}
		//}
		//m_supor_preconditioner.setFromTriplets(entries.begin(), entries.end());

		//Solver.compute(m_supor_preconditioner);
		Solver.compute(m_supor);
		if (Solver.info() != Eigen::Success)
		{
			return 0;
		}

		//Eigen::VectorXcd rho_naught_2(m_dims);
		Eigen::SparseVector<std::complex<double>> rho_naught(m_dims);
		
		for (int i = 0; i < m_dims; i++)
		{
			if (i >= rho.size())
			{
				//rho_naught.coeffRef(i) = std::complex<double>(0, 0);
				break;
			}
			else if (rho[i] != std::complex<double>(0, 0))
			{
				rho_naught.coeffRef(i) = rho[i];
			}

			//guess[i] = rho_naught[i] / m_supor.coeff(i, i);
		}

		//Eigen::VectorXcd guess(m_dims);
		//guess = Solver.solve(rho_naught);
		//
		//Solver2.compute(m_supor);
		//rho_naught_2 = Eigen::VectorXcd(rho_naught);
		//Eigen::VectorXcd observable(m_dims);
		//Solver2.setTolerance(1e-16);
		//double error = 1;
		//observable = Solver2.solveWithGuess(rho_naught_2, guess);

		Eigen::SparseVector<std::complex<double>> observable(m_dims);
		observable = Solver.solve(rho_naught);
		std::vector<double> yeilds = {};
		for (int i = 0; i < setup.num_radicals; i++)
		{
			state_type vec = {};
			for (int e = 0; e < setup.num_radicals; e++)
			{
				if (e != i)
				{
					state_type vec_temp = FlattenMatrixVec(MakeZeroOperator(setup.dims));
					vec.insert(vec.end(), vec_temp.begin(), vec_temp.end());
				}
				else
				{
					state_type projection_vec = FlattenMatrixVec(projection_operator);
					vec.insert(vec.end(), projection_vec.begin(), projection_vec.end());
				}
			}

			std::complex<double> sum(0, 0);
			for (int a = 0; a < vec.size(); a++)
			{
				sum = sum + (setup.rate_constants[i] * vec[a] * observable.coeff(a));
			}
			double yeild = sum.real();
			yeilds.push_back(yeild);
		}

		double total_yeild = 0.0;

		for (int i = 0; i < yeilds.size(); i++)
		{
			total_yeild = total_yeild + yeilds[i];
		}
		return total_yeild;
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
	int m_NumRadicals;
public:
	Simulation(TrajectoryType* traj, output_lock_guard* olg, observer* obv, Matrix* ps, int num_radicals)
		:m_traj(traj), m_olg(olg), m_obv(obv), Ps(ps), m_NumRadicals(num_radicals) {}
	
	void state_fraction()
	{
		uint32_t loop = 0;
		while (!m_olg->kill)
		{
			while (loop < m_olg->write_limit)
			{
				state_type result_all = m_obv->m_states[0];
				const size_t vec_size = result_all.size() / m_NumRadicals;

				std::vector<std::complex<double>> Trace;

				for (int i = 0; i < m_NumRadicals; i++)
				{
					state_type temp(result_all.begin() + (i * vec_size), result_all.begin() + ((i + 1) * vec_size));
					Matrix rho = ReformMatrix(temp);
					Matrix ResultMat = *Ps * rho;
					Trace.push_back(ResultMat.diagonal().sum());
				}

				double time = m_obv->m_times[0];

				m_obv->m_states.erase(m_obv->m_states.begin());
				m_obv->m_times.erase(m_obv->m_times.begin());

				//m_traj->push_back(std::make_tuple(time, result_mat_1.diagonal().sum(), result_mat_2.diagonal().sum()));
				m_traj->push_back({ time, Trace });
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

	Structure_param two_radical;
	two_radical.spins = {0.5, 0.5, 0.5};
	two_radical.num_radicals = 2;
	two_radical.DipoleBinding = { {"EED_Wc", 0}, {"EED_Wd", 1} };
	two_radical.HyperfineBinding = { {"N5_Wc", {0,2}, 0}, {"N1_Wc", {1,3}, 0}, {"N5_Wd", {0,2}, 1}, {"N1_Wc", {1,4}, 1} };

	Structure_param one_radical;
	one_radical.spins = { 1, 1 };
	//one_radical.spins = { 1 };
	one_radical.num_radicals = 1;
	one_radical.DipoleBinding = { {"EED_Wc", 0} };
	one_radical.HyperfineBinding = { {"N5_Wc", {0,2}, 0}, {"N1_Wc", {1,3}, 0} };
	//one_radical.HyperfineBinding = { {"A", {0,2}, 0} };

	Structure radical_sys(two_radical);
	radical_sys.CreateRadicalSystem();

	MATRIX3x3 ide;
	
	for (int i = 0; i < 3; i++)
	{
		for (int e = 0; e < 3; e++)
		{
			if (i == e)
			{
				ide[i][e] = 1;
				continue;
			}
			ide[i][e] = 0;
		}
	
	}
	
	auto dims = radical_sys.get_dims();
	Matrix singlet_projection_operator = 0.25 * MakeSpinOperator(dims, {}) - MakeHamiltonian(dims, 0, 1, ide);
	Matrix identity_mat = MakeSpinOperator(dims, {});
	Matrix triplet_projection_operator = identity_mat - singlet_projection_operator;

	state_type vec = {};
	
	{
		Matrix rho_c = singlet_projection_operator / singlet_projection_operator.diagonal().sum();
		//Matrix rho_d = MakeZeroOperator(dims);

		auto vec_c = FlattenMatrixVec(rho_c);
		//auto vec_d = FlattenMatrixVec(rho_d);

		vec.insert(vec.end(), std::make_move_iterator(vec_c.begin()), std::make_move_iterator(vec_c.end()));
		//vec.insert(vec.end(), std::make_move_iterator(vec_d.begin()), std::make_move_iterator(vec_d.end()));
	}

	std::vector<std::array<double, 3>> oris = {};
	std::vector<double> phi;
	double theta = M_PI_2;
	std::vector<std::array<double, 3>> yeilds;
	for (int i = 0; i < 180; i++)
	{
		phi.push_back(((M_PI_2 / 90) * i));
	}
	for (auto p : phi)
	{
		double x = std::sin(p) * std::cos(theta);
		double y = std::sin(p) * std::sin(theta);
		double z = std::cos(p);
		oris.push_back({ x,y,z });
		yeilds.push_back({ p,theta,0 });
	}
	//std::vector<std::array<double, 3>> oris = FibonacciSphere(250);
	auto multiply = [](std::array<double, 3> arr, double val) {for (int i = 0; i < 3; i++) { arr[i] = arr[i] * val; return arr; }};

	//double t_max = 1 / KF;
	double t_max = 1;

	{
		std::filesystem::path filepath = std::string("Results/test.txt");
		bool exists = std::filesystem::is_directory(filepath.parent_path());
		if (!exists)
		{
			std::filesystem::create_directory("Results");
		}

		filepath = std::string("Results/yeild/test.txt");
		exists = std::filesystem::is_directory(filepath.parent_path());
		if (!exists)
		{
			std::filesystem::create_directory("Results/yeild");
		}
	}

	output_lock_guard* olg = new output_lock_guard;
	int i = 0;
	for (auto ori : oris)
	{;
		olg->UpdateKill(false);
		int dim_size = radical_sys.get_dims(true)[0];
		Matrix Leff_data(dim_size, dim_size);

		{
			auto b0vec = multiply(ori, B0);
			radical_sys.UpdateZeemanHamiltonian(b0vec);
			//Matrix Leff_data = radical_sys.CreateSuperOperator(singlet_projection_operator);
			Leff_data = radical_sys.CreateSuperOperator(singlet_projection_operator);
			radical_sys.clear();
#if MODE == 0:
			DotProductVec(Leff_data, {}, true);
#endif
		}

		//Matrix m(1, 1);
		//QuantumMasterEquation eq(m);
		system_setup setup = { radical_sys.get_num_radicals(), {KF_C, KF_D}, radical_sys.get_dims()};
		QuantumMasterEquation eq(Leff_data, Leff_data.rows());
		double y = eq.SparseSolver(vec, singlet_projection_operator, setup);
		yeilds[i][2] = y;

		/*
		std::vector<state_type> x_vec = {};
		std::vector<double> time = {};
		observer ob(x_vec, time, olg);

		TrajectoryType traj = { {0.0, {std::complex<double>(1.0,0)}} };
		//TrajectoryType traj = {std::make_tuple(0.0, std::complex<double>(1.0,0), std::complex<double>(0.0,0.0))};
		Simulation* sim = new Simulation(&traj, olg, &ob, &singlet_projection_operator, radical_sys.get_num_radicals());

		std::thread sim_thread(&Simulation::state_fraction, sim);

		double abs_error = 1e-8, rel_error = 1e-6;
		auto steps = odeint::integrate_const(stepper_type(), eq, vec, 0.0, t_max, 0.001, ob);
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

			//traj.push_back(std::make_tuple(time[i], result_mat_1.diagonal().sum(), result_mat_2.diagonal().sum()));
		}

		std::vector<double> tlist;
		std::vector<std::vector<double>> ps;
		std::vector<std::string> FileVec;

		for (int i = 0; i < radical_sys.get_num_radicals(); i++)
		{
			std::string filename = "radical_pair_" + std::to_string((i + 1)) + "_" + time_string + ".txt";
			FileVec.push_back(filename);
		}

		std::fstream file;
		
		for (int e = 0; e < traj.size(); e++)
		{
			auto [t, tr] = traj[e];
			tlist.push_back(t);

			std::vector<double> RealTrace;
			for (auto trace : tr)
			{
				RealTrace.push_back(trace.real());
			}

			for (int i = 0; i < radical_sys.get_num_radicals(); i++)
			{
				file.open(FileVec[i], std::ios::app);
				file << tlist[e] << " , " << RealTrace[i] << "\n";
				file.close();
			}
		}

		double yeild = 0;

		auto f = [](double ps, double t, double kr) {return ps * std::exp(-kr * t); };
		
		std::vector<double> y_list = {};
		for (int i = 0; i < ps.size(); i++)
		{
			y_list.push_back(f(ps[i][1], tlist[i], KF)); //only looking at the singlet fraction of the second radical pair due to the current scheme 
		}
		yeild = KF * simpson_integration(tlist, y_list);
		*/
		i++;
	}
#if MODE == 1:
	std::string filename = "radical_pair_yeild_only_" + time_string + ".txt";
	std::fstream file;
	file.open(filename, std::ios::app);
	for (int i = 0; i < yeilds.size(); i++)
	{
		file << yeilds[i][0] << " , " << yeilds[i][1] << " , " << yeilds[i][2] << "\n";
	}
	file.close();

#endif

	delete olg;
	

}

