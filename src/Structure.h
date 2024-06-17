#pragma once
#include <vector>
#include <string>
#include <map>

#include "constants.h"
#include "Methods.h"

struct Structure_param
{
	int num_radicals;
	std::vector<double> spins;
	std::vector<std::pair<std::string, int>> DipoleBinding;
	std::vector<std::tuple<std::string, std::pair<int, int>, int>> HyperfineBinding;
};

struct system_setup
{
	int num_radicals = 0;
	std::vector<double> rate_constants = {};
	std::vector<int> dims = {};
};

class Structure
{
private:
	
	std::vector<int> m_dims;
	int m_dim;
	
	std::vector<std::tuple<std::string, std::pair<int, int>, int>> m_HyperfineNucleiBinding;
	std::vector<std::pair<std::string, int>> m_DipoleDipoleRadicalBinding;
	int m_NumRadicals;

	std::vector<Matrix> m_OrientationFreeHamiltonians;
	std::vector<Matrix> m_Hamiltonians;

private:
	Matrix SingleRadicalPair(Matrix& singlet_projection_operator);
	Matrix DoubleRadicalPair(Matrix& singlet_projection_operator);

public:
	void UpdateZeemanHamiltonian(std::array<double, 3> b0_orientaion);
	void CreateRadicalSystem();
	Matrix CreateSuperOperator(Matrix& singlet_projection_operator);

	Structure(Structure_param param);

	~Structure();

	std::vector<int> get_dims(bool val = false) {
		if (val)
		{
			return { m_dim };
		}
		return m_dims;
	}
	int get_num_radicals() {
		return m_NumRadicals;
	}

	void clear() {
		m_Hamiltonians.clear();
	}

};

