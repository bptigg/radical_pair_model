#include "Structure.h"

static std::map<std::string, std::array<std::array<double, 3>, 3>*> s_HyperfineList = {
    {"N5", &N5},
    {"N1", &N1},
    {"N5_Wc", &N5_Wc},
    {"N1_Wc", &N1_Wc},
    {"N5_Wd", &N5_Wd},
    {"N1_Wd", &N1_Wd},
    {"A", &A}

};

static std::map<std::string, std::array<std::array<double, 3>, 3>*> s_EEDList = {
    {"EED_Wd", &DDSE_Wd},
    {"EED_Wc", &DDSE_Wc}
};

Matrix Structure::SingleRadicalPair(Matrix& singlet_projection_operator)
{
    Matrix identity_mat = MakeSpinOperator(m_dims, {});
    int Liouville_dim = std::pow(m_dim, 2);

    Matrix identity_supor = identity(std::pow(m_dim, 2));

    Matrix Leff(Liouville_dim, Liouville_dim);
    {
        std::vector<Matrix> vec1 = { m_Hamiltonians[0], identity_mat };
        std::vector<Matrix> vec2 = { identity_mat, m_Hamiltonians[0].transpose() };
        std::vector<Matrix> vec3 = { singlet_projection_operator, identity_mat };
        std::vector<Matrix> vec4 = { identity_mat, singlet_projection_operator };
        Leff = (std::complex<double>(0, 1) * (MakeTensor(vec1, m_dims) - MakeTensor(vec2, m_dims))) + (0.5 * KR_S1 * (MakeTensor(vec3, m_dims) + MakeTensor(vec4, m_dims)) + (KF * identity_supor));
    }

    return Leff;
}
std::vector<Matrix> Structure::DoubleRadicalPair(Matrix& singlet_projection_operator)
{
    Matrix identity_mat = MakeSpinOperator(m_dims, {});
    int Liouville_dim = std::pow(m_dim, 2);

    Matrix identity_supor = identity(std::pow(m_dim, 2));

    Matrix Leff_c(Liouville_dim, Liouville_dim);
    {
        std::vector<Matrix> vec1 = { m_Hamiltonians[0], identity_mat };
        std::vector<Matrix> vec2 = { identity_mat, m_Hamiltonians[0].transpose() };
        std::vector<Matrix> vec3 = { singlet_projection_operator, identity_mat };
        std::vector<Matrix> vec4 = { identity_mat, singlet_projection_operator };
        Leff_c = (std::complex<double>(0, 1) * (MakeTensor(vec1, m_dims) - MakeTensor(vec2, m_dims))) + (0.5 * KR_S1 * (MakeTensor(vec3, m_dims) + MakeTensor(vec4, m_dims))) + (KF_C * identity_supor);
    }

    Matrix Leff_d(Liouville_dim, Liouville_dim);
    {
        std::vector<Matrix> vec1 = { m_Hamiltonians[1], identity_mat };
        std::vector<Matrix> vec2 = { identity_mat, m_Hamiltonians[1].transpose() };
        std::vector<Matrix> vec3 = { singlet_projection_operator, identity_mat };
        std::vector<Matrix> vec4 = { identity_mat, singlet_projection_operator };
        Leff_d = (std::complex<double>(0, 1) * (MakeTensor(vec1, m_dims) - MakeTensor(vec2, m_dims))) + (0.5 * KR_S2 * (MakeTensor(vec3, m_dims) + MakeTensor(vec4, m_dims))) + (KF_D * identity_supor);
    }

    int dim_size = 2 * Liouville_dim;
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
            if (sys.mode == 2 and (i == 1 || i == 2))
            {
                continue;
            }
            for (int j = 0; j < Leff_size; j++)
            {
                for (int k = 0; k < Leff_size; k++)
                {
                    if (sys.mode == 2 and (j == k))
                    {
                        continue;
                    }

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
        
        if (sys.mode != 2)
        {
            return { Leff_data };
        }

        entries.clear();
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < Leff_size; j++)
            {
                for (int k = 0; k < Leff_size; k++)
                {
                    if (sys.mode == 2 and (i == 0 || i == 3))
                    {
                        if (j != k)
                        {
                            continue;
                        }
                    }

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

        Matrix ExchangeMatrix(dim_size, dim_size);
        ExchangeMatrix.setFromTriplets(entries.begin(), entries.end());
        return { Leff_data , ExchangeMatrix };
    }
}
void Structure::UpdateZeemanHamiltonian(std::array<double, 3> b0_orientaion)
{
    m_Hamiltonians.clear();
    Matrix ZeemanHamiltonian = MakeHamiltonian(m_dims, 0, b0_orientaion) + MakeHamiltonian(m_dims, 1, b0_orientaion);
    for (auto H : m_OrientationFreeHamiltonians)
    {
        Matrix temp = H + ZeemanHamiltonian;
        m_Hamiltonians.push_back(temp);
    }
}

void Structure::CreateRadicalSystem()
{
    m_OrientationFreeHamiltonians.clear();

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
        };
    };

    std::vector<std::pair<Matrix, int>> eed;

    for (auto i : m_DipoleDipoleRadicalBinding)
    {
        eed.push_back({MakeHamiltonian(m_dims, 0, 1, *s_EEDList[i.first]), i.second});
    }

    std::vector<std::pair<Matrix, int>> hyp;

    for (auto i : m_HyperfineNucleiBinding)
    {
        auto [nucleus, binding, rp] = i;

        std::array<std::array<double, 3>, 3> hyperfine = *s_HyperfineList[nucleus];
        for (int i = 0; i < 9; i++)
        {
            int row = std::floor((double)i / 3.0), col = i % 3;
            hyperfine[row][col] = hyperfine[row][col] * 2 * M_PI;
        }

        hyp.push_back({ MakeHamiltonian(m_dims, binding.first, binding.second, hyperfine), rp });
    }


    for (int i = 0; i < m_NumRadicals; i++)
    {
        auto second = [i](const auto& val) {
            return val.second == i;
            };
        
        Matrix H_0(m_dim, m_dim);
        auto it = std::find_if(eed.begin(), eed.end(), second);
        H_0 = H_0 + it->first;

        for (int e = 0; e < m_HyperfineNucleiBinding.size(); e++)
        {
            if (hyp[e].second == i)
            {
                H_0 = H_0 + hyp[e].first;
            }
        }

        m_OrientationFreeHamiltonians.push_back(H_0);
    }
}

std::vector<Matrix> Structure::CreateSuperOperator(Matrix& singlet_projection_operator)
{
    switch (m_NumRadicals)
    {
    case 1:
        return { SingleRadicalPair(singlet_projection_operator) };
        break;
    case 2:
        return DoubleRadicalPair(singlet_projection_operator);
        break;
    default:
        break;
    }
}

Structure::Structure(Structure_param param, system_setup system)
    :sys(system)
{
    m_NumRadicals = param.num_radicals;
    m_DipoleDipoleRadicalBinding = param.DipoleBinding;
    m_dims = { 2,2 };
    for (int i = 0; i < param.spins.size(); i++)
    {
        m_dims.push_back(std::round((2 * param.spins[i]) + 1.0));
    }
    for (int i = 0; i < param.HyperfineBinding.size(); i++)
    {
        auto [paramatrix, binding, radical] = param.HyperfineBinding[i];
        m_HyperfineNucleiBinding.push_back(std::make_tuple(paramatrix, binding, radical));
    }

    auto product = [](int a, int b) {return a * b; };
    m_dim = std::reduce(m_dims.begin(), m_dims.end(), 1, product);
}

Structure::~Structure()
{
    m_Hamiltonians.clear();
    m_OrientationFreeHamiltonians.clear();
}
