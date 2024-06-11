#pragma once

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>

#include <array>


#define G_FACTOR_ELECTRON 2.0023193043617
#define B0 1.4*2*M_PI

#define KF 1

#define KF_C 0
#define KF_D 0

#define KR_S1 1
#define KR_S2 0

#define KCD 1
#define KDC 1

#pragma region Hyperfine_mat

static double HYPERFINE_PREFACTOR = (G_FACTOR_ELECTRON * GSL_CONST_MKSA_BOHR_MAGNETON * (1 / GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR * 1e-9));

#pragma region N5_rows
static std::array<double, 3> N5r1 = { -0.0994933 * HYPERFINE_PREFACTOR, 0.00287092 * HYPERFINE_PREFACTOR, 0.0 };
static std::array<double, 3> N5r2 = { 0.00287092 * HYPERFINE_PREFACTOR, -0.0874862 * HYPERFINE_PREFACTOR, 0.0 };
static std::array<double, 3> N5r3 = { 0.0,								0.0,	1.75687 * HYPERFINE_PREFACTOR };
#pragma endregion 
static std::array<std::array<double, 3>, 3> N5 = { N5r1, N5r2, N5r3 };
#pragma region N1_rows
static std::array<double, 3> N1r1 = { -0.0529598 * HYPERFINE_PREFACTOR, -0.0586562 * HYPERFINE_PREFACTOR, 0.0460172 * HYPERFINE_PREFACTOR };
static std::array<double, 3> N1r2 = { -0.0586562 * HYPERFINE_PREFACTOR, 0.564443 * HYPERFINE_PREFACTOR, -0.564764 * HYPERFINE_PREFACTOR };
static std::array<double, 3> N1r3 = { 0.0460172 * HYPERFINE_PREFACTOR, -0.564764 * HYPERFINE_PREFACTOR, 0.453074 * HYPERFINE_PREFACTOR };
#pragma endregion
static std::array<std::array<double, 3>, 3> N1 = { N1r1, N1r2, N1r3 };
#pragma region test
//A = np.diag([-2.6, -2.6, 49.2]) * 2*math.pi
static std::array<double, 3> A_r1 = { -2.6, 0,0 };
static std::array<double, 3> A_r2 = { 0, -2.6, 0 };
static std::array<double, 3> A_r3 = { 0, 0,49.2 };
static std::array<std::array<double, 3>, 3> A = { A_r1, A_r2, A_r3 };

#pragma endregion test
#pragma endregion Hypefine_mat

//multiply hyperfine by 2pi

#pragma region FAD/WC

static std::array<double, 3> DDSE_Wc_r1 = { 26.47042689, -55.90357828, 50.1679204 };
static std::array<double, 3> DDSE_Wc_r2 = { -55.90357828, -20.86385225, 76.13493805 };
static std::array<double, 3> DDSE_Wc_r3 = { 50.1679204, 76.13493805, -5.60657464 };
static std::array<std::array<double, 3>, 3> DDSE_Wc = { DDSE_Wc_r1, DDSE_Wc_r2, DDSE_Wc_r3 };

static std::array<double, 3>N5_Wc_r1 = {-0.36082693, -0.0702137, -1.41518116};
static std::array<double, 3>N5_Wc_r2 = {-0.0702137, -0.60153649, 0.32312139};
static std::array<double, 3>N5_Wc_r3 = {-1.41518116, 0.32312139, 50.80213093};
static std::array<std::array<double, 3>, 3> N5_Wc = { N5_Wc_r1, N5_Wc_r2, N5_Wc_r3 };

static std::array<double, 3>N1_Wc_r1 = { 2.13814981, 3.19255832, -2.48895215 };
static std::array<double, 3>N1_Wc_r2 = { 3.19255832, 15.45032887, -12.44778343};
static std::array<double, 3>N1_Wc_r3 = { -2.48895215, -12.44778343, 12.49532827 };
static std::array<std::array<double, 3>, 3> N1_Wc = { N1_Wc_r1, N1_Wc_r2, N1_Wc_r3 };

#pragma endregion 

#pragma region FAD/WD

static std::array<double, 3> DDSE_Wd_r1 = { 11.08087889, -34.6687169, 12.14623706 };
static std::array<double, 3> DDSE_Wd_r2 = { -34.6687169, -33.09039672, 22.36229081 };
static std::array<double, 3> DDSE_Wd_r3 = { 12.14623706, 22.36229081, 22.00951783 };
static std::array<std::array<double, 3>, 3> DDSE_Wd = { DDSE_Wd_r1, DDSE_Wd_r2, DDSE_Wd_r3 };

static std::array<double, 3>N5_Wd_r1 = { -2.94412424e-01, -5.68059200e-02, -1.02860888e+00 };
static std::array<double, 3>N5_Wd_r2 = { -5.68059200e-02, -5.40578469e-01, -2.67686240e-02 };
static std::array<double, 3>N5_Wd_r3 = { -1.02860888e+00, -2.67686240e-02, 5.05815320e+01 };
static std::array<std::array<double, 3>, 3> N5_Wd = { N5_Wd_r1, N5_Wd_r2, N5_Wd_r3 };

static std::array<double, 3>N1_Wd_r1 = { 0.98491908, 3.28010265, -0.53784491 };
static std::array<double, 3>N1_Wd_r2 = { 3.28010265, 25.88547678, -1.6335986 };
static std::array<double, 3>N1_Wd_r3 = { -0.53784491, -1.6335986, 1.41368001 };
static std::array<std::array<double, 3>, 3> N1_Wd = { N1_Wd_r1, N1_Wd_r2, N1_Wd_r3 };

#pragma endregion

/*

FAD / Wc:

ErC_Dee = ErC_Dee = np.array([[ 26.47042689, -55.90357828, 50.1679204 ],
    [-55.90357828, -20.86385225, 76.13493805],
    [50.1679204, 76.13493805, -5.60657464]] ) # in Mrad / s

    N5 : array([[-0.36082693, -0.0702137, -1.41518116],
        [-0.0702137, -0.60153649, 0.32312139],
        [-1.41518116, 0.32312139, 50.80213093]] ) # in MHz

    N1 : array([[ 2.13814981, 3.19255832, -2.48895215],
        [3.19255832, 15.45032887, -12.44778343],
        [-2.48895215, -12.44778343, 12.49532827]] ) # in MHz


FAD / W_D:

ErD_Dee = np.array([[ 11.08087889, -34.6687169, 12.14623706],
    [-34.6687169, -33.09039672, 22.36229081],
    [12.14623706, 22.36229081, 22.00951783]] ) #  in Mrad / s

N5 : array([[-2.94412424e-01, -5.68059200e-02, -1.02860888e+00],
        [-5.68059200e-02, -5.40578469e-01, -2.67686240e-02],
        [-1.02860888e+00, -2.67686240e-02, 5.05815320e+01]] ) # in MHz

N1 : array([[ 0.98491908, 3.28010265, -0.53784491],
        [3.28010265, 25.88547678, -1.6335986],
        [-0.53784491, -1.6335986, 1.41368001]] ) # in MHz
*/