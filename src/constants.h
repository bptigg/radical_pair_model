#pragma once

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>

#include <array>


#define G_FACTOR_ELECTRON 2.0023193043617


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

#pragma endregion Hypefine_mat
