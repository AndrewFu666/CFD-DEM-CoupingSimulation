#pragma once
#include <vector>
#include <complex>

using Var = std::vector<double>;
using Var_ = std::vector<std::complex<double>>;
using Mat = std::vector<std::vector<double>>;

constexpr double PI = 3.14159265358979323846;
constexpr std::complex<double> Im(0, 1.0);