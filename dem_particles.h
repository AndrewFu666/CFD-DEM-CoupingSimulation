#pragma once
#include <vector>

class Particles
{
public:
	Particles(double d, double m, size_t n, double COR, double COF) :
		diameter(d), mass(m), e(COR), mu(COF), count(n),
		vx(n, 0.0), vy(n, 0.0), vz(n, 0.0),
		omegaX(n, 0.0), omegaY(n, 0.0), omegaZ(n, 0.0),
		x(n, 0.0), y(n, 0.0), z(n, 0.0) { }
	const double diameter;//颗粒直径
	const double mass;//颗粒质量
	const double e;//碰撞恢复系数
	const double mu;//颗粒摩擦系数
	const size_t count;//颗粒数量
	std::vector<double> vx, vy, vz;//颗粒速度分量
	std::vector<double> omegaX, omegaY, omegaZ;//颗粒转动速度分量
	std::vector<double> x, y, z;//颗粒位置
	void output(double solutionTime);
};