#include "dem_time_integrator.h"
#include "dem_vector_operate.h"

using namespace std;

static Vec3 interpolationMethod(int i, int j, int k, const CFDData& f)
{
	const int nx = f.nx, ny = f.ny, nz = f.nz;
	auto idx = [nx, ny, nz](int ii, int jj, int kk) -> int { return kk * nx * ny + jj * nx + ii; };
	double uAtp = 0.125 * (f.u[idx(i, j, k)] + f.u[idx(i + 1, j, k)] + f.u[idx(i, j + 1, k)] + f.u[idx(i, j, k + 1)] +
		f.u[idx(i + 1, j + 1, k)] + f.u[idx(i + 1, j, k + 1)] + f.u[idx(i, j + 1, k + 1)] + f.u[idx(i + 1, j + 1, k + 1)]);
	double vAtp = 0.125 * (f.v[idx(i, j, k)] + f.v[idx(i + 1, j, k)] + f.v[idx(i, j + 1, k)] + f.v[idx(i, j, k + 1)] +
		f.v[idx(i + 1, j + 1, k)] + f.v[idx(i + 1, j, k + 1)] + f.v[idx(i, j + 1, k + 1)] + f.v[idx(i + 1, j + 1, k + 1)]);
	double wAtp = 0.125 * (f.w[idx(i, j, k)] + f.w[idx(i + 1, j, k)] + f.w[idx(i, j + 1, k)] + f.w[idx(i, j, k + 1)] +
		f.w[idx(i + 1, j + 1, k)] + f.w[idx(i + 1, j, k + 1)] + f.w[idx(i, j + 1, k + 1)] + f.w[idx(i + 1, j + 1, k + 1)]);
	return { uAtp,vAtp,wAtp };
}

static Vec3 velocityInterpolate(int pID, const CFDData& f, const Particles& p)
{
	const double dx = f.lx / (f.nx - 1), dz = f.lz / (f.nz - 1);//真实的x、z方向网格长度
	//计算与颗粒相邻的（较小）网格节点的空间索引
	int i = static_cast<int>(p.x[pID] / dx), k = static_cast<int>(p.z[pID] / dz);
	const double yCFD = p.y[pID] - 0.5 * f.h;//CFD和DEM坐标系的x-z平面不重合，需要作变换
	int j = lower_bound(f.y.begin(), f.y.end(), yCFD) - f.y.begin() - 1;
	//速度插值
	return interpolationMethod(i, j, k, f);
}

static Vec3 computeDragForce(int pID, const CFDData& f, const Particles& p)
{
	Vec3 vf = velocityInterpolate(pID, f, p);//颗粒所在位置的流体速度
	Vec3 vp = { p.vx[pID],p.vy[pID],p.vz[pID] };//颗粒速度
	Vec3 dv = { vf.x - vp.x,vf.y - vp.y,vf.z - vp.z };//速度差
	double dvMagnitude = len(dv);
	//计算曳力
	constexpr double PI = 3.14159265358979323846;
	double Re_p = p.diameter * dvMagnitude / f.viscosity;//颗粒雷诺数（流体体积分数近似为1）
	double Cd = (Re_p < 1000) ? 24 * (1 + 0.15 * pow(Re_p, 0.687)) / Re_p : 0.44;
	double k = 0.125 * Cd * PI * f.density * p.diameter * p.diameter * dvMagnitude;//曳力公式标量因子
	return{ k * dv.x,k * dv.y,k * dv.z };
}

void integrate(const CFDData& f, Particles& p, double dt_DEM)
{
	double m = p.mass;
	for (int i = 0; i < p.count; ++i)
	{
		//更新位置（显式方法）
		p.x[i] += p.vx[i] * dt_DEM;
		p.y[i] += p.vy[i] * dt_DEM;
		p.z[i] += p.vz[i] * dt_DEM;
		//更新速度（显式方法）
		Vec3 dragForce = computeDragForce(i, f, p);//计算曳力
		double ax = dragForce.x / m, ay = dragForce.y / m, az = dragForce.z / m;//计算加速度（仅考虑曳力）
		p.vx[i] += ax * dt_DEM;
		p.vy[i] += ay * dt_DEM;
		p.vz[i] += az * dt_DEM;
	}
}