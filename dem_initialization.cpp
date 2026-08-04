#include "dem_initialization.h"
#include <random>

using namespace std;

void initializeParticles(const CFDData& cfdData, Particles& p,const InitBox& box, unsigned seed)
{
	random_device rd;
	mt19937 gen(seed == 0 ? rd() : seed);
	double radius = 0.5 * p.diameter;
	//从物理流场统计速度分布
	size_t totalNodes = cfdData.nx * cfdData.ny * cfdData.nz;
	double sum = 0.0, sum2 = 0.0;
	for (size_t i = 0; i < totalNodes; ++i)
	{
		sum += cfdData.u[i];
		sum2 += cfdData.u[i] * cfdData.u[i];
	}
	double meanVx = sum / totalNodes;
	double stdVx = sqrt(max(sum2 / totalNodes - meanVx * meanVx, 0.0));
	//位置：均匀分布，保证球心在 [radius, lx-radius] 等范围内
	uniform_real_distribution<double> distX(box.xMin + radius, box.xMax - radius);
	uniform_real_distribution<double> distY(box.yMin + radius, box.yMax - radius);
	uniform_real_distribution<double> distZ(box.zMin + radius, box.zMax - radius);
	//速度分布（正态分布，均值和标准差）
	normal_distribution<double> distVx(meanVx, stdVx);
	normal_distribution<double> distVy(0.0, stdVx * 0.6); // 横向扰动较小
	normal_distribution<double> distVz(0.0, stdVx * 0.6);
	for (size_t i = 0; i < p.count; ++i)
	{
		bool overlap = true;
		while (overlap)
		{
			p.x[i] = distX(gen);
			p.y[i] = distY(gen);
			p.z[i] = distZ(gen);
			//检查与已生成颗粒是否重叠
			overlap = false;
			for (size_t j = 0; j < i; ++j)
			{
				double dx = p.x[i] - p.x[j], dy = p.y[i] - p.y[j], dz = p.z[i] - p.z[j];
				if (dx * dx + dy * dy + dz * dz < p.diameter * p.diameter)
				{
					overlap = true;
					break;
				}
			}
		}
		p.vx[i] = distVx(gen);
		p.vy[i] = distVy(gen);
		p.vz[i] = distVz(gen);
		p.omegaX[i] = p.omegaY[i] = p.omegaZ[i] = 0.0;
	}
}