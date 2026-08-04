#pragma once
#include <vector>

class Mesh
{
public:
	Mesh(int nx, int ny, int nz, double h, double lx, double lz);
	//网格维度
	const int Nx, Ny, Nz;//各方向网格点总数
	const int Nn;//物理空间网格点总数
	const int NyEx;//Chebyshev变换拓展后的y方向网格点数
	const int NnEx;//Chebyshev变换拓展后物理空间总网格点数
	//域尺寸
	const double H;//计算域高度
	const double Lx, Lz;//流向、展向域长度
	const double Alpha, Beta;//基本波数
	//预计算的网格坐标
	const std::vector<double> x;
	const std::vector<double> y;
	const std::vector<double> z;
	//一阶Ny* Ny大小的Chebyshev微分矩阵
	Mat D1;
	//二阶Ny* Ny大小的Chebyshev微分矩阵
	Mat D2;
	//索引映射
	inline size_t idx1(size_t i, size_t j, size_t k) const//物理空间坐标到一维扁平索引
	{
		return k * Nx * Ny + j * Nx + i;
	}
	inline size_t idx2(int m, int p, int n) const//谱空间坐标到一维扁平索引
	{
		return (m + Nx / 2) * Ny * Nz + (n + Nz / 2) * Ny + p;
	}
};