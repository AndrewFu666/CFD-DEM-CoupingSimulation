#include "les_types.h"
#include "les_mesh.h"
#include <vector>

using namespace std;

static vector<double> GenerateUniformGrid(size_t n, double L)
{
	vector<double> g(n, 0.0);
	for (size_t i = 0; i < n; ++i)
		g[i] = L * i / (n - 1);
	return g;
}

static vector<double> GaussLobbatoPoints(size_t n)
{
	vector<double> g(n, 0.0);
	for (size_t j = 0; j < n; ++j)
		g[j] = -cos(PI * j / (n - 1));
	return g;
}

Mesh::Mesh(int nx, int ny, int nz, double h, double lx, double lz):
    Nx(nx), Ny(ny), Nz(nz),
    Nn(nx* ny* nz),
    NyEx(2 * ny - 2),
    NnEx(nx* nz* NyEx),
    H(h), Lx(lx), Lz(lz),
    Alpha(2 * PI / lx), Beta(2 * PI / lz),
    x(GenerateUniformGrid(nx, lx)),
    y(GaussLobbatoPoints(ny)),
    z(GenerateUniformGrid(nz, lz)),
	D2(Ny, vector<double>(Ny, 0.0)), 
	D1(Ny, vector<double>(Ny, 0.0))
{
	//计算一阶Ny*Ny大小的Chebyshev微分矩阵
	for (size_t I = 0; I < Ny; ++I)
	{
		for (size_t J = I + 1; J < Ny; J += 2)
		{
			double c = (I == 0) ? 2 : 1;//存疑
			D1[I][J] = 2.0 * J / c;
		}
	}
	//计算二阶Ny*Ny大小的Chebyshev微分矩阵，算到倒数第三行，剩下两行用边界条件方程替换
	for (size_t I = 0; I < Ny - 2; I++)
	{
		for (size_t J = I; J < Ny; J++)
		{
			double s = 0.0;
			for (size_t K = I + 1; K < J; K++)
			{
				s += D1[I][K] * D1[K][J];
			}
			D2[I][J] = s;
		}
	}
	//以下初始化仅适用于边界条件为Neumann边界条件的情形，若是Dirichlet边界条件则须重新赋值
	//二阶Chebyshev微分矩阵的倒数第二行用y=1处的边界条件导出的方程代替
	for (int p = 0; p < Ny; ++p)
	{
		double s = 0.0;
		for (size_t I = 0; I < p; ++I)
			s += D1[I][p];
		D2[Ny - 2][p] = s;
	}
	//最后一行用y=-1处的边界条件导出的方程代替
	for (int p = 0; p < Ny; ++p)
	{
		double s = 0.0;
		int k = 1;
		for (size_t I = 0; I < p; ++I)
		{
			s += D1[I][p] * k;
			k = -k;
		}
		D2[Ny - 1][p] = s;
	}
}