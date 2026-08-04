#include "les_types.h"
#include "les_SGSmodel.h"
#include "les_mesh.h"

using namespace std;

SmagorinskyModel::SmagorinskyModel(const Mesh& m, const Config& con) :
	mesh(m), config(con), Cs(m.Ny, 0.0), delta(m.Ny, 0.0)
{
	for (size_t j = 0; j < mesh.Ny; ++j)
	{
		double y = mesh.y[j];
		double y_plus = (y > 0 ? 1 - y : y + 1) * config.ReTau;
		//壁面阻尼系数修正的Smagorinsky系数
		Cs[j] = 0.08 * (1 - exp(-y_plus / 26));
		//槽道湍流网格长宽比较大，故用以下方法确定当量网格长度
		double dy = 0.0, dx = mesh.Lx / (mesh.Nx - 1), dz = mesh.Lz / (mesh.Nz - 1);
		if (j == 0 || j == mesh.Ny - 1)dy = 1 - cos(PI * 1 / (mesh.Ny - 1));
		else dy = 0.5 * (cos(PI * (j - 1) / (mesh.Ny - 1)) - cos(PI * (j + 1) / (mesh.Ny - 1)));
		double a1 = dx / dy, a2 = dz / dy, f = cosh(sqrt(4 * (log(a1) * log(a1) - log(a1) * log(a2) + log(a2) * log(a2)) / 27));
		//当量网格长度
		delta[j] = f * pow(dx * dy * dz, 2.0 / 3);
	}
}

void SmagorinskyModel::ComputeSGS(const vector<vector<Var>>& du, vector<vector<Var>>& tau) const
{
	for (size_t i = 0; i < mesh.Nx; ++i)
	{
		for (size_t j = 0; j < mesh.Ny; ++j)
		{
			for (size_t k = 0; k < mesh.Nz; ++k)
			{
				size_t pos = mesh.idx1(i, j, k);
				//计算S_ij*S_ij
				double S = 0.0;
				for (int I = 0; I < 3; ++I)
				{
					for (int J = 0; J < 3; ++J)
					{
						double S_IJ = 0.5 * (du[I][J][pos] + du[J][I][pos]);
						S += S_IJ * S_IJ;
					}
				}
				S = sqrt(S);
				//亚格子涡粘系数
				double mu = Cs[j] * Cs[j] * delta[j] * S;
				//计算亚格子应力分量
				for (int I = 0; I < 3; ++I)
				{
					for (int J = 0; J <= I; ++J)
					{
						double S_IJ = du[I][J][pos] + du[J][I][pos];
						tau[I][J][pos] = mu * S_IJ;
						if (I != J)tau[J][I][pos] = tau[I][J][pos];
					}
				}
			}
		}
	}
}