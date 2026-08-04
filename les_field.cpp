#include "les_field.h"
#include <random>
#include <fstream>
#include <cassert>
using namespace std;

Field::Field(const Mesh& m, const Config& con) :
	u(m.Nn, 0.0),
	v(m.Nn, 0.0),
	w(m.Nn, 0.0),
	p(m.Nn, 0.0),
	u_(m.Nn, 0.0),
	v_(m.Nn, 0.0),
	w_(m.Nn, 0.0),
	p_(m.Nn, 0.0),
	du(3, vector<Var>(3, Var(m.Nn, 0.0))),
	o_x(m.Nn, 0.0), o_y(m.Nn, 0.0), o_z(m.Nn, 0.0),
	o_x_(m.Nn, 0.0), o_y_(m.Nn, 0.0), o_z_(m.Nn, 0.0),
	Nu(m.Nn, 0.0), Nv(m.Nn, 0.0), Nw(m.Nn, 0.0),
	Nu_(m.Nn, 0.0), Nv_(m.Nn, 0.0), Nw_(m.Nn, 0.0),
	tau(3, vector<Var>(3, Var(m.Nn, 0.0))),
	tau_(3, vector<Var_>(3, Var_(m.Nn, 0.0))),
	dtau_xy_(m.Nn, 0.0),
	dtau_yy_(m.Nn, 0.0),
	dtau_yz_(m.Nn, 0.0),
	mesh(m),
	config(con)
{}

//初始化采用如下流场：初始流向速度的平均值沿法向为抛物形分布法
//向速度和展向速度的平均值为零初始速度场是在三
//个平均初始速度分布上附加了随机的脉动值形成的。压力的初场为零。
//考虑到LES中自然转戾的不可靠性，这也许只是暂时的，或者需要改进
void Field::Initialize()
{
	double u_m = config.Um;//初始平均流场最大无量纲速度
	double a = config.A;//叠加的随机脉动（相对）幅值
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(-1.0, 1.0);
	for (size_t i = 0; i < mesh.Nx; ++i)
	{
		for (size_t j = 0; j < mesh.Ny; ++j)
		{
			for (size_t k = 0; k < mesh.Nz; ++k)
			{
				double y = mesh.y[j];
				double u_av = config.Um * (1 - y * y);
				if (j == 0 || j == mesh.Ny - 1)
				{// 壁面无滑移
					u[mesh.idx1(i, j, k)] = 0.0;
					v[mesh.idx1(i, j, k)] = 0.0;
					w[mesh.idx1(i, j, k)] = 0.0;
				}
				else
				{
					u[mesh.idx1(i, j, k)] = u_av + config.A * u_av * dis(gen);
					v[mesh.idx1(i, j, k)] = config.A * config.Um * dis(gen);
					w[mesh.idx1(i, j, k)] = config.A * config.Um * dis(gen);
				}
			}
		}
	}
}

void Field::PhysicalToSpectral(const unique_ptr<TransformEngine>& eng, TransformSpace& ts)
{
	ts.reset();
	ts.AddPhysicalField(u);
	ts.AddPhysicalField(v);
	ts.AddPhysicalField(w);
	ts.AddPhysicalField(p);
	eng->Forward(ts);
	ts.OutSpectralField(p_);
	ts.OutSpectralField(w_);
	ts.OutSpectralField(v_);
	ts.OutSpectralField(u_);
}

void Field::SpectralToPhysical(const unique_ptr<TransformEngine>& eng, TransformSpace& ts)
{
	ts.reset();
	ts.AddSpectralField(u_);
	ts.AddSpectralField(v_);
	ts.AddSpectralField(w_);
	ts.AddSpectralField(p_);
	eng->Inverse(ts);
	ts.OutPhysicalField(p);
	ts.OutPhysicalField(w);
	ts.OutPhysicalField(v);
	ts.OutPhysicalField(u);
}

void Field::ComputeVelocityGrad(const unique_ptr<TransformEngine>& eng)
{
	int Nx = mesh.Nx, Ny = mesh.Ny, Nz = mesh.Nz, Nn = mesh.Nn, NnEx = mesh.NnEx, NyEx = mesh.NyEx;
	//对x导数用谱方法计算
	for (size_t j = 0; j < Ny; ++j)
	{
		for (size_t k = 0; k < Nz; ++k)
		{
			//从物理空间中抽取平行于x轴的一排数据用于FFT，暂存于dudx等中
			vector<complex<double>> dudx(Nx, 0.0), dvdx(Nx, 0.0), dwdx(Nx, 0.0);
			size_t jk = k * Nx * Ny + j * Nx;
			for (int I = 0; I < Nx; ++I)
			{
				dudx[I] = u[jk + I];
				dvdx[I] = v[jk + I];
				dwdx[I] = w[jk + I];
			}
			eng->FFT(&dudx[0], Nx, false, true); 
			eng->FFT(&dvdx[0], Nx, false, true);
			eng->FFT(&dwdx[0], Nx, false, true);
			for (int m = -Nx / 2; m < Nx / 2; ++m)
			{
				double md = static_cast<double>(m);
				dudx[m + Nx / 2] *= Im * mesh.Alpha * md;
				dvdx[m + Nx / 2] *= Im * mesh.Alpha * md;
				dwdx[m + Nx / 2] *= Im * mesh.Alpha * md;
			}
			eng->FFT(&dudx[0], Nx, true, true);
			eng->FFT(&dvdx[0], Nx, true, true);
			eng->FFT(&dwdx[0], Nx, true, true);
			//将求得的导数存入du中
			for (size_t i = 0; i < Nx; ++i)
			{
				size_t pos = mesh.idx1(i, j, k);
				du[0][0][pos] = real(dudx[i]);
				du[1][0][pos] = real(dvdx[i]);
				du[2][0][pos] = real(dwdx[i]);
			}
		}
	}
	//对z导数也用谱方法计算
	for (size_t i = 0; i < Nx; ++i)
	{
		for (size_t j = 0; j < Ny; ++j)
		{
			//从物理场中抽取平行于z轴的一排数据用于FFT，暂存于dudz等中
			vector<complex<double>> dudz(Nz, 0.0), dvdz(Nz, 0.0), dwdz(Nz, 0.0);
			for (size_t k = 0; k < Nz; ++k)
			{
				size_t pos = mesh.idx1(i, j, k);
				dudz[k] = u[pos];
				dvdz[k] = v[pos];
				dwdz[k] = w[pos];
			}
			eng->FFT(&dudz[0], Nz, false, true); 
			eng->FFT(&dvdz[0], Nz, false, true);
			eng->FFT(&dwdz[0], Nz, false, true);
			for (int n = -Nz / 2; n < Nz / 2; ++n)
			{
				double nd = static_cast<double>(n);
				dudz[n + Nz / 2] *= Im * mesh.Beta * nd;
				dvdz[n + Nz / 2] *= Im * mesh.Beta * nd;
				dwdz[n + Nz / 2] *= Im * mesh.Beta * nd;
			}
			eng->FFT(&dudz[0], Nz, true, true);
			eng->FFT(&dvdz[0], Nz, true, true); 
			eng->FFT(&dwdz[0], Nz, true, true);
			//将求得的导数存入du中
			for (size_t k = 0; k < Nz; ++k)
			{
				size_t pos = mesh.idx1(i, j, k);
				du[0][2][pos] = real(dudz[k]);
				du[1][2][pos] = real(dvdz[k]);
				du[2][2][pos] = real(dwdz[k]);
			}
		}
	}
	//y方向用差分方法计算
	for (size_t i = 0; i < Nx; ++i)
	{
		for (size_t j = 0; j < Ny; ++j)
		{
			for (size_t k = 0; k < Nz; ++k)
			{
				if (j == 0)
				{
					double dy = 1.0 - cos(PI * (j + 1) / (Ny - 1));
					size_t p0 = mesh.idx1(i, j, k), p1 = mesh.idx1(i, j + 1, k), p2 = mesh.idx1(i, j + 2, k);
					du[0][1][p0] = (4 * u[p1] - u[p2]) / (2 * dy);
					du[1][1][p0] = (4 * v[p1] - v[p2]) / (2 * dy);
					du[2][1][p0] = (4 * w[p1] - w[p2]) / (2 * dy);
				}
				else if (j == Ny - 1)
				{
					double dy = 1.0 - cos(PI * (j - 1) / (Ny - 1));
					size_t p = mesh.idx1(i, j, k), p_1 = mesh.idx1(i, j - 1, k), p_2 = mesh.idx1(i, j - 2, k);
					du[0][1][p] = (-4 * u[p_1] + u[p_2]) / (2 * dy);
					du[1][1][p] = (-4 * v[p_1] + v[p_2]) / (2 * dy);
					du[2][1][p] = (-4 * w[p_1] + w[p_2]) / (2 * dy);
				}
				else
				{
					double dy = cos(PI * (j - 1) / (Ny - 1)) - cos(PI * (j + 1) / (Ny - 1));
					size_t p = mesh.idx1(i, j, k), p_1 = mesh.idx1(i, j - 1, k), p1 = mesh.idx1(i, j + 1, k);
					du[0][1][p] = (u[p1] - u[p_1]) / dy;
					du[1][1][p] = (v[p1] - v[p_1]) / dy;
					du[2][1][p] = (w[p1] - w[p_1]) / dy;
				}
			}
		}
	}
}

void Field::ComputeVorticity()
{
	for (size_t i = 0; i < mesh.Nx; ++i)
	{
		for (size_t j = 0; j < mesh.Ny; ++j)
		{
			for (size_t k = 0; k < mesh.Nz; ++k)
			{
				size_t pos = mesh.idx1(i, j, k);
				o_x[pos] = du[2][1][pos] - du[1][2][pos];
				o_z[pos] = du[0][2][pos] - du[2][0][pos];
				o_y[pos] = du[1][0][pos] - du[0][1][pos];
			}
		}
	}
}

void Field::ComputeNonlinearTerms()
{
	for (size_t i = 0; i < mesh.Nx; ++i)
	{
		for (size_t j = 0; j < mesh.Ny; ++j)
		{
			for (size_t k = 0; k < mesh.Nz; ++k)
			{
				size_t pos = mesh.idx1(i, j, k);
				Nu[pos] = v[pos] * o_z[pos] - w[pos] * o_y[pos] + config.Fx;//压力驱动项也算入非线性项中
				Nv[pos] = w[pos] * o_x[pos] - u[pos] * o_z[pos];
				Nw[pos] = u[pos] * o_y[pos] - v[pos] * o_x[pos];
			}
		}
	}
}

void Field::ComputeSpectralNonlinearTermAndSGS(const unique_ptr<TransformEngine>& eng, TransformSpace& ts)
{
	//将要变换的数据（非线性项+亚格子应力）载入工作空间
	ts.reset();
	ts.AddPhysicalField(Nu);
	ts.AddPhysicalField(Nv);
	ts.AddPhysicalField(Nw);
	for (int I = 0; I < 3; ++I)
		for (int J = I; J < 3; ++J)//亚格子应力是对称张量，故只对一部分分量进行变换
			ts.AddPhysicalField(tau[I][J]);
	eng->Forward(ts);
	for (int I = 2; I >= 0; --I)
		for (int J = 2; J >= I; --J)
			ts.OutSpectralField(tau_[I][J]);
	ts.OutSpectralField(Nw_);
	ts.OutSpectralField(Nv_);
	ts.OutSpectralField(Nu_);
	//用于计算的一阶Chebyshev微分矩阵
	const Mat d1(mesh.D1);
	//涉及到对y的导数的亚格子应力分量的谱分量利用Chebyshev微分矩阵求出
	//同时求出涡量分量的谱空间投影
	for (int m = -mesh.Nx / 2; m < mesh.Nx / 2; ++m)
	{
		for (int n = -mesh.Nz / 2; n < mesh.Nz / 2; ++n)
		{
			for (int p = 0; p < mesh.Ny; ++p)
			{
				size_t wav = mesh.idx2(m, p, n);
				complex<double> w_y_(0.0), u_y_(0.0);
				for (size_t J = p + 1; J < mesh.Ny; J += 2)
				{
					size_t wav1 = mesh.idx2(m, J, n);
					double d1pJ = d1[p][J];
					dtau_xy_[wav] += d1pJ * tau_[0][1][wav1];
					dtau_yy_[wav] += d1pJ * tau_[1][1][wav1];
					dtau_yz_[wav] += d1pJ * tau_[1][2][wav1];
					w_y_ += d1pJ * w_[wav1];
					u_y_ += d1pJ * u_[wav1];
				}
				o_x_[wav] = w_y_ - Im * mesh.Beta * static_cast<double>(n) * v_[wav];
				o_y_[wav] = Im * mesh.Beta * static_cast<double>(n) * u_[wav] - Im * mesh.Alpha * static_cast<double>(m) * w_[wav];
				o_z_[wav] = Im * mesh.Alpha * static_cast<double>(m) * v_[wav] - u_y_;
			}
		}
	}
}

void Field::Output(const unique_ptr<TransformEngine>& eng, TransformSpace& ts)
{
	SpectralToPhysical(eng, ts);
	ofstream file("les_result.dat");//结果会写在文件“result”中
	assert(file.is_open());
	file << "TITLE = \"Channel Flow LES\"" << endl;
	file << "VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"P\"" << endl;
	file << "ZONE I=" << mesh.Nx << ", J=" << mesh.Ny << ", K=" << mesh.Nz << ", DATAPACKING=POINT" << endl;
	for (size_t k = 0; k < mesh.Nz; ++k)
	{
		for (size_t j = 0; j < mesh.Ny; ++j)
		{
			for (size_t i = 0; i < mesh.Nx; ++i)
			{
				double x = mesh.x[i], z = mesh.z[k], y = mesh.y[j];
				size_t pos = mesh.idx1(i, j, k);
				file << x << " " << y << " " << z << " "
					<< u[pos] << " " << v[pos] << " " << w[pos] << " " <<p[pos] << endl;
			}
		}
	}
}

void History::update()
{
	history.pop_front();
	history.push_back(next);
}