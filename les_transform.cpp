#include "les_types.h"
#include "les_transform.h"
#include <cassert>
#include <complex>

using namespace std;

void TransformSpace::AddPhysicalField(const Var& var)//用于将待正向变换的数据装载入变换空间
{
	for (int I = 0; I < mesh.Nn; ++I)yz[Nvar * mesh.Nn + I] = var[I];
	++Nvar;
}

void TransformSpace::AddSpectralField(const Var_& var_)//用于将待逆向变换的数据装载入变换空间
{
	for (int m = -mesh.Nx / 2; m < mesh.Nx / 2; ++m)
	{
		for (int n = -mesh.Nz / 2; n < mesh.Nz / 2; ++n)
		{
			size_t dest = Nvar * mesh.NnEx + (m + mesh.Nx / 2) * mesh.NyEx * mesh.Nz + (n + mesh.Nz / 2) * mesh.NyEx;
			for (int p = 0; p < mesh.Ny; ++p)
				spec[dest + p] = var_[(m + mesh.Nx / 2) * mesh.Ny * mesh.Nz + (n + mesh.Nz / 2) * mesh.Ny + p];
		}
	}
	++Nvar;
}

void TransformSpace::OutPhysicalField(Var& var)//用于将已逆向变换的数据输出
{
	assert(Nvar != 0);
	--Nvar;
	for (int I = 0; I < mesh.Nn; ++I)var[I] = real(yz[Nvar * mesh.Nn + I]);
}

void TransformSpace::OutSpectralField(Var_& var_)//用于将已正向变换的数据输出
{
	assert(Nvar != 0);
	--Nvar;
	for (int m = -mesh.Nx / 2; m < mesh.Nx / 2; ++m)
	{
		for (int n = -mesh.Nz / 2; n < mesh.Nz / 2; ++n)
		{
			size_t src = Nvar * mesh.NnEx + (m + mesh.Nx / 2) * mesh.NyEx * mesh.Nz + (n + mesh.Nz / 2) * mesh.NyEx;
			for (int p = 0; p < mesh.Ny; ++p)
				var_[(m + mesh.Nx / 2) * mesh.Ny * mesh.Nz + (n + mesh.Nz / 2) * mesh.Ny + p] = spec[src + p];
		}
	}
}

void MyEngine::FFT(complex<double>* begin, size_t N, bool invert, bool center)
{//正向FFT函数接受复数数组，返回所有谱分量，invert为true时进行逆向FFT
	assert(begin != nullptr);
	assert((N & (N - 1)) == 0); // 必须是2的幂
	//频移：频率范围从0~N变为-N/2~N/2-1
	if (!invert && center)
	{
		for (size_t n = 0; n < N; ++n)
			if (n & 1) begin[n] = -begin[n]; // 乘以 (-1)^n（ebeginp(iPIn)）
	}
	//位逆序重排
	unsigned num = log2(N);
	for (size_t i = 0; i < N; ++i)
	{
		size_t idx = 0;
		size_t a = i;
		for (int j = num - 1; j >= 0; --j)
		{
			idx += (1ULL << j) * (a & 1);
			a >>= 1;
		}
		if (i < idx) swap(begin[i], begin[idx]);
	}
	// 迭代蝶形
	for (size_t len = 2; len <= N; len <<= 1)
	{
		double ang = (invert ? 2 : -2) * PI / len;
		complex<double> domega(cos(ang), sin(ang));
		for (size_t i = 0; i < N; i += len)
		{
			complex<double> omega(1);
			for (size_t j = 0; j < len / 2; ++j)
			{
				complex<double> u = begin[i + j];
				complex<double> v = begin[i + j + len / 2] * omega;
				begin[i + j] = u + v;
				begin[i + j + len / 2] = u - v;
				omega *= domega;
			}
		}
	}
	if (invert)
	{
		for (size_t i = 0; i < N; ++i) begin[i] /= N;
		if (center)
		{
			for (size_t n = 0; n < N; ++n)
				if (n & 1) begin[n] = -begin[n];
		}
	}
}

void MyEngine::Forward(TransformSpace& ts)
{
	int Nx = mesh.Nx, Ny = mesh.Ny, Nz = mesh.Nz, Nn = mesh.Nn, NnEx = mesh.NnEx, NyEx = mesh.NyEx;
	//先在x方向进行FFT，得到中间结果，存储在工作空间中
	for (size_t j = 0; j < Ny; ++j)
	{
		for (size_t k = 0; k < Nz; ++k)
		{
			size_t jk = k * Nx * Ny + j * Nx;
			for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
			{
				size_t head = cnt * Nn;//第cnt个物理量在工作空间中的存储位置
				FFT(&ts.yz[head + jk], Nx, false, true);
			}
		}
	}
	//再在z方向进行FFT，得到进一步结果，存储在工作空间中
	for (size_t j = 0; j < Ny; ++j)
	{
		for (int m = -Nx / 2; m < Nx / 2; ++m)
		{
			//从第一步结果中抽取出平行于z轴的一排数据用于FFT
			size_t jm = j * Nx * Nz + (m + Nx / 2) * Nz;
			for (size_t k = 0; k < Nz; ++k)
			{
				size_t idx_yz = k * Nx * Ny + j * Nx + m + Nx / 2;
				size_t idx_y = jm + k;
				for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
					ts.y[cnt * Nn + idx_y] = ts.yz[cnt * Nn + idx_yz];
			}
			for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
			{
				size_t head = cnt * Nn;//第cnt个物理量在工作空间中的存储位置
				FFT(&ts.y[head + jm], Nz, false, true);
			}
		}
	}
	//对y方向进行变换要先延拓
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int n = -Nz / 2; n < Nz / 2; ++n)
		{
			size_t mn = (n + Nz / 2) * Nx * NyEx + (m + Nx / 2) * NyEx;
			//从第二步结果中抽取出平行于y轴的一列数据用于FFT
			for (size_t j = 0; j < Ny; ++j)
			{
				size_t idx_y = j * Nx * Nz + (m + Nx / 2) * Nz + n + Nz / 2;
				for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
					ts.spec[cnt * NnEx + mn + j] = ts.y[cnt * Nn + idx_y];
			}
			//延拓
			for (size_t j = 1; j < Ny - 1; ++j)
			{
				for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
				{
					auto temp = ts.spec[cnt * NnEx + mn + Ny - 1 - j];
					ts.spec[cnt * NnEx + mn + Ny - 1 + j] = temp;
				}
			}
			for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
			{
				FFT(&ts.spec[cnt * NnEx + mn], NyEx, false, false);
				for (int p = 0; p < Ny; ++p)ts.spec[cnt * NnEx + mn + p] *= 1.0 / (Ny - 1);
			}
		}
	}
}

void MyEngine::Inverse(TransformSpace& ts)
{
	int Nx = mesh.Nx, Ny = mesh.Ny, Nz = mesh.Nz, Nn = mesh.Nn, NnEx = mesh.NnEx, NyEx = mesh.NyEx;
	//先在p方向进行IDCT，得到中间结果，存储在工作空间中
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int n = -Nz / 2; n < Nz / 2; ++n)
		{
			size_t mn = (n + Nz / 2) * Nx * NyEx + (m + Nx / 2) * NyEx;
			//延拓
			for (int p = 1; p < Ny - 1; ++p)
			{
				for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
				{
					auto temp = ts.spec[cnt * NnEx + mn + Ny - 1 - p];
					ts.spec[cnt * NnEx + mn + Ny - 1 + p] = temp;
				}
			}
			for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
			{
				FFT(&ts.spec[cnt * NnEx + mn], NyEx, true, false);
				for (size_t j = 0; j < Ny; ++j)ts.spec[cnt * NnEx + mn + j] *= Ny - 1;
			}
		}
	}
	//再在n方向进行IFFT，得到进一步结果
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int j = 0; j < Ny; ++j)
		{//从第一步结果中抽取出平行于n轴的一排数据用于IFFT
			size_t jm = j * Nx * Nz + (m + Nx / 2) * Nz;
			for (int n = -Nz / 2; n < Nz / 2; ++n)
			{
				size_t idx_y = jm + n + Nz / 2;
				size_t idx_spec = (n + Nz / 2) * Nx * NyEx + (m + Nx / 2) * NyEx + j;
				for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
					ts.y[cnt * Nn + idx_y] = ts.spec[cnt * NnEx + idx_spec];
			}
			for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
				FFT(&ts.y[cnt * Nn + jm], Nz, true, true);
		}
	}
	//最后对m方向进行IFFT
	for (int j = 0; j < Ny; ++j)
	{
		for (size_t k = 0; k < Nz; ++k)
		{
			//从第二步结果中抽取出平行于m轴的一列数据用于IFFT
			size_t jk = k * Nx * Ny + j * Nx;
			for (int m = -Nx / 2; m < Nx / 2; ++m)
			{
				size_t idx_yz = jk + m + Nx / 2;
				size_t idx_y = j * Nx * Nz + (m + Nx / 2) * Nz + k;
				for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
					ts.yz[cnt * Nn + idx_yz] = ts.y[cnt * Nn + idx_y];
			}
			for (size_t cnt = 0; cnt < ts.Nvar; ++cnt)
				FFT(&ts.yz[cnt * Nn + jk], Nx, true, true);
		}
	}
}