#include "les_time_splitting_solver.h"
#include "les_field.h"
#include "les_SGSmodel.h" 

using namespace std;

void NonlinearStep(History& data,const unique_ptr<TransformEngine>& eng, TransformSpace& ts, const unique_ptr<SGSModel>& sgs)
{
	auto& s = data.s(), & s_1 = data.s_1(), & s_2 = data.s_2(), & next = data.next;
	int Nx = data.mesh.Nx, Ny = data.mesh.Ny, Nz = data.mesh.Nz;
	double Alpha = data.mesh.Alpha, Beta = data.mesh.Beta, dt = data.config.dt, q = data.config.q;
	//先更新s时间步物理空间的数据
	s.SpectralToPhysical(eng, ts);
	//计算速度梯度
	s.ComputeVelocityGrad(eng);
	//先计算s时间步上的涡量，后面要用
	s.ComputeVorticity();
	//计算s时间步物理空间非线性项
	s.ComputeNonlinearTerms();
	//计算亚格子应力
	sgs->ComputeSGS(s.du, s.tau);
	//计算s时间步亚格子应力和非线性项在谱空间的投影
	s.ComputeSpectralNonlinearTermAndSGS(eng, ts);
	//时间推进到s+1/3
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int n = -Nz / 2; n < Nz / 2; ++n)
		{
			for (int p = 0; p < Ny; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				double md = static_cast<double>(m), nd = static_cast<double>(n);
				complex<double> uRHSs = s.Nu_[wav] + q * (Im * Alpha * md * s.tau_[0][0][wav] + s.dtau_xy_[wav] + Im * Beta * nd * s.tau_[0][2][wav]),
					uRHSs_1 = s_1.Nu_[wav] + q * (Im * Alpha * md * s_1.tau_[0][0][wav] + s_1.dtau_xy_[wav] + Im * Beta * nd * s_1.tau_[0][2][wav]),
					uRHSs_2 = s_2.Nu_[wav] + q * (Im * Alpha * md * s_2.tau_[0][0][wav] + s_2.dtau_xy_[wav] + Im * Beta * nd * s_2.tau_[0][2][wav]);
				next.u_[wav] = 3.0 * s.u_[wav] - 1.5 * s_1.u_[wav] + 1.0 * s_2.u_[wav] / 3.0 + dt * (3.0 * uRHSs - 3.0 * uRHSs_1 + uRHSs_2);
				complex<double> vRHSs = s.Nv_[wav] + q * (Im * Alpha * md * s.tau_[0][1][wav] + s.dtau_yy_[wav] + Im * Beta * nd * s.tau_[1][2][wav]),
					vRHSs_1 = s_1.Nv_[wav] + q * (Im * Alpha * md * s_1.tau_[0][1][wav] + s_1.dtau_yy_[wav] + Im * Beta * nd * s_1.tau_[1][2][wav]),
					vRHSs_2 = s_2.Nv_[wav] + q * (Im * Alpha * md * s_2.tau_[0][1][wav] + s_2.dtau_yy_[wav] + Im * Beta * nd * s_2.tau_[1][2][wav]);
				next.v_[wav] = 3.0 * s.v_[wav] - 1.5 * s_1.v_[wav] + 1.0 * s_2.v_[wav] / 3.0 + dt * (3.0 * vRHSs - 3.0 * vRHSs_1 + vRHSs_2);
				complex<double> wRHSs = s.Nw_[wav] + q * (Im * Alpha * md * s.tau_[0][2][wav] + s.dtau_yz_[wav] + Im * Beta * nd * s.tau_[2][2][wav]),
					wRHSs_1 = s_1.Nw_[wav] + q * (Im * Alpha * md * s_1.tau_[0][2][wav] + s_1.dtau_yz_[wav] + Im * Beta * nd * s_1.tau_[2][2][wav]),
					wRHSs_2 = s_2.Nw_[wav] + q * (Im * Alpha * md * s_2.tau_[0][2][wav] + s_2.dtau_yz_[wav] + Im * Beta * nd * s_2.tau_[2][2][wav]);
				next.w_[wav] = 3.0 * s.w_[wav] - 1.5 * s_1.w_[wav] + 1.0 * s_2.w_[wav] / 3.0 + dt * (3.0 * wRHSs - 3.0 * wRHSs_1 + wRHSs_2);
			}
		}
	}
}

void Helmholtz(Mat A, vector<complex<double>>& b, complex<double>* x, const Mesh& m)
{
	int Ny = m.Ny;
	//列主元高斯消元
	for (int i = 0; i < Ny; ++i)
	{
		//选主元：只需比较最后两行的同列元素，其他中间位置上的元素都是0
		int pivotRow = i;
		double maxAbs = abs(A[i][i]);
		if (abs(A[Ny - 2][i]) > maxAbs && i < Ny - 1)
		{
			maxAbs = abs(A[Ny - 2][i]);
			pivotRow = Ny - 2;
		}
		if (abs(A[Ny - 1][i]) > maxAbs)
		{
			maxAbs = abs(A[Ny - 1][i]);
			pivotRow = Ny - 1;
		}
		// 交换行
		if (pivotRow != i)
		{
			swap(A[i], A[pivotRow]);
			swap(b[i], b[pivotRow]);
		}
		//消去最后两行同列元素
		int start = Ny - 2;
		if (i == Ny - 2 || i == Ny - 1)start = i + 1;
		for (int k = start; k < Ny; ++k)
		{
			double factor = A[k][i] / A[i][i];
			for (int j = i; j < Ny; ++j)
				A[k][j] -= factor * A[i][j];
			b[k] -= factor * b[i];
		}
	}
	// 回代求解（此时A为上三角）
	for (int i = Ny - 1; i >= 0; --i)
	{
		complex<double> s = b[i];
		for (int j = i + 1; j < Ny; ++j)
			s -= A[i][j] * x[j];
		x[i] = s / A[i][i];
	}
}

void PressureStep(History& data)
{
	auto& s = data.s(), & s_1 = data.s_1(), & s_2 = data.s_2(), & next = data.next;
	int Nx = data.mesh.Nx, Ny = data.mesh.Ny, Nz = data.mesh.Nz, Nn = data.mesh.Nn;
	double Alpha = data.mesh.Alpha, Beta = data.mesh.Beta, dt = data.config.dt, Re = data.config.ReTau;
	//构建压力方程系数矩阵，用二阶Chebyshev微分矩阵初始化
	Mat A(data.mesh.D2);
	//用于计算的一阶Chebyshev微分矩阵
	Mat d1(data.mesh.D1);
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int n = -Nz / 2; n < Nz / 2; ++n)
		{
			if (m == 0 && n == 0)//m和n都为0时直接将p_置为0，若按原先方法计算会导致系数矩阵奇异
			{
				size_t mn = Nn / 2 + Nz * Ny / 2;
				for (int p = 0; p < Ny; ++p)next.p_[mn + p] = 0;
				continue;
			}
			double md = static_cast<double>(m), nd = static_cast<double>(n);
			//对于每一对m、n系数矩阵对角项为相应的波数
			double a_pp = -pow(Alpha * m, 2) - pow(Beta * n, 2);
			for (int p = 0; p < Ny - 2; ++p)
				A[p][p] = a_pp;
			//构建压力方程的右端项
			vector<complex<double>> RHS(Ny, 0.0);
			//计算前Ny-1个方程右端项，为此需要先计算s+1/3时间步上v对y的导数在谱空间的投影
			vector<complex<double>> dvdy_(Ny, 0.0);
			//利用一阶Chebyshev矩阵计算v对y的导数的谱分量
			for (int p = 0; p < Ny; ++p)
			{
				for (size_t J = p + 1; J < Ny; J += 2)
				{
					size_t wav = data.mesh.idx2(m, J, n);
					dvdy_[p] += d1[p][J] * next.v_[wav];
				}
			}
			for (int p = 0; p < Ny - 2; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				complex<double> r_p = (Im * Alpha * md * next.u_[wav] + dvdy_[p] + Im * Beta * nd * next.w_[wav]) / dt;
				RHS[p] = r_p;
			}
			//边界条件方程的右端项由边界条件确定
			vector<complex<double>> G_up(3, complex<double>(0.0)), G_low(3, complex<double>(0.0));//上/下边界对应的右端项
			//计算边界条件对应的右端项
			double k = 1.0;
			for (int p = 0; p < Ny; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				G_up[0] += Im * Alpha * md * s.o_x_[wav] - Im * Beta * nd * s.o_z_[wav];
				G_low[0] += (Im * Alpha * md * s.o_x_[wav] - Im * Beta * nd * s.o_z_[wav]) * k;
				G_up[1] += Im * Alpha * md * s_1.o_x_[wav] - Im * Beta * nd * s_1.o_z_[wav];
				G_low[1] += (Im * Alpha * md * s_1.o_x_[wav] - Im * Beta * nd * s_1.o_z_[wav]) * k;
				G_up[2] += Im * Alpha * md * s_2.o_x_[wav] - Im * Beta * nd * s_2.o_z_[wav];
				G_low[2] += (Im * Alpha * md * s_2.o_x_[wav] - Im * Beta * nd * s_2.o_z_[wav]) * k;
				k = -k;
			}
			for (int I = 0; I < 3; ++I)
			{
				G_up[I] = G_up[I] / Re;
				G_low[I] = G_low[I] / Re;
			}
			RHS[Ny - 2] = 3.0 * G_up[0] - 3.0 * G_up[1] + G_up[2];
			RHS[Ny - 1] = 3.0 * G_low[0] - 3.0 * G_low[1] + G_low[2];
			//求解方程组
			size_t mn = (m + Nx / 2) * Ny * Nz + (n + Nz / 2) * Ny;
			Helmholtz(A, RHS, &next.p_[mn], data.mesh);
			//推进到s+2/3时间步
			for (int p = 0; p < Ny; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				next.u_[wav] += -Im * Alpha * md * next.p_[wav] * dt;
				next.w_[wav] += -Im * Beta * nd * next.p_[wav] * dt;
				//对速度分量v推进需要计算П^(1)
				complex<double> P1 = { 0,0 };
				for (size_t J = p + 1; J < Ny; J += 2) P1 += d1[p][J] * next.p_[data.mesh.idx2(m, J, n)];
				next.v_[wav] += -P1 * dt;
			}
		}
	}
}

inline void BuildCoefficientMatrixOfVelocity(Mat& A, const Mesh& m)
{
	int Ny = m.Ny;
	//系数矩阵的倒数第二行用y=1处的边界条件导出的方程代替
	for (size_t j = 0; j < Ny; ++j)
		A[Ny - 2][j] = 1.0;
	A[Ny - 2][0] *= 0.5; A[Ny - 2][Ny - 1] *= 0.5;
	//系数矩阵的最后一行用y=-1处的边界条件导出的方程代替
	int a = 1;
	for (size_t j = 0; j < Ny; ++j)
	{
		A[Ny - 1][j] = a;
		a = -a;
	}
	A[Ny - 1][0] *= 0.5; A[Ny - 1][Ny - 1] *= 0.5;
}

void ViscosityStep(History& data)
{
	auto& s = data.s(), & next = data.next;
	Mat Au(data.mesh.D2);
	BuildCoefficientMatrixOfVelocity(Au, data.mesh);
	int Nx = data.mesh.Nx, Ny = data.mesh.Ny, Nz = data.mesh.Nz;
	double Alpha = data.mesh.Alpha, Beta = data.mesh.Beta, dt = data.config.dt, Re = data.config.ReTau, q = data.config.q;
	double gamma = 11.0 / 6;
	//粘性步需要解三个速度方程组，其形式与压力方程非常相似，三个方向速度分量方程的系数矩阵完全一致
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int n = -Nz / 2; n < Nz / 2; ++n)
		{
			//对于每一对m、n设定系数矩阵对角项
			double a_pp = -pow(Alpha * m, 2) - pow(Beta * n, 2) - Re / (gamma * dt);
			for (size_t p = 0; p < Ny - 2; ++p)
				Au[p][p] = a_pp;
			Mat Av(Au), Aw(Au);
			//构建压力方程的右端项
			vector<complex<double>> ru(Ny, 0.0), rv(Ny, 0.0), rw(Ny, 0.0);
			for (int p = 0; p < Ny - 2; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				double md = static_cast<double>(m), nd = static_cast<double>(n);
				complex<double> sgs_u = (1 - q) * (Im * Alpha * md * s.tau_[0][0][wav] + s.dtau_xy_[wav] + Im * Beta * nd * s.tau_[0][2][wav]) / gamma,
					sgs_v = (1 - q) * (Im * Alpha * md * s.tau_[0][1][wav] + s.dtau_yy_[wav] + Im * Beta * nd * s.tau_[1][2][wav]) / gamma,
					sgs_w = (1 - q) * (Im * Alpha * md * s.tau_[0][2][wav] + s.dtau_yz_[wav] + Im * Beta * nd * s.tau_[2][2][wav]) / gamma;
				ru[p] = -Re * next.u_[wav] / dt + sgs_u;
				rv[p] = -Re * next.v_[wav] / dt + sgs_v;
				rw[p] = -Re * next.w_[wav] / dt + sgs_w;
			}
			//上下壁面无滑移条件（不用作任何处理）
			//推进到下一子时间步
			size_t mn = (m + Nx / 2) * Ny * Nz + (n + Nz / 2) * Ny;
			Helmholtz(Au, ru, &next.u_[mn], data.mesh);
			Helmholtz(Av, rv, &next.v_[mn], data.mesh);
			Helmholtz(Aw, rw, &next.w_[mn], data.mesh);
		}
	}
}

//以下是前两次时间推进的程序
void FirstStep(History& data, const unique_ptr<TransformEngine>& eng, TransformSpace& ts, const unique_ptr<SGSModel>& sgs)
{
	auto& s_2 = data.s_2(), & s_1 = data.s_1();
	int Nx = data.mesh.Nx, Ny = data.mesh.Ny, Nz = data.mesh.Nz, Nn = data.mesh.Nn;
	double Alpha = data.mesh.Alpha, Beta = data.mesh.Beta, dt = data.config.dt, Re = data.config.ReTau, q = data.config.q;
	double gamma = 11.0 / 6;
	//************************非线性步************************
	//先更新0时间步物理空间的数据
	s_2.PhysicalToSpectral(eng, ts);
	//计算速度梯度
	s_2.ComputeVelocityGrad(eng);
	//先计算0时间步上的涡量，后面要用
	s_2.ComputeVorticity();
	//计算0时间步物理空间非线性项
	s_2.ComputeNonlinearTerms();
	//计算亚格子应力
	sgs->ComputeSGS(s_2.du, s_2.tau);
	//计算0时间步亚格子应力和非线性项在谱空间的投影
	s_2.ComputeSpectralNonlinearTermAndSGS(eng, ts);
	//时间推进到0+1/3
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int n = -Nz / 2; n < Nz / 2; ++n)
		{
			for (int p = 0; p < Ny; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				double md = static_cast<double>(m), nd = static_cast<double>(n);
				complex<double> uRHS = s_2.Nu_[wav] + q * (Im * Alpha * md * s_2.tau_[0][0][wav] + s_2.dtau_xy_[wav] + Im * Beta * nd * s_2.tau_[0][2][wav]);
				s_1.u_[wav] = s_2.u_[wav] + dt * uRHS;
				complex<double> vRHS = s_2.Nv_[wav] + q * (Im * Alpha * md * s_2.tau_[0][1][wav] + s_2.dtau_yy_[wav] + Im * Beta * nd * s_2.tau_[1][2][wav]);
				s_1.v_[wav] = s_2.v_[wav] + dt * vRHS;
				complex<double> wRHS = s_2.Nw_[wav] + q * (Im * Alpha * md * s_2.tau_[0][2][wav] + s_2.dtau_yz_[wav] + Im * Beta * nd * s_2.tau_[2][2][wav]);
				s_1.w_[wav] = s_2.w_[wav] + dt * wRHS;
			}
		}
	}
	//**************************压强步************************
	//构建压力方程系数矩阵，用二阶Chebyshev微分矩阵初始化
	Mat A(data.mesh.D2);
	//用于计算的一阶Chebyshev微分矩阵
	Mat d1(data.mesh.D1);
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int n = -Nz / 2; n < Nz / 2; ++n)
		{
			if (m == 0 && n == 0)
			{
				size_t mn = Nn / 2 + Nz * Ny / 2;
				for (int p = 0; p < Ny; ++p)s_1.p_[mn + p] = 0;
				continue;
			}
			double md = static_cast<double>(m), nd = static_cast<double>(n);
			//对于每一对m、n系数矩阵对角项为相应的波数
			double a_pp = -pow(Alpha * m, 2) - pow(Beta * n, 2);
			for (int p = 0; p < Ny - 2; ++p)
				A[p][p] = a_pp;
			//构建压力方程的右端项
			vector<complex<double>> RHS(Ny, 0.0);
			//计算前Ny-1个方程右端项，为此需要先计算0+1/3时间步上v对y的导数在谱空间的投影
			vector<complex<double>> dvdy_(Ny, 0.0);
			//利用一阶Chebyshev矩阵计算v对y的导数的谱分量
			for (int p = 0; p < Ny; ++p)
			{
				for (size_t J = p + 1; J < Ny; J += 2)
				{
					size_t wav = data.mesh.idx2(m, J, n);
					dvdy_[p] += d1[p][J] * s_1.v_[wav];
				}
			}
			for (int p = 0; p < Ny - 2; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				complex<double> r_p = (Im * Alpha * md * s_1.u_[wav] + dvdy_[p] + Im * Beta * nd * s_1.w_[wav]) / dt;
				RHS[p] = r_p;
			}
			//边界条件方程的右端项由边界条件确定
			complex<double> G_up(0.0), G_low(0.0);//上/下边界对应的右端项
			double k = 1.0;
			for (int p = 0; p < Ny; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				G_up += Im * Alpha * md * s_2.o_x_[wav] - Im * Beta * nd * s_2.o_z_[wav];
				G_low += (Im * Alpha * md * s_2.o_x_[wav] - Im * Beta * nd * s_2.o_z_[wav]) * k;
				k = -k;
			}
			G_up = G_up / Re;
			G_low = G_low / Re;
			//计算右端项
			RHS[Ny - 2] = G_up;
			RHS[Ny - 1] = G_low;
			//求解方程组
			size_t mn = (m + Nx / 2) * Ny * Nz + (n + Nz / 2) * Ny;
			Helmholtz(A, RHS, &s_1.p_[mn], data.mesh);
			//推进到下一子时间步
			for (int p = 0; p < Ny; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				s_1.u_[wav] += -Im * Alpha * md * s_1.p_[wav] * dt;
				s_1.w_[wav] += -Im * Beta * nd * s_1.p_[wav] * dt;
				//对速度分量v推进需要计算П^(1)
				complex<double> P1 = 0.0;
				for (size_t J = p + 1; J < Ny; J += 2) P1 += d1[p][J] * s_1.p_[data.mesh.idx2(m, J, n)];
				s_1.v_[wav] += -P1 * dt;
			}
		}
	}
	//**************************粘性步************************
	Mat Au(data.mesh.D2);
	BuildCoefficientMatrixOfVelocity(Au, data.mesh);
	//粘性步需要解三个速度方程组，其形式与压力方程非常相似，三个方向速度分量方程的系数矩阵完全一致
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int n = -Nz / 2; n < Nz / 2; ++n)
		{
			//对于每一对m、n设定系数矩阵对角项
			double a_ii = -pow(Alpha * m, 2) - pow(Beta * n, 2) - gamma * Re / dt;
			for (size_t i = 0; i < Ny - 2; ++i)
				Au[i][i] = a_ii;
			Mat Av(Au), Aw(Au);
			//构建压力方程的右端项，所需的亚格子应力在谱空间的投影已经计算过
			vector<complex<double>> ru(Ny, 0.0), rv(Ny, 0.0), rw(Ny, 0.0);
			for (int p = 0; p < Ny - 2; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				double md = static_cast<double>(m), nd = static_cast<double>(n);
				complex<double> sgs_u = (1 - q) * (Im * Alpha * md * s_2.tau_[0][0][wav] + s_2.dtau_xy_[wav] + Im * Beta * nd * s_2.tau_[0][2][wav]) / gamma,
					sgs_v = (1 - q) * (Im * Alpha * md * s_2.tau_[0][1][wav] + s_2.dtau_yy_[wav] + Im * Beta * nd * s_2.tau_[1][2][wav]) / gamma,
					sgs_w = (1 - q) * (Im * Alpha * md * s_2.tau_[0][2][wav] + s_2.dtau_yz_[wav] + Im * Beta * nd * s_2.tau_[2][2][wav]) / gamma;
				ru[p] = -Re * s_1.u_[wav] / dt + sgs_u;
				rv[p] = -Re * s_1.v_[wav] / dt + sgs_v;
				rw[p] = -Re * s_1.w_[wav] / dt + sgs_w;
			}
			//推进到下一子时间步
			size_t mn = (m + Nx / 2) * Ny * Nz + (n + Nz / 2) * Ny;
			Helmholtz(Au, ru, &s_1.u_[mn], data.mesh);
			Helmholtz(Av, rv, &s_1.v_[mn], data.mesh);
			Helmholtz(Aw, rw, &s_1.w_[mn], data.mesh);
		}
	}
}

void SecondStep(History& data, const unique_ptr<TransformEngine>& eng, TransformSpace& ts, const unique_ptr<SGSModel>& sgs)
{
	auto& s_2 = data.s_2(), & s_1 = data.s_1(), & s = data.s();
	int Nx = data.mesh.Nx, Ny = data.mesh.Ny, Nz = data.mesh.Nz, Nn = data.mesh.Nn;
	double Alpha = data.mesh.Alpha, Beta = data.mesh.Beta, dt = data.config.dt, Re = data.config.ReTau, q = data.config.q;
	double gamma = 11.0 / 6;
	//************************非线性步************************
	//先更新1时间步物理空间的数据
	s_1.SpectralToPhysical(eng, ts);
	//计算速度梯度
	s_1.ComputeVelocityGrad(eng);
	//先计算1时间步上的涡量，后面要用
	s_1.ComputeVorticity();
	//计算1时间步物理空间非线性项
	s_1.ComputeNonlinearTerms();
	//计算亚格子应力
	sgs->ComputeSGS(s_1.du, s_1.tau);
	//计算1时间步亚格子应力和非线性项在谱空间的投影
	s_1.ComputeSpectralNonlinearTermAndSGS(eng, ts);
	//时间推进到1+1/3时间步
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int n = -Nz / 2; n < Nz / 2; ++n)
		{
			for (int p = 0; p < Ny; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				double md = static_cast<double>(m), nd = static_cast<double>(n);
				complex<double> uRHS0 = s_2.Nu_[wav] + q * (Im * Alpha * md * s_2.tau_[0][0][wav] + s_2.dtau_xy_[wav] + Im * Beta * nd * s_2.tau_[0][2][wav]),
					uRHS1 = s_1.Nu_[wav] + q * (Im * Alpha * md * s_1.tau_[0][0][wav] + s_1.dtau_xy_[wav] + Im * Beta * nd * s_1.tau_[0][2][wav]);
				s.u_[wav] = s_1.u_[wav] + dt * (1.5 * uRHS1 - 0.5 * uRHS0);
				complex<double> vRHS0 = s_2.Nv_[wav] + q * (Im * Alpha * md * s_2.tau_[0][1][wav] + s_2.dtau_yy_[wav] + Im * Beta * nd * s_2.tau_[1][2][wav]),
					vRHS1 = s_1.Nv_[wav] + q * (Im * Alpha * md * s_1.tau_[0][1][wav] + s_1.dtau_yy_[wav] + Im * Beta * nd * s_1.tau_[1][2][wav]);
				s.v_[wav] = s_1.v_[wav] + dt * (1.5 * vRHS1 - 0.5 * vRHS0);
				complex<double> wRHS0 = s_2.Nw_[wav] + q * (Im * Alpha * md * s_2.tau_[0][2][wav] + s_2.dtau_yz_[wav] + Im * Beta * nd * s_2.tau_[2][2][wav]),
					wRHS1 = s_1.Nw_[wav] + q * (Im * Alpha * md * s_1.tau_[0][2][wav] + s_1.dtau_yz_[wav] + Im * Beta * nd * s_1.tau_[2][2][wav]);
				s.w_[wav] = s_1.w_[wav] + dt * (1.5 * wRHS1 - 0.5 * wRHS0);
			}
		}
	}
	//**************************压强步************************
	//构建压力方程系数矩阵，用二阶Chebyshev微分矩阵初始化
	Mat A(data.mesh.D2);
	//用于计算的一阶Chebyshev微分矩阵
	Mat d1(data.mesh.D1);
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int n = -Nz / 2; n < Nz / 2; ++n)
		{
			if (m == 0 && n == 0)
			{
				size_t mn = Nn / 2 + Nz * Ny / 2;
				for (int p = 0; p < Ny; ++p)s.p_[mn + p] = 0;
				continue;
			}
			double md = static_cast<double>(m), nd = static_cast<double>(n);
			//对于每一对m、n系数矩阵对角项为相应的波数
			double a_pp = -pow(Alpha * m, 2) - pow(Beta * n, 2);
			for (int p = 0; p < Ny - 2; ++p)
				A[p][p] = a_pp;
			//构建压力方程的右端项
			vector<complex<double>> RHS(Ny, 0.0);
			//计算前Ny-1个方程右端项，为此需要先计算s+1/3时间步上v对y的导数在谱空间的投影
			vector<complex<double>> dvdy_(Ny, 0.0);
			//利用一阶Chebyshev矩阵计算v对y的导数的谱分量
			for (int p = 0; p < Ny; ++p)
			{
				for (size_t J = p + 1; J < Ny; J += 2)
				{
					size_t wav = data.mesh.idx2(m, J, n);
					dvdy_[p] += d1[p][J] * s_1.v_[wav];
				}
			}
			for (int p = 0; p < Ny - 2; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				complex<double> r_p = (Im * Alpha * md * s.u_[wav] + dvdy_[p] + Im * Beta * nd * s.w_[wav]) / dt;
				RHS[p] = r_p;
			}
			//边界条件方程的右端项由边界条件确定
			//需要先计算涡量在谱空间的投影，因此还需要计算dv/dz和dv/dx
			vector<complex<double>> G_up(2, complex<double>(0.0)), G_low(2, complex<double>(0.0));//上/下边界对应的右端项
			double k = 1.0;
			for (int p = 0; p < Ny; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				G_up[0] += Im * Alpha * md * s_2.o_x_[wav] - Im * Beta * nd * s_2.o_z_[wav];
				G_low[0] += (Im * Alpha * md * s_2.o_x_[wav] - Im * Beta * nd * s_2.o_z_[wav]) * k;
				G_up[1] += Im * Alpha * md * s_1.o_x_[wav] - Im * Beta * nd * s_1.o_z_[wav];
				G_low[1] += (Im * Alpha * md * s_1.o_x_[wav] - Im * Beta * nd * s_1.o_z_[wav]) * k;
				k = -k;
			}
			for (int I = 0; I < 2; ++I)
			{
				G_up[I] = G_up[I] / Re;
				G_low[I] = G_low[I] / Re;
			}
			RHS[Ny - 2] = 1.5 * G_up[1] - 0.5 * G_up[0];
			RHS[Ny - 1] = 1.5 * G_low[1] - 0.5 * G_low[0];
			//求解方程组
			size_t mn = (m + Nx / 2) * Ny * Nz + (n + Nz / 2) * Ny;
			Helmholtz(A, RHS, &s.p_[mn], data.mesh);
			//推进到下一子时间步
			for (int p = 0; p < Ny; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				s.u_[wav] += -Im * Alpha * md * s.p_[wav] * dt;
				s.w_[wav] += -Im * Beta * nd * s.p_[wav] * dt;
				//对速度分量v推进需要计算П^(1)
				complex<double> P1 = 0.0;
				for (size_t J = p + 1; J < Ny; J += 2) P1 += d1[p][J] * s.p_[data.mesh.idx2(m, J, n)];
				s.v_[wav] += -P1 * dt;
			}
		}
	}
	//**************************粘性步************************
	Mat Au(data.mesh.D2);
	BuildCoefficientMatrixOfVelocity(Au, data.mesh);
	//粘性步需要解三个速度方程组，其形式与压力方程非常相似，三个方向速度分量方程的系数矩阵完全一致
	for (int m = -Nx / 2; m < Nx / 2; ++m)
	{
		for (int n = -Nz / 2; n < Nz / 2; ++n)
		{
			//对于每一对m、n设定系数矩阵对角项
			double a_ii = -pow(Alpha * m, 2) - pow(Beta * n, 2) - Re / (gamma * dt);
			for (size_t i = 0; i < Ny - 2; ++i)
				Au[i][i] = a_ii;
			Mat Av(Au), Aw(Au);
			//构建压力方程的右端项，所需的亚格子应力在谱空间的投影已经计算过
			vector<complex<double>> ru(Ny, 0.0), rv(Ny, 0.0), rw(Ny, 0.0);
			for (int p = 0; p < Ny - 2; ++p)
			{
				size_t wav = data.mesh.idx2(m, p, n);
				double md = static_cast<double>(m), nd = static_cast<double>(n);
				complex<double> sgs_u = (1 - q) * (Im * Alpha * md * s_1.tau_[0][0][wav] + s_1.dtau_xy_[wav] + Im * Beta * nd * s_1.tau_[0][2][wav]) / gamma,
					sgs_v = (1 - q) * (Im * Alpha * md * s_1.tau_[0][1][wav] + s_1.dtau_yy_[wav] + Im * Beta * nd * s_1.tau_[1][2][wav]) / gamma,
					sgs_w = (1 - q) * (Im * Alpha * md * s_1.tau_[0][2][wav] + s_1.dtau_yz_[wav] + Im * Beta * nd * s_1.tau_[2][2][wav]) / gamma;
				ru[p] = -Re * s.u_[wav] / dt + sgs_u;
				rv[p] = -Re * s.v_[wav] / dt + sgs_v;
				rw[p] = -Re * s.w_[wav] / dt + sgs_w;
			}
			//推进到下一子时间步
			size_t mn = (m + Nx / 2) * Ny * Nz + (n + Nz / 2) * Ny;
			Helmholtz(Au, ru, &s.u_[mn], data.mesh);
			Helmholtz(Av, rv, &s.v_[mn], data.mesh);
			Helmholtz(Aw, rw, &s.w_[mn], data.mesh);
		}
	}
}