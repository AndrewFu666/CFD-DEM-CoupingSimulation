#pragma once
#include "les_types.h"
#include "les_config.h"
#include "les_mesh.h"
#include "les_transform.h"
#include <vector>
#include <complex>
#include <deque>

class Field
{
public:
	Field(const Mesh& m, const Config& con);
	//物理空间流场初始化函数
	void Initialize();
	//快速傅里叶（逆）变换实现物理空间和谱空间之间的转换
	void PhysicalToSpectral(const std::unique_ptr<TransformEngine>& eng, TransformSpace& ts);//物理空间到谱空间
	void SpectralToPhysical(const std::unique_ptr<TransformEngine>& eng, TransformSpace& ts);//谱空间到物理空间
	void ComputeVelocityGrad(const std::unique_ptr<TransformEngine>& eng);//根据速度计算速度梯度张量
	void ComputeVorticity();//计算涡量
	void ComputeNonlinearTerms();//计算非线性项
	void ComputeSpectralNonlinearTermAndSGS(const std::unique_ptr<TransformEngine>& eng, TransformSpace& ts);
	void Output(const std::unique_ptr<TransformEngine>& eng, TransformSpace& ts);//输出计算结果
	Var u, v, w, p;//物理空间速度和压强
	Var_ u_, v_, w_, p_;//速度和压强的谱空间分量
	std::vector<std::vector<Var>> du;//速度梯度张量
	Var o_x;//涡量的x分量
	Var o_y;//涡量的y分量
	Var o_z;//涡量的z分量
	Var_ o_x_;//涡量的x分量的谱空间投影
	Var_ o_y_;//涡量的y分量的谱空间投影
	Var_ o_z_;//涡量的z分量的谱空间投影
	Var Nu;//u方程的非线性项
	Var Nv;//v方程的非线性项
	Var Nw;//w方程的非线性项
	Var_ Nu_;//Nu的谱空间投影
	Var_ Nv_;//Nv的谱空间投影
	Var_ Nw_;//Nw的谱空间投影
	std::vector<std::vector<Var>> tau;//亚格子应力张量
	std::vector<std::vector<Var_>> tau_;//亚格子应力的谱空间投影
	Var_ dtau_xy_;//需要对y求导的亚格子应力分量对y的导数的谱空间投影
	Var_ dtau_yy_;
	Var_ dtau_yz_;
private:
	const Mesh& mesh;
	const Config& config;
};

class History//用来存储显式时间推进所需要的当前时间步、上一时间步以及上上时间步的流场数据
{
public:
	History(const Mesh& m, const Config& con) :history(3, Field(m, con)), next(m, con), mesh(m), config(con) {}
	//三个函数分别调出s、s-1和s-2时间步的流场数据用于时间分裂推进计算
	Field& s() { return history[2]; }
	Field& s_1() { return history[1]; }
	Field& s_2() { return history[0]; }
	Field next;
	void update();
	const Mesh& mesh;
	const Config& config;
private:
	std::deque<Field> history;
};