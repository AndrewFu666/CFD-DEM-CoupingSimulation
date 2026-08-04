#pragma once

class Config
{
public:
	Config();
	//物性参数
	const double ReTau;//摩擦雷诺数
	const double Visc;//运动粘度
	const double Fx;//流向驱动压力梯度
	//时间控制
	const double dt;//时间推进步长
	const size_t MaxStp;//最大时间步数
	const size_t StatsStp;//统计起始步数
	//初始化参数
	const double Um;//初始流向平均速度
	const double A;//初始随机扰动幅值
	//决定分配到非线性步的亚格子应力项占比的系数
	const double q;
};