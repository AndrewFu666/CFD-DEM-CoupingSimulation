#pragma once
#include <vector>

struct CFDData 
{
    double lx, h, lz;//根据长度标尺和无量纲几何参数重建的真实计算域尺寸
    double density, viscosity;//流体的密度和动力粘性系数
    int nx, ny, nz;//x，y和z方向的网格数量
    std::vector<double> y;//真实的网格节点y坐标
    std::vector<double> u, v, w;//真实的速度分量
};