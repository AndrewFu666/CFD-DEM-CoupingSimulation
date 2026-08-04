#pragma once
#include "dem_cfd_data.h"

struct Offset { int dx, dy, dz; };

class Grid
{
public:
	Grid(int nx, int ny, int nz, double cellSize, const CFDData& domain);
    int nCellsX() const, nCellsY() const, nCellsZ() const;
    const std::vector<Offset>& neighborOffsets() const;
	int getCellIndex(double x, double y, double z) const;
private:
    int nx_, ny_, nz_;//三个方向DEM网格数量
    double cellSize_;//网格边长
    double overshootX_, overshootY_, overshootZ_;//边界网格超出物理边界的量
    const std::vector<Offset> neighbors_ =
    {
        {1,-1,-1}, {1,-1,0}, {1,-1,1},
        {1, 0,-1}, {1, 0,0}, {1, 0,1},
        {1, 1,-1}, {1, 1,0}, {1, 1,1},
        {0, 1,-1}, {0, 1,0}, {0, 1,1},
        {0, 0, 1}
    };//邻居网格索引偏移量
};