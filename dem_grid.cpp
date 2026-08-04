#include "dem_grid.h"

Grid::Grid(int nx, int ny, int nz, double cellSize, const CFDData& domain):
	nx_(nx), ny_(ny), nz_(nz), cellSize_(cellSize),
	overshootX_(0.5 * (nx * cellSize - domain.lx)),
	overshootY_(0.5 * (ny * cellSize - domain.h)),
	overshootZ_(0.5 * (nz * cellSize - domain.lz)){ }


int Grid::getCellIndex(double x, double y, double z) const
{
	int ix = static_cast<int>((x + overshootX_) / cellSize_);
	int iy = static_cast<int>((y + overshootY_) / cellSize_);
	int iz = static_cast<int>((z + overshootZ_) / cellSize_);
	return iz * nx_ * ny_ + iy * nx_ + ix;
}

int Grid::nCellsX() const { return nx_; }
int Grid::nCellsY() const { return ny_; }
int Grid::nCellsZ() const { return nz_; }

const std::vector<Offset>& Grid::neighborOffsets() const { return neighbors_; }