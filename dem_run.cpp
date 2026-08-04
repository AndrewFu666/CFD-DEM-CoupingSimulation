#include "dem_run.h"

void DEMRun::run(double dt_DEM, const CFDData& cfddata)
{
	auto candidates = findPotentialCollisions(grid, particles, cfddata);
	detectAndResolveOverlaps(particles, cfddata, candidates);
	integrate(cfddata, particles, dt_DEM);
}