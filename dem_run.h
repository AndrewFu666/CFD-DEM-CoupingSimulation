#pragma once
#include "dem_particles.h"
#include "dem_cfd_data.h"
#include "dem_grid.h"
#include "dem_collision_predict.h"
#include "dem_collision_solve.h"
#include "dem_time_integrator.h"


class DEMRun
{
public:
	DEMRun(Particles& p, const Grid& g) : particles(p), grid(g) { }
	void run(double dt_DEM, const CFDData& cfddata);
private:
	Particles& particles;
	const Grid& grid;
};