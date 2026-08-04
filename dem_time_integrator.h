#pragma once
#include "dem_cfd_data.h"
#include "dem_particles.h"

void integrate(const CFDData& f, Particles& p, double dt_DEM);
