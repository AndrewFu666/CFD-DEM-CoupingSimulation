#pragma once
#include "dem_cfd_data.h"
#include "dem_particles.h"

struct InitBox { double xMin, xMax, yMin, yMax, zMin, zMax; };

void initializeParticles(const CFDData& cfddata, Particles& p,const InitBox& box, unsigned seed = 0);
