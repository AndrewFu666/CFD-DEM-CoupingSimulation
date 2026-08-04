#pragma once
#include "dem_grid.h"
#include "dem_particles.h"
#include "dem_collision_list.h"
#include "dem_cfd_data.h"

CollisionList findPotentialCollisions(const Grid& grid, const Particles& p, const CFDData& cfd);