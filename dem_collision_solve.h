#pragma once
#include "dem_particles.h"
#include "dem_cfd_data.h"
#include "dem_collision_list.h"

void detectAndResolveOverlaps(Particles& p, const CFDData& f, const CollisionList& candidates);
