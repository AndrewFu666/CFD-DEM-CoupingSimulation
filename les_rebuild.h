#pragma once
#include "les_field.h"
#include "les_mesh.h"
#include "les_types.h"
#include "dem_cfd_data.h"

struct Attribute
{
	double density;//流体密度
	double lenScale;//半槽高度
	double vScale;//摩擦速度
	double viscosity;//流体动力粘度
};

CFDData rebuildFlowFieldForDEM(const Field& field, const Mesh& mesh, const Attribute& attribute);

void updateCFDData(const Var& u, const Var& v, const Var& w, const Attribute& attribute, CFDData& cfddata);