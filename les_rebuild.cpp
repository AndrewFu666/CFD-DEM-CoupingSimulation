#include "les_rebuild.h"
#include "les_types.h"

using namespace std;

CFDData rebuildFlowFieldForDEM(const Field& field, const Mesh& mesh, const Attribute& attribute)
{
	CFDData cfddata;
	cfddata.density = attribute.density;
	cfddata.nx = mesh.Nx;
	cfddata.ny = mesh.Ny;
	cfddata.nz = mesh.Nz;
	cfddata.lx = 4.0 * PI * 2 * attribute.lenScale;
	cfddata.h = 2 * attribute.lenScale;
	cfddata.lz = 2.0 * PI * 2 * attribute.lenScale;
	cfddata.viscosity = attribute.viscosity;
	cfddata.y.resize(mesh.Ny);
	cfddata.u.resize(mesh.Nn);
	cfddata.v.resize(mesh.Nn);
	cfddata.w.resize(mesh.Nn);
	for (int i = 0; i < mesh.Ny; ++i)
		cfddata.y[i] = mesh.y[i] * attribute.lenScale;
	for (int i = 0; i < mesh.Nn; ++i)
	{
		cfddata.u[i] = field.u[i] * attribute.vScale;
		cfddata.v[i] = field.v[i] * attribute.vScale;
		cfddata.w[i] = field.w[i] * attribute.vScale;
	}
	return cfddata;
}

void updateCFDData(const Var& u, const Var& v, const Var& w, const Attribute& attribute, CFDData& cfddata)
{
	for (int i = 0; i < u.size(); ++i)
	{
		cfddata.u[i] = u[i] * attribute.vScale;
		cfddata.v[i] = v[i] * attribute.vScale;
		cfddata.w[i] = w[i] * attribute.vScale;
	}
}