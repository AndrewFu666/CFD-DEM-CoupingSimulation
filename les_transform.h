#pragma once
#include <vector>
#include <complex>
#include "les_mesh.h"

class TransformSpace
{
public:
	TransformSpace(const Mesh& m) :yz(27 * m.Nn, 0.0), y(27 * m.Nn, 0.0), spec(27 * m.NnEx, 0.0), Nvar(0), mesh(m) {}
	int Nvar;//传入变换空间的物理量个数
	std::vector<std::complex<double>> yz;//所有物理量扁平化存储
	std::vector<std::complex<double>> y;
	std::vector<std::complex<double>> spec;
	void reset() { Nvar = 0; }
	void AddPhysicalField(const Var& var);
	void AddSpectralField(const Var_& var_);
	void OutPhysicalField(Var& var);
	void OutSpectralField(Var_& var_);
private:
	const Mesh& mesh;
};

class TransformEngine//物理空间-谱空间变换引擎接口
{
public:
	virtual void FFT(std::complex<double>* begin, size_t N, bool invert, bool center) = 0;
	virtual void Forward(TransformSpace& ts) = 0;
	virtual void Inverse(TransformSpace& ts) = 0;
	virtual ~TransformEngine() = default;
};

class MyEngine :public TransformEngine
{
public:
	MyEngine(const Mesh& m) : mesh(m) {}
	void FFT(std::complex<double>* begin, size_t N, bool invert, bool center)override;
	void Forward(TransformSpace& ts)override;
	void Inverse(TransformSpace& ts)override;
private:
	const Mesh& mesh;
};
