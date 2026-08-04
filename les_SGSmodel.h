#pragma once
#include "les_mesh.h"
#include "les_config.h"
#include <vector>
#include <complex>

class SGSModel//亚格子应力模型接口
{
public:
	virtual void ComputeSGS(const std::vector<std::vector<Var>>& du, std::vector<std::vector<Var>>& tau) const = 0;
	virtual ~SGSModel() = default;
};

class SmagorinskyModel :public SGSModel
{
public:
	SmagorinskyModel(const Mesh& m, const Config& con);
	void ComputeSGS(const std::vector<std::vector<Var>>& du, std::vector<std::vector<Var>>& tau) const override;
private:
	const Mesh& mesh;
	const Config& config;
	std::vector<double> Cs;//亚格子涡粘系数
	std::vector<double> delta;//当量网格尺寸的平方
};
