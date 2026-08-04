#pragma once
#include "les_field.h"
#include "les_SGSmodel.h"

void NonlinearStep(History& data, const std::unique_ptr<TransformEngine>& eng,
	TransformSpace& ts, const std::unique_ptr<SGSModel>& sgs);

void Helmholtz(Mat A, std::vector<std::complex<double>>& b, std::complex<double>* x, const Mesh& m);

void PressureStep(History& data);

void ViscosityStep(History& data);

void FirstStep(History& data, const std::unique_ptr<TransformEngine>& eng,
	TransformSpace& ts, const std::unique_ptr<SGSModel>& sgs);

void SecondStep(History& data, const std::unique_ptr<TransformEngine>& eng,
	TransformSpace& ts, const std::unique_ptr<SGSModel>& sgs);