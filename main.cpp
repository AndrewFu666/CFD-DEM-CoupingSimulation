#include "les_config.h"
#include "les_time_splitting_solver.h"
#include "les_rebuild.h"
#include "dem_run.h"
#include "dem_initialization.h"
#include <iostream>

using namespace std;

int main()
{
	//初始化流场属性
	const Config config;
	const Mesh mesh(64, 65, 64, 2, 4.0 * PI * 2, 2.0 * PI * 2);
	unique_ptr<SGSModel> sgsmodel = make_unique<SmagorinskyModel>(mesh, config);
	unique_ptr<TransformEngine> trans_eng = make_unique<MyEngine>(mesh);
	History data(mesh, config);
	TransformSpace trans_sp(mesh);
	data.s_2().Initialize();
	//初始化颗粒集群
	Attribute attribute{ 1.225,0.02,0.135,1.5e-5 };//空气物理属性
	CFDData cfddata = rebuildFlowFieldForDEM(data.s_2(), mesh, attribute);
	Grid DEMGrid(1008, 82, 504, 5e-4, cfddata);
	Particles particles(1e-4, 1.3e-9, 10000, 1.0, 0.8);//颗粒直径、质量、数量、碰撞恢复系数、摩擦系数
	InitBox box{ 0.25 * cfddata.lx,0.75 * cfddata.lx,0.25 * cfddata.h,
		0.75 * cfddata.h,0.25 * cfddata.lz,0.75 * cfddata.lz };//生成颗粒的区域为流场中间的一块区域
	initializeParticles(cfddata, particles, box, 123);//固定种子123
	DEMRun driver(particles, DEMGrid);//DEM计算驱动器
	double solutionTime = 0;//dem求解时间
	double dt_DEM = 0.05 * config.dt * attribute.lenScale / attribute.vScale;
	//开始迭代
	FirstStep(data, trans_eng, trans_sp, sgsmodel);
	cout << 1 << endl;
	SecondStep(data, trans_eng, trans_sp, sgsmodel);
	cout << 2 << endl;
	for (int stp = 3; stp <= config.MaxStp; ++stp)
	{
		NonlinearStep(data, trans_eng, trans_sp, sgsmodel);
		PressureStep(data);
		ViscosityStep(data);
		if (stp > 5000)//流场迭代5000步后开始与DEM耦合计算
		{
			updateCFDData(data.s().u, data.s().v, data.s().w, attribute, cfddata);
			for (int i = 1; i < 21; ++i)
			{
				driver.run(dt_DEM, cfddata);
				if (i % 5 == 0)//每5次DEM迭代计算输出一次结果
				{
					particles.output(solutionTime);
					solutionTime += 5 * dt_DEM;
				}
			}
		}
		cout << stp << endl;
	}
	data.s().Output(trans_eng, trans_sp);
}

