#include "dem_collision_solve.h"
#include "dem_vector_operate.h"

void detectAndResolveOverlaps(Particles& p, const CFDData& f, const CollisionList& candidates)
{
	//边界碰撞监测
	for (const auto& wallContact : candidates.wallContacts)
	{
		int IDp = wallContact.first;//球的ID
		auto IDb = wallContact.second;//边界的ID
		double radius = 0.5 * p.diameter;
		double safety = 0.005 * radius;//为了避免舍入误差而引入的余量
		switch (IDb)
		{
			case Boundaries::yMax:
				if (f.h - p.y[IDp] - radius < 0)
				{
					p.vy[IDp] = -p.vy[IDp];
					p.y[IDp] = f.h - radius - safety;
				}
				break;
			case Boundaries::yMin:
				if (p.y[IDp] - radius < 0)
				{
					p.vy[IDp] = -p.vy[IDp];
					p.y[IDp] = radius + safety;
				}
				break;
			case Boundaries::zMin:
				if (p.z[IDp] - radius < 0)
				{
					p.vz[IDp] = -p.vz[IDp];
					p.z[IDp] = radius + safety;
				}
				break;
			case Boundaries::zMax:
				if (f.lz - p.z[IDp] - radius < 0)
				{
					p.vz[IDp] = -p.vz[IDp];
					p.z[IDp] = f.lz - radius - safety;
				}
				break;
			case Boundaries::xMax:
				if (f.lx - p.x[IDp] - radius < 0)
				{
					p.vx[IDp] = -p.vx[IDp];
					p.x[IDp] = f.lx - radius - safety;
				}
				break;
			case Boundaries::xMin:
				if (p.x[IDp] - radius < 0)
				{
					p.vx[IDp] = -p.vx[IDp];
					p.x[IDp] = radius + safety;
				}
				break;
		}
	}
	//颗粒间碰撞监测
	for (const auto& particlePair : candidates.particlePairs)
	{
		//第一个球的参数
		int ID1 = particlePair.first;
		double radius1 = 0.5 * p.diameter;
		double mass1 = p.mass;
		Vec3 v1{ p.vx[ID1], p.vy[ID1], p.vz[ID1] };
		Vec3 w1{ p.omegaX[ID1], p.omegaY[ID1], p.omegaZ[ID1] };
		Vec3 r1{ p.x[ID1], p.y[ID1], p.z[ID1] };
		//第二个球的参数
		int ID2 = particlePair.second;
		double radius2 = 0.5 * p.diameter;
		double mass2 = p.mass;
		Vec3 v2{ p.vx[ID2], p.vy[ID2], p.vz[ID2] };
		Vec3 w2{ p.omegaX[ID2], p.omegaY[ID2], p.omegaZ[ID2] };
		Vec3 r2{ p.x[ID2], p.y[ID2], p.z[ID2] };
		//颗粒间距离检测
		double overlap = len(r1 - r2) - radius1 - radius2;
		//若两球接触则计算碰撞
		if (overlap < 0)
		{
			Vec3 n = norm(r1 - r2);//两颗粒发生碰撞时的法向单位矢量
			Vec3 g = v1 - v2;//颗粒质心碰前相对速度
			Vec3 g_ct = g + cross(w1, n) * radius1 + cross(w2, n) * radius2;//碰前碰撞点的相对速度
			Vec3 t = norm(g_ct);//切向单位矢量
			double ng = dot(n, g), gct = len(g_ct);
			//颗粒间发生切向滑动（True），颗粒间不发生切向滑动（False）
			bool slipCriterion = (ng / gct) < (2 / (7 * (1 + p.e)));
			if (slipCriterion)
			{
				//计算碰撞后的速度
				double k1v = (1 + p.e) * mass2 * ng / (mass1 + mass2), k2v = (1 + p.e) * mass1 * ng / (mass1 + mass2);//标量因子
				Vec3 dv0 = n - t * p.mu;//矢量因子
				Vec3 dv1 = dv0 * k1v, dv2 = dv0 * k2v;//速度变化量
				p.vx[ID1] = (v1 - dv1).x; p.vy[ID1] = (v1 - dv1).y; p.vz[ID1] = (v1 - dv1).z;
				p.vx[ID2] = (v2 + dv2).x; p.vy[ID2] = (v2 + dv2).y; p.vz[ID2] = (v2 + dv2).z;
				//计算碰撞后角速度
				double k1w = (1 + p.e) * mass2 * ng * 2.5 / (radius1 * (mass1 + mass2)),
					k2w = (1 + p.e) * mass1 * ng * 2.5 / (radius2 * (mass1 + mass2));//标量因子
				Vec3 dw0 = cross(n, t);//矢量因子
				Vec3 dw1 = dw0 * k1w, dw2 = dw0 * k2w;//角速度变化量
				p.omegaX[ID1] = (w1 - dw1).x; p.omegaY[ID1] = (w1 - dw1).y; p.omegaZ[ID1] = (w1 - dw1).z;
				p.omegaX[ID2] = (w2 - dw2).x; p.omegaY[ID2] = (w2 - dw2).y; p.omegaZ[ID2] = (w2 - dw2).z;
			}
			else
			{
				//计算碰撞后速度
				double k1n = (1 + p.e) * mass2 * ng / (mass1 + mass2), k1t = 2.0 * gct * mass2 / (7 * (mass1 + mass2)),
					k2n = (1 + p.e) * mass1 * ng / (mass1 + mass2), k2t = 2.0 * gct * mass1 / (7 * (mass1 + mass2));//标量因子
				Vec3 dv1 = n * k1n - t * k1t, dv2 = n * k2n - t * k2t;
				p.vx[ID1] = (v1 - dv1).x; p.vy[ID1] = (v1 - dv1).y; p.vz[ID1] = (v1 - dv1).z;
				p.vx[ID2] = (v2 + dv2).x; p.vy[ID2] = (v2 + dv2).y; p.vz[ID2] = (v2 + dv2).z;
				//计算碰撞后角速度
				double k1 = mass2 * gct * 5.0 / (7 * radius1 * (mass1 + mass2)),
					k2 = mass1 * gct * 5.0 / (7 * radius2 * (mass1 + mass2));//标量因子
				Vec3 dw0 = cross(n, t);//矢量因子
				Vec3 dw1 = dw0 * k1, dw2 = dw0 * k2;//角速度变化量
				p.omegaX[ID1] = (w1 - dw1).x; p.omegaY[ID1] = (w1 - dw1).y; p.omegaZ[ID1] = (w1 - dw1).z;
				p.omegaX[ID2] = (w2 - dw2).x; p.omegaY[ID2] = (w2 - dw2).y; p.omegaZ[ID2] = (w2 - dw2).z;
			}
			//将两球分开微小距离
			Vec3 dr = n * (0.005 * std::max(radius1, radius2) - overlap);
			p.x[ID1] += dr.x; p.y[ID1] += dr.y; p.z[ID1] += dr.z;
			p.x[ID2] -= dr.x; p.y[ID2] -= dr.y; p.z[ID2] -= dr.z;
		}
	}
}