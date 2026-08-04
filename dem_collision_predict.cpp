#include "dem_collision_predict.h"
#include <unordered_map>

using namespace std;

CollisionList findPotentialCollisions(const Grid& grid, const Particles& p, const CFDData& cfd)
{
	//统计颗粒所在网格ID算法
	vector<int> cellData(p.count, 0);
	unordered_map<int, int> cellCount, cellStart;
	cellCount.reserve(p.count);//预估容量，避免rehash
	//第一次遍历所有颗粒统计各个非空网格中颗粒的数量
	vector<int> particleCellID(p.count, -1);//存储颗粒所在网格的ID，避免重复计算
	for (unsigned particleID = 0; particleID < p.count; ++particleID)
	{
		particleCellID[particleID] = grid.getCellIndex(p.x[particleID], p.y[particleID], p.z[particleID]);
		cellCount[particleCellID[particleID]]++;
	}
	int pos = 0;
	for (const auto& ID_count : cellCount)
	{
		cellStart[ID_count.first] = pos;
		pos += ID_count.second;
	}
	//cellStart中存储各个非空网格中的颗粒ID集在cellData中存储的位置
	//writePos是它的副本，用来在第二次遍历中向cellData中写入
	unordered_map<int, int> writePos(cellStart);
	for (unsigned particleID = 0; particleID < p.count; particleID++)
	{
		cellData[writePos[particleCellID[particleID]]] = particleID;
		writePos[particleCellID[particleID]]++;
	}
	//邻居搜索算法
	CollisionList candidates;
	vector<int> neighborP;//用来存储每一个网格的邻居网格中的颗粒
	for (const auto& ID_count : cellCount)//遍历每一个有颗粒的网格
	{
		int cellID = ID_count.first, count = ID_count.second, start = cellStart[cellID];
		//先根据当前网格ID反向计算其空间索引
		int iz = cellID / (grid.nCellsX() * grid.nCellsY()), 
			rem = cellID % (grid.nCellsX() * grid.nCellsY()), 
			iy = rem / grid.nCellsX(), ix = rem % grid.nCellsX();
		//搜索邻居网格（一半）并记录其中所有颗粒
		neighborP.clear();
		for (const auto& neighbor : grid.neighborOffsets())
		{
			int ix_n = ix + neighbor.dx, iy_n = iy + neighbor.dy, iz_n = iz + neighbor.dz;
			if (ix_n < 0 || ix_n >= grid.nCellsX() || iy_n < 0 || iy_n >= grid.nCellsY() || iz_n < 0 || iz_n >= grid.nCellsZ())
				continue;
			int ID_n = iz_n * grid.nCellsX() * grid.nCellsY() + iy_n * grid.nCellsX() + ix_n;
			auto it = cellStart.find(ID_n);
			if (it != cellStart.end())
			{
				for (int j = 0; j < cellCount[ID_n]; ++j)
					neighborP.push_back(cellData[it->second + j]);
			}
		}
		for (int i = 0; i < count; ++i)//遍历当前网格中的每一个颗粒
		{
			int p1 = cellData[start + i];
			//先进行壁面检测（规定代表壁面的对象写在后）
			if (p.x[p1] < 0.5 * p.diameter)candidates.wallContacts.push_back({ p1,Boundaries::xMin });
			else if (cfd.lx - p.x[p1] < 0.5 * p.diameter)candidates.wallContacts.push_back({ p1,Boundaries::xMax });
			if (p.y[p1] < 0.5 * p.diameter)candidates.wallContacts.push_back({ p1,Boundaries::yMin });
			else if (cfd.h - p.y[p1] < 0.5 * p.diameter)candidates.wallContacts.push_back({ p1,Boundaries::yMax });
			if (p.z[p1] < 0.5 * p.diameter)candidates.wallContacts.push_back({ p1,Boundaries::zMin });
			else if (cfd.lz - p.z[p1] < 0.5 * p.diameter)candidates.wallContacts.push_back({ p1,Boundaries::zMax });
			//进行网格内配对
			for (int j = i + 1; j < count; ++j)
			{
				int p2 = cellData[start + j];
				candidates.particlePairs.push_back({ p1,p2 });
			}
			//进行邻居配对
			for (const auto& p2 : neighborP)candidates.particlePairs.push_back({ p1,p2 });
		}
	}
	return candidates;
}