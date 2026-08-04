#pragma once

//定义代表各个边界的值
enum class Boundaries { xMin, xMax, yMin, yMax, zMin, zMax };

struct CollisionList 
{
    std::vector<std::pair<int, int>> particlePairs;
    std::vector<std::pair<int, Boundaries>> wallContacts;
};