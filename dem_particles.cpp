#include "dem_particles.h"
#include <fstream>
#include <string>
#include <cassert>

using namespace std;

void Particles::output(double solutionTime)
{
	static bool firstCall = true;
	string name = "dem_result.dat";//结果会写在dem_result中
	ofstream file(name, firstCall ? ios::trunc : ios::app);
	assert(file.is_open());
	if (firstCall)
	{
		file << "TITLE = \"Particles\"" << endl;
		file << "VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"D\"" << endl;
		firstCall = false;
	}
	file << "ZONE I=" << count << ", F=POINT, SOLUTIONTIME=" << solutionTime <<", STRANDID=1" << endl;
	for (int i = 0; i < count; ++i)
	{
		file << x[i] << " " << y[i] << " " << z[i] << " "
			<< vx[i] << " " << vy[i] << " " << vz[i] << " " << diameter << endl;
	}
}