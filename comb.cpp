// lotka-volterra dynamics with eigen array class, RK4 instead of euler
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <fstream>
#include <omp.h>
#include <Eigen/Dense>
#include <iomanip>
using namespace Eigen;
using namespace std;

random_device generator;
mt19937 twist(generator());


int main()
{

	ofstream heat;
	heat.open("allcolour3.txt");
	
	ofstream five;
	five.open("allfive3.txt");
	
	for (int mi = 0; mi < 12; mi ++)
	{
		for (int gi = 0; gi < 6; gi ++)
		{
			for (int si = 0; si < 21; si ++)
			{
				string measure;
				ifstream mfile("five3_" + to_string(mi) + "_" + to_string(gi) + "_" + to_string(si) + ".txt");
				string colour;
				ifstream cfile("colour3_" + to_string(mi) + "_" + to_string(gi) + "_" + to_string(si) + ".txt");
				
				if (cfile.is_open())
				{
					getline (cfile, colour);
					cfile.close();
				}
				else {cout << "badness" << endl;}
				if (mfile.is_open())
				{
					getline (cfile, measure);
					mfile.close();
				}
				else {cout << "badness" << endl;}
				
				heat << colour;
				five << measure;
				
				if (si != 20) {heat << ","; five << ",";}
			}
			if (gi != 5) {heat << ","; five << ",";}
		}
		if (mi != 11) {heat << ","; five << ",";}
		heat.close(); five.close();
	}
	return 0;
}
	
				
			
			


















