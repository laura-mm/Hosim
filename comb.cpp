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
	heat.open("l.txt");	

	for (int g = 0; g <= 100; g++)
	{
		for (int s = 0; s <= 100; s++)
		{
			string colour;
			string filename = "l_" + to_string(g) + "_" + to_string(s) + ".txt";
			ifstream file(filename);
			if (file.is_open())
			{
				getline (file, colour);
				file.close();
			}
			else {cout << filename << endl; colour = "1,1,1";}
			heat << colour;
			if (s != 100) heat << ",";
		}
		if (g != 100) heat << ",";
	}
	heat.close();

	return 0;
}
	
				
			
			


















