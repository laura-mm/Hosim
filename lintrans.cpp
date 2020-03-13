// transforming colours for the linear condor
// black to white
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <fstream>
#include <omp.h>
#include <iomanip>
using namespace std;




int main()
{
	ofstream heat;
	heat.open("ho2-_lt.txt");
	ifstream file("ho2-.txt");
	
	for (int i = 0; i <= 10200; i ++)
	{
		string r;
		getline(file, r, ',');
		double R = stod(r);
		
		string g;
		getline(file, g, ',');
		double G = stod(g);

		string b;
		getline(file, b, ',');
		double B = stod(b);

		if (R + G + B == 3.0) {R = 0.0; G = 0.0; B = 0.0;}

		else
		{
			double div = 1.0 - (R + G + B);

			R += div;
			G += div;
			B += div;
		}
		

		heat << R << "," << G << "," << B;
		if (i != 10200) heat << ",";
	}
	heat.close();

	return 0;
}
		
	
				
			
			


















