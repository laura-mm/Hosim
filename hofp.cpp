// changed z2 to match the paper
// also simplified parameters
// z1 <= z2
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <iomanip>
#include <Eigen/Dense>
//#include <eigen3/Eigen/Dense>
using namespace Eigen;
using namespace std;

IOFormat pp(20, DontAlignCols, ",", ",", "", "", "", "");
IOFormat cc(6, DontAlignCols, ",", ",", "", "", "", "");

string carryon;
bool info = false;
double p;
double gama;
double mu;
double k;

// first plan - to be able to find fixed point values of M, q, X for given values of sigma, mu, gamma, for a fixed value of p
// but I think what we have at the moment is a map from z1, mu, gamma to the rest.
// not sure if this is the best way to find it but I can look at that later
// for multiple values of p can we use vector notation?

double phi(double z1)
{
	return (1.0 - erf(z1/sqrt(2.0)))/2.0;
}
// this becomes zero when z1 = infinity, which means sigma = 0

double second(double z1)
{
	return (((z1*z1) + 1.0)*(1.0 - erf(z1/sqrt(2.0)))/2.0) - (z1*exp(-z1*z1/2.0)/sqrt(2.0*M_PI));
	// this is the second order integral, this is never zero, im not sure if this is true!
}

double X(double z1)
{
	double pphi = phi(z1);
	return (gama*(p - 1.0)*pphi*pphi/second(z1)) + pphi;
}
// this is zero when phi = 0, so when sigma = 0, also for some gamma that probably is not in the range?

double help(double z1)
{
	return phi(z1)/X(z1);
}
// this is undefined for sigma = 0

double first(double z1)
{
	return (z1*(erf(z1/sqrt(2.0)) - 1.0)/2.0) + (exp(-z1*z1/2.0)/sqrt(2.0*M_PI));
}
// this is the first order integral, this can not be zero

double Mcoeff(double z1) // this is the coefficient of first order M in its polynomial
{
	//if (z1 == 1.0/0.0) return -1.0;
	//else 
	return z1*help(z1)/first(z1);
}

double Mgrad(double z1, double MM, double small)
{
	double MMcoeff = Mcoeff(z1);
	double grad = (MMcoeff*(MM + small)) + (mu*pow((MM + small), p - 1.0)) + k;
	grad -= (MMcoeff*(MM - small)) + (mu*pow((MM - small), p - 1.0)) + k;
	return grad/(2.0*small);
}

double M(double Mstart, double z1) //this is a solution to a p-order polynomial, so will probably do newton-raphson?
{
	double MMcoeff = Mcoeff(z1);
	double MM = Mstart; // maybe we dont need this? I think this should be set to k
	
	if (info == true)
	{
	ofstream mm; mm.open("zz.txt");
	ofstream ab; ab.open("abs.txt");

	cout << "now plotting function for M for z1 = " << z1 << endl;
	
	for (int i = 0; i <= 200000; i++)
	{
		double m = (0.00001*(double)i);
		mm << m;
		ab << (MMcoeff*m) + (mu*pow(m, p - 1.0)) + k;
		if (i != 200000) {mm << ", "; ab << ", ";}
	}
	mm.close(); ab.close();
	cout << "finished plotting M" << endl;
	}
	

	double small = pow(10.0, -6.0); // dont know about this?
	double y = (MMcoeff*MM) + (mu*pow(MM, p - 1.0)) + k;
	double grad;
	//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double Mnew;
	while (abs(y) > 0.0)
	{
		grad = Mgrad(z1, MM, small);
		if (grad == 0.0) break; //{cout << "grad 0 " << g(z1, zz + small) << ", " << g(z1, zz - small) << endl;}
		//if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		Mnew = MM - (y/grad);
		// if (z1 + zz < 0.0) cout << "badd" << endl;
		ynew = (MMcoeff*Mnew) + (mu*pow(Mnew, p - 1.0)) + k;
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		MM = Mnew;
		//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	
	//if (!(abs(sigma(z1, z2)) >= 0.0  && abs(gamma(z1, z2)) >= 0.0)) // found the wrong solution

	return MM;
}

double melp(double z1) // (k + mutot), maybe like help but for M
{
	double Mstart = k;
	double melp = k + (mu*pow(M(Mstart, z1), p - 1.0));
	if (melp == 0.0) cout << "melp = 0!! :(" << endl;
	return melp;
}
// this could be zero for a particular mu, but thats another problem

double sigtot(double z1)
{
	//if (z1 == 1.0/0.0) return 0.0;
	if (z1 == 0.0) cout << "z1 = 0 in sigtot function!! :(" << endl; // I dont know about this!
	return -melp(z1)/z1;
}
// this can be zero if melp is zero, but this can be undefined if z1 = 0, and thats probably a worse problem

double q(double z1)
{
	//if (z1 == 1.0/0.0) return pow(M(k, z1), 2.0);
	if (z1 == 0.0) cout << "z1 = 0 in q function!! :(" << endl;
	return second(z1)*pow(X(z1)*melp(z1)/(z1*phi(z1)), 2.0);
}

double sigma(double z1)
{
	//if (z1 == 1.0/0.0) return 0;
	return sqrt(2.0/(p*pow(q(z1), p - 1.0)))*sigtot(z1);
}
	

bool compare(ArrayXd i, ArrayXd j) {return (i(0) < j(0));}
void fiveplot(int grid)
{
	//string critstr = "HOcrit_" + to_string((int)p) + "_" + to_string((int)(10*mu)) + ".txt";
	//ofstream critfile; critfile.open(critstr);
	// so i think the range of gamma is between -1/(p-1) and 1, but I need to check this!!
	//
	
	//for (int gi = 0; gi <= 2; gi++)
	//{
		int gi = 2;
		
		string filename = "measfp_" + to_string((int)p) + "_" + to_string((int)(10*mu)) + "_" + to_string(gi) + ".txt"; // changed to hist for sigma != 0
		ofstream file; file.open(filename);


		gama = ((p - 2.0)*gi*gi/(2.0*(p - 1.0))) + ((4.0 - p)*gi/(2.0*(p - 1.0))) -1.0/(p - 1.0);
		cout << gama << endl;
		//critfile << crit(gamma)[0]; // i think this is critical sigma

		vector<ArrayXd> sigmat;
		double start = -40.0;
		double end = 40.0;
		double range = end - start;

		for (int i = 0; i <= grid; i++)
		{
			double z1 = (k - mu)/(((double)i*range/(double)grid) + start); // where did i get this from?
			if (M(k, z1) < 0.0) continue; //{cout << "bad solution 1 " << z1 << endl; continue;}
			if (q(z1) < 0.0) continue; //{cout << "bad solution 2 " << z1 <<endl; continue;}
			if (sigma(z1) < 0.0) continue; //{cout << "bad solution 3 " << z1 << endl; continue;}
			if (!(abs(sigma(z1)) >= 0.0)) continue; //{cout << "bad solution 4 " << z1 << endl; continue;}
			if (sigma(z1) > 10.0 || sigma(z1) < 0.1) continue;
			//if (sig(z1, z2) == 0.0) continue; // this is for histograms only

			ArrayXd entry = ArrayXd::Zero(6); // sigma, z1, phi, M, q, help
			entry(0) = sigma(z1);
			entry(1) = z1;
			entry(2) = phi(z1);
			entry(3) = M(k, z1);
			entry(4) = q(z1);
			entry(5) = help(z1);
			sigmat.push_back(entry);
			cout << "i " << i << ", z1 " << z1 << ", phi " << phi(z1) << ", M " << M(k, z1) << ", sig " << sigma(z1) << ", q " << q(z1) << endl;
			//cin >> carryon;
			//info = true;
		}
		sort(sigmat.begin(), sigmat.end(), compare);

		cout << gi << " " << sigmat.size() << " sorted!" << endl;
	
		for (int i = 0; i < sigmat.size(); i++)
		{
			file << sigmat[i].format(cc);
			//cout << sigmat[i].format(c) << endl;
			if (i != sigmat.size()-1) file << ",";
		}
		file.close();
		//if (gi != 2) critfile << ",";
	//}
	//critfile.close();
}
	


// what about if z1 = 0 and phi = 1/2 ??



int main()
{
	int grids = 400000; // for graphs with varying sigma
	p = 2.0;
	//gama;
	mu = 0.0;
	k = 1.0;
	
	fiveplot(grids);
	
	return 0;
}
/*
int main()
{
	p = 2.0;
	mu = 0.0;
	k = 1.0;
	gama = 1.0;
	double z1 = -0.745379;
	cout << phi(z1) << endl;
	cout << second(z1) << endl;
	cout << X(z1) << endl;
	cout << help(z1) << endl;
	cout << first(z1) << endl;
	cout << Mcoeff(z1) << endl;
	cout << M(k, z1) << endl;
	
	return 0;
}
*/
	




