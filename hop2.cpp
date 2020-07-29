
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
double p = 2.0;
double gama;
double mu;
double k;

vector<double> muval;
vector<double> gammaval;
vector<double> sigmaval;

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
	//return phi(z1)/X(z1);
	double sec = second(z1);
	return sec/(sec + (gama*(p - 1.0)*phi(z1))); // may be better because doesnt rely on other parameters
}
// this is undefined for sigma = 0

double first(double z1)
{
	return (z1*(erf(z1/sqrt(2.0)) - 1.0)/2.0) + (exp(-z1*z1/2.0)/sqrt(2.0*M_PI));
}
// this is the first order integral, this can not be zero


double M(double z1) //this is a solution to a p-order polynomial, so will probably do newton-raphson?
{
	return -k/((z1*help(z1)/first(z1)) + mu); // this could give -infinity

}


double sigtot(double z1)
{
	return help(z1)*M(z1)/first(z1);
}

double q(double z1)
{
	//if (z1 == 1.0/0.0) return pow(M(k, z1), 2.0);
	//if (z1 == 0.0) cout << "z1 = 0 in q function!! :(" << endl;
	//return second(z1)*pow(X(z1)*melp(z1)/(z1*phi(z1)), 2.0);
	return second(z1)*pow(M(z1)/first(z1), 2.0);
}

double sigma(double z1)
{
	//if (z1 == 1.0/0.0) return 0;
	//return sqrt(2.0/(p*pow(q(z1), p - 1.0)))*sigtot(z1);
	//return sqrt(2.0/(second(z1)*p*pow(q(z1), p - 2.0)))*help(z1);
	
	
	return help(z1)*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/M(z1), p - 2.0); // this one!
	//return help(z1)*sqrt(2.0/(p*pow(second(z1), p - 1.0)));
}

double mud2(double z1) // for p = 2 this is the value of mu where we find divergence
{
	return -z1*help(z1)/first(z1);
}


///////////////////////////////// end of equation type stuff


double crit_z1() // depends on p only
{	
	double z1 = -10.0;
	double small = pow(10.0, -6.0);
	double y = second(z1) - (phi(z1)*(p - 1.0));
	double grad;
	//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = ((second(z1 + small) - (phi(z1 + small)*(p - 1.0))) - (second(z1 - small) - (phi(z1 - small)*(p - 1.0))))/(2.0*small);
		if (grad == 0.0) break;
		//if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		znew = z1 - (y/grad);
		ynew = second(znew) - (phi(znew)*(p - 1.0));
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		z1 = znew;
		//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	return z1;	
}


double div2_z1() // for p = 2, depends on mu
{
	if (info == true)
	{
	ofstream ab; ab.open("abs.txt");

	cout << "now plotting function for zdiv for mu = " << mu << endl;
	
	for (int i = 0; i <= 400000; i++)
	{
		double z = (10.0*(double)i/400000.0) - 5.0;
		ab << z << "," << z*help(z)/first(z) + mu;
		if (i != 400000) {ab << ", ";}
	}
	ab.close();
	cout << "finished plotting z" << endl;
	}
	
	
	double z1 = 0.0;
	double small = pow(10.0, -6.0);
	double y = z1*help(z1)/first(z1) + mu;
	double grad;
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		grad = (((z1 + small)*help(z1 + small)/first(z1 + small)) - ((z1 - small)*help(z1 - small)/first(z1 - small)))/(2.0*small);
		//cout << "grad " << grad << endl;
		//cout << "grad1 " << ((z1 + small)*help(z1 + small)/first(z1 + small)) << endl;
		//cout << "z1 + small " << z1 + small << endl;
		//cout << "help(z1 + small) " << help(z1 + small) << endl;
		//cout << "first(z1 + small) " << first(z1 + small) << endl;
		if (grad == 0) {if (abs(y) < pow(10.0, -5.0)) return z1; else return -1.0/0.0;}
		znew = z1 - (y/grad);
		ynew = znew*help(znew)/first(znew) + mu;
		if (info == true) {cout << "znew " << znew << ", ynew " << ynew << endl;}
		//cout << "abs(ynew) " << abs(ynew) << ", abs(y) " << abs(y) << endl;
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break; //{cout << "thisun" << endl; break;}
		y = ynew;
		z1 = znew;
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	if (gama == -0.5 && mu == 1.5) return -1.0/0.0;
	else return z1;
}
	


void fiveplot2(int grid, int mi)
{
	mu = muval[mi];

	
	double zcrit = crit_z1();
	//cout << "crit " << zcrit << endl;
	
	for (int gi = 0; gi < 6; gi++)
	{
		cout << mi << "," << gi << endl;
		
		string filename = "measfp2_" + to_string(mi) + "_" + to_string(gi) + ".txt";
		ofstream file; file.open(filename);

		gama = gammaval[gi];
	

		for (int i = 0; i <= grid; i++)
		{
			double z1 = ((double)i*100.0/(double)grid) - 50.0;
			if (sigma(z1) < 0.0) continue; //{cout << "bad solution 3 " << z1 << endl; continue;}
			if (i != 0) file << ",";
			//cout << "i " << i << ", z1 " << z1 << ", phi " << phi(z1) << ", M " << entry(3) << ", sig " << sigma(z1) << ", q " << q(z1) << endl;
			if (M(z1) >= 0.0) file << sigma(z1) << "," << z1 << "," << phi(z1) << "," << M(z1) << "," << q(z1) << "," << help(z1);
			else file << sigma(z1) << "," << z1 << "," << phi(z1) << "," << 1.0/0.0 << "," << q(z1) << "," << help(z1);
		}
		file.close();
	}
}





void bunin2(int grid, int gi)
{
	gama = gammaval[gi];
	double critz1 = crit_z1();
	double critstop = mud2(critz1);
	
	ofstream fil; fil.open("mu2_" + to_string(gi) + ".txt");
	for (int i = 0; i <= grid; i++)
	{
		mu = ((min(critstop, 2.0) + 5.0)*(double)i/(double)grid) - 5.0;
		if (i != 0) fil << ",";
		fil << mu << "," << sigma(critz1); // critical line
		//cout << "mu = " << mu << ", zstart = " << zstart << ", sigma = " << sigma(zstart) << endl;
	}
	fil.close();
	
	ofstream file; file.open("bunin2_" + to_string(gi) + ".txt");
	for (int i = 0; i <= grid; i++)
	{
		double z1 = 100.0*(double)i/(double)grid - 50.0;
		mu = mud2(z1);
		if (!(sigma(z1) >= 0.0)) break; // this one!
		if (i != 0) file << ",";
		file << mu << "," << sigma(z1);
		//cout << z1 << ", " << mu << "," << sigma(crit_z1()) << "," << sigma(z1) << endl;
	}
	file.close();
	
	ofstream points; points.open("points2_" + to_string(gi) + ".txt");
	for (int mi = 0; mi < 12; mi ++)
	{
		//if (mi == 11) info = true;
		mu = muval[mi];
		if (mu <= critstop) points << sigma(critz1) << "," << sigma(div2_z1()) << ",";
		else points << 1.0/0.0 << "," << sigma(div2_z1()) << ",";
		double zlow;
		for (int i = 0; i <= grid; i++)
		{
			double z1 = 100.0*(double)i/(double)grid - 50.0;
			if (M(z1) > 0.0 && M(z1) < 1.0/0.0) {zlow = z1; break;}
		}
		//cout << "mu = " << mu << ", zlow = " << zlow << endl;
		if (abs(mud2(zlow) - mu) < 0.001) points << sigma(zlow);
		else points << 1.0/0.0;
		if (mi != 11) points << ",";
	}
	points.close();
}


void testing(int grids)
{
	ofstream file; file.open("testing.txt");
	cout << crit_z1() << endl;
	double prev = 1.0;
	for (int i = 0; i <= grids; i ++)
	{
		double z1 = (10.0*double(i)/(double)grids) - 9.0;
		file << z1 << "," << sigma(z1) << "," << pow(z1*help(z1)/first(z1), 2.0)/(4.0*k);
		//cout << z1 << ", " << -z1/first(z1) << endl;
		//if (z1/second(z1) > prev) cout << z1 << ", " << prev << endl << endl;
		//prev = z1/second(z1);
		if (i != grids) file << ",";
	}
	file.close();
}
		

int main()
{
	int grids = 400000; //400000; // for graphs with varying sigma
	//gama = -0.0;
	//mu = 0.0;
	k = 1.0;
	
	muval = vector<double>{-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5};
	gammaval = vector<double>{-1.0, -0.8, -0.5, 0.0, 0.5, 1.0};
	sigmaval = vector<double>(21);
	for (int i = 0; i <= 20; i++) sigmaval[i] = pow(10.0, (0.1*(double)i) - 1.0);
	
	//for (int mi = 0; mi < 12; mi++) fiveplot2(grids, mi);
	for (int gi = 0; gi < 6; gi++) bunin2(1000000, gi);
	//bunin2(1000000, 2);
	
	return 0;
}
/*
int main()
{
	info = true;
	p = 3.0;
	mu = 0.0;
	k = 1.0;
	gama = 0.0;
	double z1 = -0.83;
	cout << phi(z1) << endl;
	cout << second(z1) << endl;
	cout << X(z1) << endl;
	cout << help(z1) << endl;
	cout << first(z1) << endl;
	cout << M(k, z1) << endl;
	cout << q(z1) << endl;
	cout << sigma(z1) << endl;
	
	return 0;
}
*/




