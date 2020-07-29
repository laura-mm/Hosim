
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
double p = 3.0;
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
	if (mu == 0.0) return -k*first(z1)/(z1*help(z1));
	else
	{
		//if (pow(z1*help(z1)/first(z1), 2.0) < 4.0*mu*k) return 1.0/0.0;
		//cout << "(z1*help(z1)/first(z1)) " << (z1*help(z1)/first(z1)) << endl;
		//cout << "sqrt(pow(z1*help(z1)/first(z1), 2.0) - (4.0*mu*k)))/(-2.0*mu) " << sqrt(pow(z1*help(z1)/first(z1), 2.0) - (4.0*mu*k))/(-2.0*mu) << endl;
		return ((z1*help(z1)/first(z1)) + sqrt(pow(z1*help(z1)/first(z1), 2.0) - (4.0*mu*k)))/(-2.0*mu);
	}

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
	//cout << "help " << help << endl;
	//cout << "sqrt " << sqrt(2.0/(p*pow(second(z1), p - 1.0))) << endl;
	//cout << "pow " << pow(first(z1)/M(k, z1), p - 2.0) << endl;
	
	return help(z1)*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/M(z1), p - 2.0); // this one!
	//return help(z1)*sqrt(2.0/(p*pow(second(z1), p - 1.0)));
}

double mud3(double z1) // for p = 3, value of mu where no solutions for M, not sure if divergence
{
	return pow(z1*help(z1)/first(z1), 2.0)/(4.0*k);
}

double otherM(double z1) // for p = 3 only
{
	return ((z1*help(z1)/first(z1)) - sqrt(pow(z1*help(z1)/first(z1), 2.0) - (4.0*mu*k)))/(-2.0*mu);
}

double othersig(double z1) // other sigma for p = 3
{
	return help(z1)*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/otherM(z1), p - 2.0);
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

double maxsig_z1(double zstart) // returns the z1 where sigma reaches a stationary point, depends on mu and gamma
{
	ofstream ab; ab.open("abs.txt");
	int grid = 1000000;
	double small = pow(10.0, -6.0);

	/*cout << "now plotting function for maxsig for mu = " << mu << " gamma = " << gama << endl;
	
	for (int i = 0; i <= grid; i++)
	{
		double z = 100.0*(double)i/(double)grid - 90.0;
		double sigrad = (sigma(z + small) - sigma(z - small))/(2.0*small);
		ab << z << "," << sigma(z) << "," << sigrad;
		if (i != grid) {ab << ",";}
	}
	ab.close();
	cout << "finished plotting z" << endl;*/
	
	double z1 = zstart;
	//double small = pow(10.0, -6.0);
	double y = (sigma(z1 + small) - sigma(z1 - small))/(2.0*small);
	double grad;
	//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = (sigma(z1 + small) + sigma(z1 - small) - (2.0*sigma(z1)))/(small*small);
		if (grad == 0.0) break;
		//if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		znew = z1 - (y/grad);
		ynew = (sigma(znew + small) - sigma(znew - small))/(2.0*small);
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		z1 = znew;
		//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	return z1;	
}

void fiveplot3(int grid, int mi)
{
	mu = muval[mi];
	
	double zcrit = crit_z1();
	//cout << "crit " << zcrit << endl;
	
	for (int gi = 0; gi < 6; gi++)
	{
		cout << mi << "," << gi << endl;
		
		string filename = "measfp3_" + to_string(mi) + "_" + to_string(gi)+ ".txt";
		ofstream file; file.open(filename);

		gama = gammaval[gi];
	
		double sig_prev = -1.0; // next sigma has to be bigger than previous

		for (int i = 0; i <= grid; i++)
		{
			double z1 = ((double)i*100.0/(double)grid) - 50.0;
			if (sigma(z1) < 0.0) continue; //{cout << "bad solution 3 " << z1 << endl; continue;}
			//if (!(abs(sigma(z1)) >= 0.0)) continue; //{cout << "bad solution 4 " << z1 << endl; continue;}
			if (sigma(z1) < sig_prev) break; // this one!
			sig_prev = sigma(z1);
			if (i != 0) file << ",";
			//cout << "i " << i << ", z1 " << z1 << ", phi " << phi(z1) << ", M " << entry(3) << ", sig " << sigma(z1) << ", q " << q(z1) << endl;
			file << sigma(z1) << "," << z1 << "," << phi(z1) << "," << M(z1) << "," << q(z1) << "," << help(z1);
		}
		file.close();
	}
}

void phase(int grid) // for gamma against sigma
{
	string filename = "phline_" + to_string((int)p) + "_" + to_string((int)(10*mu)) + ".txt";
	ofstream file; file.open(filename);
	double zcrit = crit_z1();
	for (int ig = 0; ig <= grid; ig++)
	{
		double gi = (double)ig/(double)grid;
		gama = (2.0*(p - 2.0)*gi*gi/(p - 1.0)) + ((4.0 - p)*gi/(p - 1.0)) - 1.0/(p - 1.0);
		file << gama << "," << sigma(zcrit);
		if (ig != grid) file << ",";
	}
	file.close();
}

void bunin3(int grid, int gi)
{
	gama = gammaval[gi];
	double critz1 = crit_z1();
	double zstart = 0.0;
	double zstop;
	
	ofstream points; points.open("points3_" + to_string(gi) + ".txt");
	
	ofstream fil; fil.open("mu3_" + to_string(gi) + ".txt");
	for (int i = 0; i <= grid; i++)
	{
		//cout << "i " << i << endl;
		mu = (4.0*(double)i/(double)grid) - 3.0;
		//cout << mu << endl;
		if (i != 0) fil << ",";
		zstart = maxsig_z1(zstart);
		if (zstart < critz1) fil << mu << "," << 1.0/0.0 << "," << sigma(zstart); // no crit and maxsig line
		else fil << mu << "," << sigma(critz1) << "," << sigma(zstart); // critical line and maxsig line
		//cout << "mu = " << mu << ", zstart = " << zstart << ", maxsig = " << sigma(zstart) << ", sigc = " << sigma(z1) << endl;
		for (int mi = 0; mi < 12; mi ++)
		{
			//if (abs(mu - muval[mi]) < pow(10.0, -15.0)) cout << mu << ", " << abs(mu - muval[mi]) << endl;
			if (abs(mu - muval[mi]) < pow(10.0, -15.0))
			{
				//cout << mu << endl;
				if (zstart < critz1) points << 1.0/0.0 << "," << sigma(zstart) << ",";
				else points << sigma(critz1) << "," << sigma(zstart) << ",";
				double zlow;
				for (int i = 0; i <= grid; i++)
				{
					double z1 = 100.0*(double)i/(double)grid - 50.0;
					if (sigma(z1) > 0.0 && sigma(z1) < 1.0/0.0) {zlow = z1; break;}
				}
				if (abs(mud3(zlow) - mu) < 0.001) points << sigma(zlow);
				else points << 1.0/0.0;
				if (mi != 11) points << ",";
			}
		}
				
		if (sigma(zstart) > 0.0) zstop = zstart;
		//cout << "sig(zstop) = " << sigma(zstop) << endl;
	}
	fil.close();
	points.close();
	//cout << "zstop = " << zstop << endl;
	ofstream file; file.open("bunin3_" + to_string(gi) + ".txt");
	for (int i = 0; i <= grid; i++)
	{
		double z1 = (50.0 + zstop)*(double)i/(double)grid - 50.0;
		mu = mud3(z1);
		double MM = (z1*help(z1)/first(z1))/(-2.0*mu);
		double sigmaa = help(z1)*sqrt(2.0/(3.0*pow(second(z1), 2.0)))*first(z1)/MM;
		if (!(sigmaa >= 0.0)) break; // this one!
		if (i != 0) file << ",";
		file << mu << "," << sigmaa;
		//cout << z1 << ", " << mu << "," << sigma(crit_z1()) << "," << sigma(z1) << endl;
	}
	file.close();
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
	
	muval = vector<double>{-0.25, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6, 1.0};
	gammaval = vector<double>{-0.5, -0.4, -0.2, 0.0, 0.4, 1.0};
	sigmaval = vector<double>(21);
	for (int i = 0; i <= 20; i++) sigmaval[i] = pow(10.0, (0.1*(double)i) - 1.5);
	
	for (int mi = 0; mi < 12; mi++) fiveplot3(grids, mi);
	//for (int gi = 0; gi < 6; gi++) bunin3(1000000, gi);
	//bunin3(10000, 1);
	

	
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




