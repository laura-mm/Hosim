
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



double Mgrad(double z1, double MM, double small)
{
	double Mcoeff = help(z1)*z1/first(z1); // may be better because doesnt rely on other parameters
	double grad = (Mcoeff*(MM + small)) + (mu*pow((MM + small), p - 1.0)) + k;
	grad -= (Mcoeff*(MM - small)) + (mu*pow((MM - small), p - 1.0)) + k;
	return grad/(2.0*small);
}

double M(double Mstart, double z1) //this is a solution to a p-order polynomial, so will probably do newton-raphson?
{
	double Mcoeff = help(z1)*z1/first(z1); // may be better because doesnt rely on other parameters
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
		ab << (Mcoeff*m) + (mu*pow(m, p - 1.0)) + k;
		if (i != 200000) {mm << ", "; ab << ", ";}
	}
	mm.close(); ab.close();
	cout << "finished plotting M" << endl;
	}
	

	double small = pow(10.0, -6.0); // dont know about this?
	double y = (Mcoeff*MM) + (mu*pow(MM, p - 1.0)) + k;
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
		ynew = (Mcoeff*Mnew) + (mu*pow(Mnew, p - 1.0)) + k;
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		MM = Mnew;
		//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	
	//if (!(abs(sigma(z1, z2)) >= 0.0  && abs(gamma(z1, z2)) >= 0.0)) // found the wrong solution

	return MM;
}
/*
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
*/

double sigtot(double z1)
{
	return help(z1)*M(k, z1)/first(z1);
}

double q(double z1)
{
	//if (z1 == 1.0/0.0) return pow(M(k, z1), 2.0);
	//if (z1 == 0.0) cout << "z1 = 0 in q function!! :(" << endl;
	//return second(z1)*pow(X(z1)*melp(z1)/(z1*phi(z1)), 2.0);
	return second(z1)*pow(M(k, z1)/first(z1), 2.0);
}

double sigma(double z1)
{
	//if (z1 == 1.0/0.0) return 0;
	return sqrt(2.0/(p*pow(q(z1), p - 1.0)))*sigtot(z1);
	//return sqrt(2.0/(second(z1)*p*pow(q(z1), p - 2.0)))*help(z1);
}

///////////////////////////////// end of equation type stuff

double gradient(double z1, double small)
{
	double grad = second(z1 + small) - (phi(z1 + small)*(p - 1.0));
	grad -= (second(z1 - small) - (phi(z1 - small)*(p - 1.0)));
	return grad/(2.0*small);
}

double crit_z1() // depends on p only
{	
	double z1 = -10.0;
	double small = pow(10.0, -6.0);
	double y = second(z1) - (phi(z1)*(p - 1.0));
	double grad;
	if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = gradient(z1, small);
		if (grad == 0.0) break;
		if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		znew = z1 - (y/grad);
		ynew = second(znew) - (phi(znew)*(p - 1.0));
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		z1 = znew;
		if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	
	/*if (!(abs(sigma(z1, z2)) >= 0.0  && abs(gamma(z1, z2)) >= 0.0)) // found the wrong solution
	{
		small = pow(10.0, -6.0);
		grad = gradient(z1, z2, small);
		//z1 -= 2.0*y/grad;
		
		z1 *= -1.0;
		y = (pow(param3(z1, z2), 2.0)*q(z1, z2)) - phi(z1, z2);
		if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		
		while (abs(y) > 0.0) // assume there are only 2 roots
		{
			grad = gradient(z1, z2, small);
			if (grad == 0.0) break; //{cout << "grad 0 " << g(z1, zz + small) << ", " << g(z1, zz - small) << endl;}
			if (info == true) {cout << "small " << small << " grad " << grad << endl;}
			znew = z1 - (y/grad);
			// if (z1 + zz < 0.0) cout << "badd" << endl;
			ynew = (pow(param3(znew, z2), 2.0)*q(znew, z2)) - phi(znew, z2);
			if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
			y = ynew;
			z1 = znew;
			if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
			if (small > pow(10.0, -10.0)) small /= 10.0;
		}
	}*/
	
	return z1;	
}
	

bool compare(ArrayXd i, ArrayXd j) {return (i(0) < j(0));}
void fiveplot(int grid)
{
	string critstr = "crit_" + to_string((int)p) + "_" + to_string((int)(10*mu)) + ".txt";
	ofstream critfile; critfile.open(critstr);
	// so i think the range of gamma is between -1/(p-1) and 1, but I need to check this!!
	
	double zcrit = crit_z1();
	cout << "crit " << zcrit << endl;
	
	
	
	for (int gi = 0; gi <= 2; gi++)
	{
		//int gi = 1;
		
		string filename = "measfp_" + to_string((int)p) + "_" + to_string((int)(10*mu)) + "_" + to_string(gi) + ".txt";
		ofstream file; file.open(filename);

		gama = ((p - 2.0)*gi*gi/(2.0*(p - 1.0))) + ((4.0 - p)*gi/(2.0*(p - 1.0))) -1.0/(p - 1.0);
		cout << gama << endl;
		critfile << sigma(zcrit); // critical sigma

		vector<ArrayXd> sigmat;
		double start = -40.0;
		double end = 40.0;
		double range = end - start;
		
		double Mstart = k;
		double sig_prev = -1.0; // next sigma has to be bigger than previous

		for (int i = 0; i <= grid; i++)
		{
			//double z1 = (k - mu)/(((double)i*range/(double)grid) + start); // where did i get this from?
			double z1 = ((double)i*50.0/(double)grid) - 50.0;
			//if (M(Mstart, z1) < 0.0) continue; //{cout << "bad solution 1 " << z1 << endl; continue;}
			//if (q(z1) < 0.0) continue; //{cout << "bad solution 2 " << z1 <<endl; continue;}
			if (sigma(z1) < 0.0) continue; //{cout << "bad solution 3 " << z1 << endl; continue;}
			//if (!(abs(sigma(z1)) >= 0.0)) continue; //{cout << "bad solution 4 " << z1 << endl; continue;}
			if (sigma(z1) > 10.0) continue;
			if (sigma(z1) < sig_prev) break;
			//if (sig(z1, z2) == 0.0) continue; // this is for histograms only

			ArrayXd entry = ArrayXd::Zero(6); // sigma, z1, phi, M, q, help
			entry(0) = sigma(z1);
			entry(1) = z1;
			entry(2) = phi(z1);
			entry(3) = M(Mstart, z1);
			entry(4) = q(z1);
			entry(5) = help(z1);
			
			Mstart = entry(3);
			sig_prev = entry(0);
			
			sigmat.push_back(entry);
			cout << "i " << i << ", z1 " << z1 << ", phi " << phi(z1) << ", M " << entry(3) << ", sig " << sigma(z1) << ", q " << q(z1) << endl;
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
		if (gi != 2) critfile << ",";
	}
	critfile.close();
}

void phase(int grid)
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

/*
// this one is for testing
bool compare(ArrayXd i, ArrayXd j) {return (i(1) < j(1));}
void fiveplot(int grid)
{
	string critstr = "Tcrit_" + to_string((int)p) + "_" + to_string((int)(10*mu)) + ".txt";
	ofstream critfile; critfile.open(critstr);
	// so i think the range of gamma is between -1/(p-1) and 1, but I need to check this!!
	
	double zcrit = crit_z1();
	cout << "crit " << zcrit << endl;
	
	//for (int gi = 0; gi <= 2; gi++)
	//{
		int gi = 1;
		
		string filename = "Tmeasfp_" + to_string((int)p) + "_" + to_string((int)(10*mu)) + "_" + to_string(gi) + ".txt";
		ofstream file; file.open(filename);

		gama = ((p - 2.0)*gi*gi/(2.0*(p - 1.0))) + ((4.0 - p)*gi/(2.0*(p - 1.0))) -1.0/(p - 1.0);
		cout << gama << endl;
		critfile << sigma(zcrit); // critical sigma

		vector<ArrayXd> sigmat;
		double start = -40.0;
		double end = 40.0;
		double range = end - start;

		for (int i = 0; i <= grid; i++)
		{
			double z1 = 20.0*(double)i/(double)grid - 10.0; //(k - mu)/(((double)i*range/(double)grid) + start); // where did i get this from?
			//if (M(k, z1) < 0.0) continue; //{cout << "bad solution 1 " << z1 << endl; continue;}
			//if (q(z1) < 0.0) continue; //{cout << "bad solution 2 " << z1 <<endl; continue;}
			//if (sigma(z1) < 0.0) continue; //{cout << "bad solution 3 " << z1 << endl; continue;}
			//if (!(abs(sigma(z1)) >= 0.0)) continue; //{cout << "bad solution 4 " << z1 << endl; continue;}
			//if (sigma(z1) > 10.0 || sigma(z1) < 0.1) continue;
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
		if (gi != 2) critfile << ",";
	//}
	critfile.close();
}
*/

int main()
{
	int grids = 400000; //400000; // for graphs with varying sigma
	p = 2.0;
	//gama = -0.1;
	mu = -2.0;
	k = 1.0;
	
	//fiveplot(grids);
	
	phase(grids);
	
	info = false;
	
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




