// this is for multiple orders, multiple values of p
// this one is for p = 2 and 3
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
double k;
ArrayXd mu(4);
ArrayXd gama(4);
ArrayXd sigma(4);

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

double X(double z1, double help)
{
	return phi(z1)/help;
}
// this is zero when phi = 0, so when sigma = 0, also for some gamma that probably is not in the range?

double first(double z1)
{
	return (z1*(erf(z1/sqrt(2.0)) - 1.0)/2.0) + (exp(-z1*z1/2.0)/sqrt(2.0*M_PI));
}
// this is the first order integral, this can not be zero

double mutot(double M)
{
	double m = 0.0;
	for (int i = 0; i < 4; i++)
	{
		m += mu(i)*pow(M, i + 1);
	}
	return m;
}

double Mgrad(double z1, double help, double MM, double small)
{
	double Mcoeff = help*z1/first(z1);
	double grad = (Mcoeff*(MM + small)) + mutot(MM + small) + k;
	grad -= (Mcoeff*(MM - small)) + mutot(MM - small) + k;
	return grad/(2.0*small);
}

double M(double z1, double help) //this is a solution to a p-order polynomial, so will probably do newton-raphson?
{
	double Mcoeff = help*z1/first(z1); // may be better because doesnt rely on other parameters
	//double MM = Mstart; // maybe we dont need this? I think this should be set to k
	double MM = k;
	
	if (info == true)
	{
	ofstream mm; mm.open("zz.txt");
	ofstream ab; ab.open("abs.txt");

	cout << "now plotting function for M for z1 = " << z1 << endl;
	
	for (int i = 0; i <= 200000; i++)
	{
		double m = (0.00001*(double)i);
		mm << m;
		ab << (Mcoeff*m) + mutot(m) + k;
		if (i != 200000) {mm << ", "; ab << ", ";}
	}
	mm.close(); ab.close();
	cout << "finished plotting M" << endl;
	}
	

	double small = pow(10.0, -6.0); // dont know about this?
	double y = (Mcoeff*MM) + mutot(MM) + k;
	double grad;
	//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double Mnew;
	while (abs(y) > 0.0)
	{
		grad = Mgrad(z1, help, MM, small);
		if (grad == 0.0) {cout << "grad 0 " << endl; return 1.0/0.0;}
		//if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		Mnew = MM - (y/grad);
		ynew = (Mcoeff*Mnew) + mutot(Mnew) + k;
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		MM = Mnew;
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	
	//if (!(abs(sigma(z1, z2)) >= 0.0  && abs(gamma(z1, z2)) >= 0.0)) // found the wrong solution

	return MM;
}

double sigtot(double z1, double help)
{
	return help*M(z1, help)/first(z1);
}

double q(double z1, double help)
{
	return second(z1)*pow(M(z1, help)/first(z1), 2.0);
}

double gamatot(double z1, double help)
{
	return (1.0 - help)/X(z1, help);
}

// above here is general for all combinations of p, below is for p = 2 and 3 only

double sig3(double z1, double help, double sig2)
{
	return sqrt(2.0*(pow(sigtot(z1, help), 2.0) - (q(z1, help)*pow(sig2, 2.0)))/(3.0*pow(q(z1, help), 2.0)));
}

double sgrad(double z1, double help, double sig2, double small)
{
	double grad = pow((sig2 + small), 2.0) + (3.0*q(z1, help)*pow(sig3(z1, help, sig2 + small), 2.0)) - (pow(help, 2.0)/phi(z1));
	grad -= pow((sig2 - small), 2.0) + (3.0*q(z1, help)*pow(sig3(z1, help, sig2 - small), 2.0)) - (pow(help, 2.0)/phi(z1));
	return grad/(2.0*small);
}

double sig2_crit(double z1, double help) // this should output the critical value of of sigma_2 for given z1 and help
{	
	if (info == true)
	{
	ofstream s2; s2.open("zz.txt");
	ofstream ab; ab.open("abs.txt");

	cout << "now plotting function for critical sigma_2 for z1 = " << z1 << " help = " << help << endl;
	
	for (int i = 0; i <= 10000; i++)
	{
		double sig2 = (double)i*sigtot(z1, help)/(sqrt(q(z1, help))*10000.0);
		s2 << sig2;
		ab << pow(sig2, 2.0) + (3.0*q(z1, help)*pow(sig3(z1, help, sig2), 2.0)) - (pow(help, 2.0)/phi(z1));
		if (i != 10000) {s2 << ", "; ab << ", ";}
	}
	s2.close(); ab.close();
	cout << "finished plotting critical sigma" << endl;
	}
			
	double sig2 = 0.5*sigtot(z1, help)/sqrt(q(z1, help));
	double small = pow(10.0, -6.0); // dont know about this?
	double y = pow(sig2, 2.0) + (3.0*q(z1, help)*pow(sig3(z1, help, sig2), 2.0)) - (pow(help, 2.0)/phi(z1));
	double grad;
	//if (info == true) 
	{cout << "sig2 " << sig2 << ", y " << y << endl;}
	double ynew;
	double sig2new;
	while (abs(y) > 0.0)
	{
		grad = sgrad(z1, help, sig2, small);
		cout << "a" << endl;
		if (grad == 0.0) {cout << "grad 0 " << pow((sig2 + small), 2.0) + (3.0*q(z1, help)*pow(sig3(z1, help, sig2 + small), 2.0)) << ", " << pow((sig2 - small), 2.0) + (3.0*q(z1, help)*pow(sig3(z1, help, sig2 - small), 2.0)) << endl; break;}
		//if (info == true)
		{cout << "small " << small << " grad " << grad << endl;}
		cout << "b" << endl;
		sig2new = sig2 - (y/grad);
		cout << "sig2new " << sig2new << endl;
		cout << "c" << endl;
		if (sig2new > sigtot(z1, help)/sqrt(q(z1, help))) {sig2new = sigtot(z1, help)/sqrt(q(z1, help)) - small; cout << "upper bound" << endl;}
		cout << "d" << endl;
		if (sig2new < 0.0) {sig2new = 0.0; cout << "lower bound" << endl;}
		cout << "e" << endl;
		ynew = pow(sig2new, 2.0) + (3.0*q(z1, help)*pow(sig3(z1, help, sig2new), 2.0)) - (pow(help, 2.0)/phi(z1));
		if (sig2new == 0.0 && ynew < 0.0) return -1.0;
		cout << "f" << endl;
		cout << "ynew " << ynew << endl;
		cout << "g" << endl;
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) {cout << "broke" << endl; break;}
		cout << "h" << endl;
		y = ynew;
		cout << "i" << endl;
		sig2 = sig2new;
		cout << "j" << endl;
		//if (info == true) {cout << "sig2 " << sig2 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	cout << "k" << endl;
	return sig2;
}

// plan is to first try with mu and gama = 0
// I dont think this works!!

void cycle(int grid1, int grid2)
{
	
	double help = 1.0;

	for (int i = 0; i <= grid1; i++)
	{
		double z1 = ((double)i/(double)grid1) - 1.0;
		//double z1 = -0.2;
		
		cout << "z1 " << z1 << endl;
		cout << "phi " << phi(z1) << endl;
		cout << "first " << first(z1) << endl;
		cout << "second " << second(z1) << endl;
		cout << "X " << X(z1, help) << endl;
		cout << "gamatot " << gamatot(z1, help) << endl;
		cout << "M " << M(z1, help) << endl;
		cout << "mutot " << mutot(M(z1, help)) << endl;
		cout << "sistot " << sigtot(z1, help) << endl;
		cout << "q " << q(z1, help) << endl << endl;
		
		ofstream file; file.open("testing.txt");
		
		//sigma(0) = sig2_crit(z1, help);
		//sigma(1) = sqrt(2.0*(pow(sigtot(z1, help), 2.0) - (q(z1, help)*pow(sigma(0), 2.0)))/(3.0*pow(q(z1, help), 2.0)));
		
		//cout << sig2_crit(z1, help) << endl;
		double sig2crit;
		if (phi(z1) == second(z1)) 
		
		
		sig2crit = sig2_crit(z1, help);
		bool solution = (sig2crit >= 0.0);
		if (solution) cout << sig2crit << endl;
		else cout << "no solution" << endl;
		
		
		// now what are the gammas?
		
		/*for (int j = 0; j <= grid2; j++)
		{
			gama(0) = (2.0*(double)j/(double)grid2) - 1.0;
			gama(1) = (gamatot(z1, help) - (gama(0)*pow(sigma(0), 2.0)))/(3.0*q(z1, help)*pow(sigma(1), 2.0));
			cout << "gama_2 " << gama(0) << endl;
			cout << "gama_3 " << gama(1) << endl;
			if (abs(sigma(1)) >= 0)
			{
				if (j != 0) file << ",";
				file << sigma(0) << "," << sigma(1) << "," << diff;
				//cout << sigma(0) << "," << sigma(1) << "," << diff;
			}
		}
		file.close();*/
		
		cin >> carryon;
	}
}

int main()
{
	gama = ArrayXd::Zero(4);
	mu = ArrayXd::Zero(4);
	sigma = ArrayXd::Zero(4);
	k = 1.0;
	
	//mu(0) = -2.0;
	//mu(1) = -2.0;
	
	//cycle(10, 10000);
	
	info = true;
	cout << M(-0.0000001, 1.0) << endl;
	
	return 0;
}
	
	
	
