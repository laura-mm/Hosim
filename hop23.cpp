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
	if (mu(1) == 0.0) return -k/((z1*help/first(z1)) + mu(0));
	else return ((z1*help/first(z1)) + mu(0) + sqrt(pow((z1*help/first(z1)) + mu(0), 2.0) - (4.0*mu(1)*k)))/(-2.0*mu(1));
	

	/*double small = pow(10.0, -6.0); // dont know about this?
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

	return MM;*/
}

double sigtot(double z1, double help)
{
	return help*M(z1, help)/first(z1);
}

double q(double z1, double help)
{
	return second(z1)*pow(M(z1, help)/first(z1), 2.0);
}

double sigmap(double z1, double help, double p) // for a single value of p
{
	return help*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/M(z1, help), p - 2.0); // this one!
}


double gamatot(double z1, double help)
{
	return (1.0 - help)/X(z1, help);
}

double crit_z1(double p) // depends on p only
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

double z2c;
double z3c;

double sig3(double z1, double help, double sig2) // for where there is a combination of p = 2 & 3
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
	//{cout << "sig2 " << sig2 << ", y " << y << endl;}
	double ynew;
	double sig2new;
	while (abs(y) > 0.0)
	{
		grad = sgrad(z1, help, sig2, small);
		//cout << "a" << endl;
		if (grad == 0.0) break; //{cout << "grad 0 " << pow((sig2 + small), 2.0) + (3.0*q(z1, help)*pow(sig3(z1, help, sig2 + small), 2.0)) << ", " << pow((sig2 - small), 2.0) + (3.0*q(z1, help)*pow(sig3(z1, help, sig2 - small), 2.0)) << endl; break;}
		//if (info == true)
		//{cout << "small " << small << " grad " << grad << endl;}
		//cout << "b" << endl;
		sig2new = sig2 - (y/grad);
		//cout << "sig2new " << sig2new << endl;
		//cout << "c" << endl;
		if (sig2new > sigtot(z1, help)/sqrt(q(z1, help))) {sig2new = sigtot(z1, help)/sqrt(q(z1, help)) - small;} // cout << "upper bound" << endl;}
		//cout << "d" << endl;
		if (sig2new < 0.0) {sig2new = 0.0;} // cout << "lower bound" << endl;}
		//cout << "e" << endl;
		ynew = pow(sig2new, 2.0) + (3.0*q(z1, help)*pow(sig3(z1, help, sig2new), 2.0)) - (pow(help, 2.0)/phi(z1));
		if (sig2new == 0.0 && ynew < 0.0) return -1.0;
		//cout << "f" << endl;
		//cout << "ynew " << ynew << endl;
		//cout << "g" << endl;
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break; //{cout << "broke" << endl; break;}
		//cout << "h" << endl;
		y = ynew;
		//cout << "y = " << y << endl;
		//cout << "i" << endl;
		sig2 = sig2new;
		//cout << "j" << endl;
		//if (info == true) {cout << "sig2 " << sig2 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	//cout << "k" << endl;
	return sig2;
}



// code to find divergent line through 3d space of mu2, mu3, sigtot
// help = 1, set z1, set mu2, start at 1 and work down? then from that find mu3?
// find critical mu3 from critical condition, then can find sigtot?
// cant really do much without having both mu2 and mu3

double maxsig_z1(double zstart, double help) // returns the z1 where sigma3 reaches a stationary point, depends on mu and gamma
{
	int grid = 1000000;
	double small = pow(10.0, -6.0);
	
	/*if (abs(mu(1)) < 0.001)
	{
	ofstream ab; ab.open("abs.txt");
	cout << "now plotting function for maxsig for mu = " << mu(1) << endl;
	
	for (int i = 0; i <= grid; i++)
	{
		double z = 100.0*(double)i/(double)grid - 90.0;
		double sigrad = (sigmap(z + small, help, 3.0) - sigmap(z - small, help, 3.0))/(2.0*small);
		ab << z << "," << sigmap(z, help, 3.0) << "," << sigrad;
		if (i != grid) {ab << ",";}
	}
	ab.close();
	cout << "finished plotting z" << endl;
	cin >> carryon;
	}*/
	
	
	double z1 = zstart;
	//double small = pow(10.0, -6.0);
	double y = (sigmap(z1 + small, help, 3.0) - sigmap(z1 - small, help, 3.0))/(2.0*small);
	double grad;
	//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = (sigmap(z1 + small, help, 3.0) + sigmap(z1 - small, help, 3.0) - (2.0*sigmap(z1, help, 3.0)))/(small*small);
		if (grad == 0.0) break;
		//if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		znew = z1 - (y/grad);
		ynew = (sigmap(znew + small, help, 3.0) - sigmap(znew - small, help, 3.0))/(2.0*small);
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		z1 = znew;
		//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	return z1;	
}

double noneg(double s) // if negitive then returns infinity
{
	if (s >= 0.0) return s;
	else return 1.0/0.0;
}

void bunin(int grid1, int grid2, int grid3)
{
	double help = 1.0;
	for (int i = 0; i <= grid1; i++)
	{
	
		//int i = 41;
		
		mu(0) = 10.0*(double)i/(double)grid1 - 4.0;
		double zstart = 0.0;
		double zstop;
		
		ofstream crit; crit.open("crit23_" + to_string(i) + ".txt");
		bool flag = false;
		for (int k = 0; k <= grid3; k++)
		{
			mu(1) = 11.0*(double)k/(double)grid3 - 4.0;
			//if (!(sigmap(z2c, help, 2.0) >= 0.0 && sigmap(z3c, help, 3.0) >= 0.0)) continue;
			//zstart = maxsig_z1(zstart);
			//fil << mu << "," << sigma(z1) << "," << sigma(zstart); // critical line and maxsig line
			//cout << "mu = " << mu << ", zstart = " << zstart << ", sigma = " << sigma(zstart) << endl;
			if (flag == true) crit << ",";
			if (sigmap(zstart, help, 3.0) > 0.0) zstart = maxsig_z1(zstart, help);
			//cout << zstart << endl;
			crit << mu(0) << "," << mu(1) << "," << noneg(sigmap(z2c, help, 2.0)) << "," << noneg(sigmap(z3c, help, 3.0)) << "," << noneg(sigmap(zstart, help, 3.0)) << "," << M(z2c, help);
			//if (i == 100 && k == 9000) cout << mu(0) << "," << mu(1) << "," << noneg(sigmap(z2c, help, 2.0)) << "," << noneg(sigmap(z3c, help, 3.0)) << "," << noneg(sigmap(zstart, help, 3.0)) << endl;
			if (sigmap(zstart, help, 3.0) > 0.0) zstop = zstart;
			flag = true;
		}
		crit.close();
		
		ofstream div; div.open("div23_" + to_string(i) + ".txt");
		flag = false;
		double intersection2 = 0.0;
		for (int j = 0; j <= grid2; j++)
		{
			//double z1 = (50.0 + max(z3c, zstop))*(double)j/(double)grid2 - 50.0;
			double z1 = 60.0*(double)j/(double)grid2 - 50.0;
			mu(1) = pow((z1*help/first(z1)) + mu(0), 2.0)/(4.0*k); // this is divergent condition
			//if (i == 41) cout << j << ", " << z1 << ", " << mu(1) << ", " << sigmap(z1, help, 2.0) << endl;
			//if (mu(1) > 6.0) continue;
			if (sigmap(z1, help, 2.0) >= 0.0 && sigmap(z1, help, 3.0) >= 0.0)
			{
				if (mu(1) > 7.0) break;
				if (flag == true) div << ",";
				if (abs(z2c - z1) < pow(10.0, -15.0)) {intersection2 = mu(1);} // cout << mu(0) << ", " << mu(1) << ", " << sigmap(z1, help, 3.0) << endl;}
				if (z1 > zstop) div << mu(0) << "," << mu(1) << "," << z1 << "," << sigtot(z1, help) << "," << sigmap(z1, help, 2.0) << "," << 1.0/0.0;
				else div << mu(0) << "," << mu(1) << "," << z1 << "," << sigtot(z1, help) << "," << sigmap(z1, help, 2.0) << "," << sigmap(z1, help, 3.0);
				//cout << mu(0) << "," << mu(1) << "," << z1 << "," << sigtot(z1, help) << "," << sigmap(z1, help, 2.0) << "," << sigmap(z1, help, 3.0) << endl;
				flag = true;
			}
			else
			{
				mu(1) = 0.0;
				if (sigmap(z1, help, 2.0) > 10.0 || sigmap(z1, help, 2.0) < 0.0) break;
				if (flag == true) div << ",";
				div << mu(0) << "," << mu(1) << "," << z1 << "," << sigtot(z1, help) << "," << sigmap(z1, help, 2.0) << "," << 1.0/0.0;
				flag = true;
			}	
		}
		div.close();
		
		mu(1) = intersection2;
		ofstream critp2; critp2.open("critp2_" + to_string(i) + ".txt");
		critp2 << mu(0) << "," << mu(1) << "," << sigmap(z2c, help, 2.0);
		critp2.close();
	}
}


			


/*
void critical(int grid1, int grid2)
{
	double help = 1.0;
	
	cout << "sig2 " << sigmap(z2c, help, 2.0) << endl;
	cout << "sig3 " << sigmap(z3c, help, 3.0) << endl;
	
	ofstream file; file.open("testing.txt");
	for (int i = 0; i <= grid1; i++)
	{
		double z1 = (z2c - z3c)*(double)i/(double)grid1 + z3c;
		
		cout << "z1 " << z1 << endl;
		cout << "phi " << phi(z1) << endl;
		cout << "first " << first(z1) << endl;
		cout << "second " << second(z1) << endl;
		cout << "X " << X(z1, help) << endl;
		cout << "gamatot " << gamatot(z1, help) << endl;
		cout << "M " << M(z1, help) << endl;
		cout << "mutot " << mutot(M(z1, help)) << endl;
		cout << "sigtot " << sigtot(z1, help) << endl;
		cout << "q " << q(z1, help) << endl << endl;
		
		
		
		//if (abs(z1 - z2c) < pow(10.0, -10.0)) {sigma(0) = help/sqrt(second(z1)); sigma(1) = 0.0;}
		//else if (abs(z1 - z3c) < pow(10.0, -10.0)) {sigma(0) = 0.0; sigma(1) = sig3(z1, help, sigma(0));}
		//else {sigma(0) = sig2_crit(z1, help); sigma(1) = sig3(z1, help, sigma(0));}
		
		
		
		//sigma(1) = sqrt(2.0*(pow(sigtot(z1, help), 2.0) - (q(z1, help)*pow(sigma(0), 2.0)))/(3.0*pow(q(z1, help), 2.0)));
		
		//cout << "sig2crit " << sigma(0) << endl;
		//cout << "sig3crit " << sigma(1) << endl;
		
		double sig2crit;
		if (phi(z1) == second(z1)) {sigma(0) = sigmap(z1, help, 2.0); sigma(1) = 0.0;} // this will change for other parameters
		else
		{
			sig2crit = sig2_crit(z1, help);
			bool solution = (sig2crit >= 0.0);
			if (solution) {sigma(0) = sig2crit; sigma(1) = sig3(z1, help, sigma(0));}
			//else cout << "no solution" << endl;
		}
		
		//cout << sigma(0) << ", " << sigma(1) << endl << endl;
		
		file << z1 << "," << sigma(0) << "," << sigma(1) << "," << M(z1, help);
		if (i != grid1) file << ",";
		
		
		
		// now what are the gammas?
		
		for (int j = 0; j <= grid2; j++)
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
		file.close();
		
		//cin >> carryon;
	}
	file.close();
}*/

void testing(int grid1)
{
	double help = 1.0;
	
	cout << "sig2 " << sigmap(z2c, help, 2.0) << endl;
	cout << "sig3 " << sigmap(z3c, help, 3.0) << endl;
	
	ofstream file; file.open("testing.txt");
	for (int i = 0; i <= grid1; i++)
	{
		//double z1 = (z2c - z3c)*(double)i/(double)grid1 + z3c;
		double z1 = 100.0*(double)i/(double)grid1 - 50.0;
		//if (z1 < 0.0 && z1 > -1.0){
		//cout << endl << z1 << ", " << z1*help/first(z1) << ", " << sigmap(z1, help, 2.0) << endl;
		// given a value of z1 between the two single p critial points, find sig2 and sig3 at the critical point
		file << z1 << "," << help*sqrt((2.0/second(z1)) - (1.0/phi(z1))) << "," << (help*first(z1)/M(z1, help))*sqrt((2.0/(3.0*second(z1)))*((1.0/phi(z1)) - (1.0/second(z1)))) << "," << M(z1, help);
		//cout << z1 << "," << help*sqrt((2.0/second(z1)) - (1.0/phi(z1))) << "," << (help*first(z1)/M(z1, help))*sqrt((2.0/(3.0*second(z1)))*((1.0/phi(z1)) - (1.0/second(z1)))) << "," << M(z1, help) << endl;
		if (i != grid1) file << ",";
	}
	file.close();
}

void s3test(int grid1, int grid2, int grid3) // investigate M with mu2 and mu3 for sigma3 = 0
{
	double help = 1.0; // both gamma = 0
	
	
	
	//for (int i = 0; i <= grid1; i++)
	//{
		
		//mu(0) = 10.0*(double)i/(double)grid1 - 4.0;
		//cout << "i " << i << endl;
		
		mu(0) = 2.0;
		//for (int k = 0; k <= grid3; k++)
		//{
			//ofstream file; file.open("s2test_" + to_string(i) + "_" + to_string(k) + ".txt");
			ofstream file; file.open("s2test.txt");
			//mu(1) = 11.0*(double)k/(double)grid3 - 4.0;
			mu(1) = -1.0;
			//cout << "k " << k << endl;
			for (int j = 0; j <= grid2; j++)
			{
				//cout << "j " << j << endl;
				double z1 = 60.0*(double)j/(double)grid2 - 50.0;
				//if (!(sigmap(z1, help, 2.0) >= 0.0)) continue;
				file << z1 << "," << mu(0) << "," << mu(1) << "," << sigmap(z1, help, 3.0) << "," << M(z1, help);
				//((z1*help/first(z1)) + mu(0) + sqrt(pow((z1*help/first(z1)) + mu(0), 2.0) - (4.0*mu(1)*k)))/(-2.0*mu(1));
				if (j != grid2) file << ",";
			}
			
			file.close();
		//}
	//}
}
	


int main()
{
	gama = ArrayXd::Zero(4);
	mu = ArrayXd::Zero(4);
	sigma = ArrayXd::Zero(4);
	k = 1.0;
	z2c = crit_z1(2.0);
	z3c = crit_z1(3.0);
	cout << "z2c " << z2c << endl;
	cout << "z3c " << z3c << endl;
	
	

	
	mu(0) = -3.0;
	mu(1) = 2.5;
	
	//cout << 1.0*sqrt(2.0/(3.0*pow(second(z3c), 2.0))) << ", " << pow(first(z3c)/M(z3c, 1.0), 1.0) << endl;
	
	
	testing(1000000); // for critical point between the two orders
	
	//bunin(100, 12000, 10000);
	
	//s3test(10, 1000000, 10);
	
	/*double help = 1.0;
	mu(0) = -4.0;
	double zstart = 0.0;
	double zstop;
	int grid3 = 100000;
	
	for (int k = 0; k <= grid3; k++)
	{
			mu(1) = 10.0*(double)k/(double)grid3 - 4.0;
			//if (!(sigmap(z2c, help, 2.0) >= 0.0 && sigmap(z3c, help, 3.0) >= 0.0)) continue;
			//zstart = maxsig_z1(zstart);
			//fil << mu << "," << sigma(z1) << "," << sigma(zstart); // critical line and maxsig line
			//cout << "mu = " << mu << ", zstart = " << zstart << ", sigma = " << sigma(zstart) << endl;
			//if (flag == true) crit << ",";
			zstart = maxsig_z1(zstart, help);
			//cout << zstart << endl;
			//crit << mu(0) << "," << mu(1) << "," << sigmap(z2c, help, 2.0) << "," << sigmap(z3c, help, 3.0) << "," << sigmap(zstart, help, 3.0);
			//cout << mu(0) << "," << mu(1) << "," << sigmap(z2c, help, 2.0) << "," << sigmap(z3c, help, 3.0) << "," << sigmap(zstart, help, 3.0) << endl;
			if (sigmap(zstart, help, 3.0) > 0.0) zstop = zstart;
			//flag = true;
	}*/
	
	//info = true;
	//cout << M(-0.0000001, 1.0) << endl;
	
	return 0;
}
	
	
	
