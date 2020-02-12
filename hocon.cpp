
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <fstream>
#include <omp.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iomanip>
using namespace Eigen;
using namespace std;
random_device generator;
mt19937 twist(generator());
//IOFormat p(20, DontAlignCols, ",", ",", "", "", "", "");
IOFormat c(6, DontAlignCols, ",", ",", "", "", "", "");
normal_distribution<double> distributionn(0.0, 1.0);
uniform_real_distribution<double> distributionr(0.0, 1.0);


int N;
double Nd; // the double version of N
int T;
double dt;
// these have index starting at 2
ArrayXd mu(4);
ArrayXd gama(4); // should this maybe be a Vector or a Tensor?
ArrayXd sigma(4); // its ok for now but may need to be different in numerical solution
int p;

ArrayXd functiom(MatrixXd m) {return m.array();}

template<typename T>
using  ArrayType = Eigen::Array<T, Eigen::Dynamic, 1>;

template<typename Scalar, typename sizeType>
auto Tarr(const Eigen::Tensor<Scalar, 1> &tensor,const sizeType rows) // this function maps tensor to array
{
    return Eigen::Map<const ArrayType<Scalar>> (tensor.data(), rows);
}

typedef Eigen::IndexPair<int> IP;
Eigen::array<IP, 1> d10 = { IP(1, 0) }; // second index of T

template<typename T>
using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar, typename sizeType>
auto Tmat(const Eigen::Tensor<Scalar, 2> &tensor,const sizeType rows, const sizeType cols)
{
    return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows, cols);
}

template<typename Scalar, typename... Dims>
auto Maten(const MatrixType<Scalar> &matrix, Dims... dims)
{
    constexpr int rank = sizeof... (Dims);
    return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), {dims...});
}

MatrixXd cholesky(int p)
{
	MatrixXd low = MatrixXd::Zero(p, p);

	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			double sum = 0.0;
			if (j == i)
			{
				for (int k = 0; k < j; k++) sum += pow(low(j, k), 2.0);
				low(j, j) = sqrt(1 - sum);
				//cout << "(" << i << "," << j << ") " << low(i,j) << endl;
			}
			else
			{
				for (int k = 0; k < j; k++) sum += (low(i, k)*low(j,k));
				low(i, j) = (gama(p-2) - sum)/low(j, j);
				//cout << "(" << i << "," << j << ") " << low(i,j) << endl;
			}
		}
	}

	return low;
}


Tensor<double, 2> order2()
{
	Tensor<double, 2> A(N, N);
	A.setZero();
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			MatrixXd cho = cholesky(2);
			VectorXd z(2);
			for (int zi = 0; zi < 2; zi++) z(zi) = distributionn(twist);
			MatrixXd list = cho*z;
			A(i,j) = (mu(0)/Nd) + ((sigma(0)*list(0))/sqrt(Nd));
			A(j,i) = (mu(0)/Nd) + ((sigma(0)*list(1))/sqrt(Nd));
		}
	}
	return A;
}
/*
Tensor<double, 2> order2() // do it like functionals
{
	MatrixXd AA(N,N);
	AA(0,0) = 0.0;
	for (int i = 1; i < N; i++)
	{
		AA(i,i) = 0.0;
		for (int j = 0; j < i; j++)
		{
			double z1 = distributionn(twist);
			double z2 = distributionn(twist);
			double y1 = (gama(0)*z1) + (sqrt(1.0-(gama(0)*gama(0)))*z2);
			AA(i,j) = (mu(0)/N) + ((sigma(0)*z1)/sqrt(N));
			AA(j,i) = (mu(0)/N) + ((sigma(0)*y1)/sqrt(N));
		}
	}
	return Maten(AA, N, N);
}
*/	
/*
Tensor<double, 3> order3()
{
	Tensor<double, 3> A(N, N, N);
	A.setZero();
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			for (int k = 0; k < j; k++)
			{
				MatrixXd cho = cholesky(3);
				VectorXd z(3);
				for (int zi = 0; zi < 3; zi++) z(zi) = distributionn(twist);
				MatrixXd list = cho*z;
				A(i,j,k) = (2.0*mu(1)/pow(Nd, 2.0)) + (sqrt(3.0)*sigma(1)*list(0)/Nd);
				A(j,k,i) = (2.0*mu(1)/pow(Nd, 2.0)) + (sqrt(3.0)*sigma(1)*list(1)/Nd);
				A(k,i,j) = (2.0*mu(1)/pow(Nd, 2.0)) + (sqrt(3.0)*sigma(1)*list(2)/Nd);
			}
		}
	}
	return A;
}

Tensor<double, 4> order4()
{
	cout << "aac" << endl;
	Tensor<double, 4> A(N, N, N, N);
	cout << "aa" << endl;
	A.setZero();
	cout << "ab" << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			for (int k = 0; k < j; k++)
			{
				for (int l = 0; l < k; l++)
				{
					MatrixXd cho = cholesky(4);
					VectorXd z(4);
					for (int zi = 0; zi < 4; zi++) z(zi) = distributionn(twist);
					MatrixXd list = cho*z;
					A(i,j,k,l) = (6.0*mu(2)/pow(Nd, 3.0)) + (sqrt(12.0/pow(Nd, 3.0))*sigma(2)*list(0));
					A(j,k,l,i) = (6.0*mu(2)/pow(Nd, 3.0)) + (sqrt(12.0/pow(Nd, 3.0))*sigma(2)*list(1));
					A(k,l,i,j) = (6.0*mu(2)/pow(Nd, 3.0)) + (sqrt(12.0/pow(Nd, 3.0))*sigma(2)*list(2));
					A(l,i,j,k) = (6.0*mu(2)/pow(Nd, 3.0)) + (sqrt(12.0/pow(Nd, 3.0))*sigma(2)*list(3));
				}
			}
		}
	}
	return A;
}
Tensor<double, 5> order5()
{
	Tensor<double, 5> A(N, N, N, N, N);
	A.setZero();
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			for (int k = 0; k < j; k++)
			{
				for (int l = 0; l < k; l++)
				{
					for (int m = 0; m < l; m++)
					{
						MatrixXd cho = cholesky(5);
						VectorXd z(5);
						for (int zi = 0; zi < 5; zi++) z(zi) = distributionn(twist);
						MatrixXd list = cho*z;
						A(i,j,k,l,m) = (24.0*mu(3)/pow(Nd, 4.0)) + (sqrt(60.0)*sigma(3)*list(0)/pow(Nd, 2.0));
						A(j,k,l,m,i) = (24.0*mu(3)/pow(Nd, 4.0)) + (sqrt(60.0)*sigma(3)*list(1)/pow(Nd, 2.0));
						A(k,l,m,i,j) = (24.0*mu(3)/pow(Nd, 4.0)) + (sqrt(60.0)*sigma(3)*list(2)/pow(Nd, 2.0));
						A(l,m,i,j,k) = (24.0*mu(3)/pow(Nd, 4.0)) + (sqrt(60.0)*sigma(3)*list(3)/pow(Nd, 2.0));
						A(m,i,j,k,l) = (24.0*mu(3)/pow(Nd, 4.0)) + (sqrt(60.0)*sigma(3)*list(4)/pow(Nd, 2.0));
					}
				}
			}
		}
	}
	return A;
}
		
*/

class simulation
{
	public:
	MatrixXd AF;
	ArrayXd xF; ArrayXd yF;
	bool fixedF; bool uniqueF; bool divergeF;
	vector<ArrayXd> trajxF; vector<ArrayXd> trajyF; //
	Tensor<double, 2> A; //Tensor<double, 3> B; Tensor<double, 4> C; Tensor<double, 5> D;
	Tensor<double, 1> x; Tensor<double, 1> y; Tensor<double, 1> k;
	bool fixed; bool unique; bool diverge;
	vector<ArrayXd> trajx; vector<ArrayXd> trajy; //
	simulation()
	{
		diverge = false;
		divergeF = false;
		x = Tensor<double, 1>(N); y = Tensor<double, 1>(N); k = Tensor<double, 1>(N);
		xF = ArrayXd(N); yF = ArrayXd(N);
		for (int i = 0; i < N; i++) {x(i) = distributionr(twist); xF(i) = x(i); y(i) = distributionr(twist); yF(i) = y(i); k(i) = 1.0;}
		//trajx.push_back(Tarr(x, N)); trajy.push_back(Tarr(y, N)); //
		//trajxF.push_back(xF); trajyF.push_back(yF); //
		A = order2(); //B = order3(); C = order4(); D = order5();
		//B.setZero(); C.setZero(); D.setZero();
		AF = Tmat(A, N, N);
	}
	~simulation() {}
	void check() // this is for traj only last 1% and this does each species induvidually
	{
		fixed = true; fixedF = true;
		for (int t = 2; t < (0.01*T)+2; t++)
		{
			if (((abs(trajx[t] - trajx[t-1]) < 0.0001) || (abs((trajx[t] - trajx[t-1])/(trajx[t-1] - trajx[t-2])) < 1.0)).minCoeff() == 0) {cout << "fixed " << endl; fixed = false; break;}
			if (((abs(trajy[t] - trajy[t-1]) < 0.0001) || (abs((trajy[t] - trajy[t-1])/(trajy[t-1] - trajy[t-2])) < 1.0)).minCoeff() == 0) {cout << "fixed " << endl; fixed = false; break;}
		}
		for (int t = 2; t < (0.01*T)+2; t++)
		{	
			if (((abs(trajxF[t] - trajxF[t-1]) < 0.0001) || (abs((trajxF[t] - trajxF[t-1])/(trajxF[t-1] - trajxF[t-2])) < 1.0)).minCoeff() == 0) {cout << "fixedF " << endl; fixedF = false; break;}
			if (((abs(trajyF[t] - trajyF[t-1]) < 0.0001) || (abs((trajyF[t] - trajyF[t-1])/(trajyF[t-1] - trajyF[t-2])) < 1.0)).minCoeff() == 0) {cout << "fixedF " << endl; fixedF = false; break;}
		}
		unique = true;
		uniqueF = true;
		for (int t = 2; t < (0.01*T)+2; t++)
		{
			if (((abs(trajx[t] - trajy[t]) < 0.0001) || ((abs(trajx[t] - trajy[t])/(trajx[t-1] - trajy[t-1])) < 1.0)).minCoeff() == 0) {cout << "unique " << t << endl; unique = false; break;}
		}
		for (int t = 2; t < (0.01*T)+2; t++)
		{
			if (((abs(trajxF[t] - trajyF[t]) < 0.0001) || ((abs(trajxF[t] - trajyF[t])/(trajxF[t-1] - trajyF[t-1])) < 1.0)).minCoeff() == 0) {cout << "uniqueF " << t << endl; uniqueF = false; break;}
		}
	}
	void run()
	{
		for (int t = 0; t < T; t++)
		{	
			//cout << t << endl;

			// should comment out a lot of this stuff when isolating orders

			// for x ////////

			Tensor<double, 1> k1 = dt*x*(k - x
			+ A.contract(x, d10) 
			//+ B.contract(x, d10).contract(x, d10)
			//+ C.contract(x, d10).contract(x, d10).contract(x, d10)
			//+ D.contract(x, d10).contract(x, d10).contract(x, d10).contract(x, d10)
			);
					
			Tensor<double, 1> k2 = dt*(x + (k1/2.0))*(k - (x + (k1/2.0))
			+ A.contract(x + (k1/2.0), d10)
			//+ B.contract(x + (k1/2.0), d10).contract(x + (k1/2.0), d10)
			//+ C.contract(x + (k1/2.0), d10).contract(x + (k1/2.0), d10).contract(x + (k1/2.0), d10)
			//+ D.contract(x + (k1/2.0), d10).contract(x + (k1/2.0), d10).contract(x + (k1/2.0), d10).contract(x + (k1/2.0), d10)
			);
	
			Tensor<double, 1> k3 = dt*(x + (k2/2.0))*(k - (x + (k2/2.0))
			+ A.contract(x + (k2/2.0), d10)
			//+ B.contract(x + (k2/2.0), d10).contract(x + (k2/2.0), d10)
			//+ C.contract(x + (k2/2.0), d10).contract(x + (k2/2.0), d10).contract(x + (k2/2.0), d10)
			//+ D.contract(x + (k2/2.0), d10).contract(x + (k2/2.0), d10).contract(x + (k2/2.0), d10).contract(x + (k2/2.0), d10)
			);
	
			Tensor<double, 1> k4 = dt*(x + k3)*(k - (x + k3)
			+ A.contract(x + k3, d10)
			//+ B.contract(x + k3, d10).contract(x + k3, d10)
			//+ C.contract(x + k3, d10).contract(x + k3, d10).contract(x + k3, d10)
			// + D.contract(x + k3, d10).contract(x + k3, d10).contract(x + k3, d10).contract(x + k3, d10)
			);
	
			x += (k1 + (2.0*k2) + (2.0*k3) + k4)/6.0;
	
				
			// for y /////////////
	
			k1 = dt*y*(k - y
			+ A.contract(y, d10) 
			//+ B.contract(y, d10).contract(y, d10)
			//+ C.contract(y, d10).contract(y, d10).contract(y, d10)
			//+ D.contract(y, d10).contract(y, d10).contract(y, d10).contract(y, d10)
			);
	
			k2 = dt*(y + (k1/2.0))*(k - (y + (k1/2.0))
			+ A.contract(y + (k1/2.0), d10)
			//+ B.contract(y + (k1/2.0), d10).contract(y + (k1/2.0), d10)
			//+ C.contract(y + (k1/2.0), d10).contract(y + (k1/2.0), d10).contract(y + (k1/2.0), d10)
			//+ D.contract(y + (k1/2.0), d10).contract(y + (k1/2.0), d10).contract(y + (k1/2.0), d10).contract(y + (k1/2.0), d10)
			);
	
			k3 = dt*(y + (k2/2.0))*(k - (y + (k2/2.0))
			+ A.contract(y + (k2/2.0), d10)
			//+ B.contract(y + (k2/2.0), d10).contract(y + (k2/2.0), d10)
			//+ C.contract(y + (k2/2.0), d10).contract(y + (k2/2.0), d10).contract(y + (k2/2.0), d10)
			//+ D.contract(y + (k2/2.0), d10).contract(y + (k2/2.0), d10).contract(y + (k2/2.0), d10).contract(y + (k2/2.0), d10)
			);
	
			k4 = dt*(y + k3)*(k - (y + k3)
			+ A.contract(y + k3, d10)
			//+ B.contract(y + k3, d10).contract(y + k3, d10)
			//+ C.contract(y + k3, d10).contract(y + k3, d10).contract(y + k3, d10)
			//+ D.contract(y + k3, d10).contract(y + k3, d10).contract(y + k3, d10).contract(y + k3, d10)
			);

			y += (k1 + (2.0*k2) + (2.0*k3) + k4)/6.0;

			//////////////////////////////////////////////
			
			
			ArrayXd k1F = dt*xF*(1.0 - xF + functiom(AF*xF.matrix())); //
			ArrayXd k2F = dt*(xF + (k1F/2.0))*(1.0 - (xF + (k1F/2.0)) + functiom(AF*(xF + (k1F/2.0)).matrix())); //
			ArrayXd k3F = dt*(xF + (k2F/2.0))*(1.0 - (xF + (k2F/2.0)) + functiom(AF*(xF + (k2F/2.0)).matrix())); //
			ArrayXd k4F = dt*(xF + k3F)*(1.0 - (xF + k3F) + functiom(AF*(xF + k3F).matrix())); //
			xF += (k1F + (2.0*k2F) + (2.0*k3F) + k4F)/6.0;

			k1F = dt*yF*(1.0 - yF + functiom(AF*yF.matrix())); //
			k2F = dt*(yF + (k1F/2.0))*(1.0 - (yF + (k1F/2.0)) + functiom(AF*(yF + (k1F/2.0)).matrix())); //
			k3F = dt*(yF + (k2F/2.0))*(1.0 - (yF + (k2F/2.0)) + functiom(AF*(yF + (k2F/2.0)).matrix())); //
			k4F = dt*(yF + k3F)*(1.0 - (yF + k3F) + functiom(AF*(yF + k3F).matrix())); //
			yF += (k1F + (2.0*k2F) + (2.0*k3F) + k4F)/6.0;


			/////////////////////////////////////////////
				
			if (Tarr(x, N).maxCoeff() > pow(10.0, 5.0) || Tarr(y, N).maxCoeff() > pow(10.0, 5.0)) {diverge = true; cout << "diverged at " << t << endl;}// break;}
			if (xF.maxCoeff() > pow(10.0, 5.0) || yF.maxCoeff() > pow(10.0, 5.0)) {divergeF = true; cout << "divergedF at " << t << endl;}// break;}
	
			// store data //////////////////

			if (t >= (0.99*T)-2) {trajx.push_back(Tarr(x, N)); trajy.push_back(Tarr(y, N));} //
			if (t >= (0.99*T)-2) {trajxF.push_back(xF); trajyF.push_back(yF);} //
			//trajx.push_back(Tarr(x, N)); trajy.push_back(Tarr(y, N));
			//trajxF.push_back(xF); trajyF.push_back(yF);
		}
		//if (diverge == false) check();
		check();
	}
};


		
				



void plot(int plots)
{
	ofstream datax;
	datax.open ("trajectoriesx.txt");
	ofstream datay;
	datay.open ("trajectoriesy.txt");
	
	ofstream dataxF;
	dataxF.open ("trajectoriesxF.txt");
	ofstream datayF;
	datayF.open ("trajectoriesyF.txt");
	
	simulation sim;
	sim.run();
	if (sim.diverge == false) cout << "tensor " << sim.fixed << sim.unique << endl;
	if (sim.divergeF == false) cout << "matrix " << sim.fixedF << sim.uniqueF << endl;
	
	//cout << sim.fixed << endl;
	//cout << sim.unique << endl;
		
	for (int i = 0; i < plots; i++)
	{
		for (int t = 0; t < sim.trajx.size(); t++)
		{
			//cout << "con " << t << endl;
			datax << sim.trajx[t](i);
			datay << sim.trajy[t](i);
			//cout << sim.trajx[t](i);
			//cout << sim.trajy[t](i);
			if (i != plots-1 || t != sim.trajx.size()-1) {datax << ","; datay << ",";}
		}
		for (int t = 0; t < sim.trajxF.size(); t++)
		{
			//cout << "F " << t << endl;
			dataxF << sim.trajxF[t](i);
			datayF << sim.trajyF[t](i);
			if (i != plots-1 || t != sim.trajxF.size()-1) {dataxF << ","; datayF << ",";}
		}
	}
}



int main(int argc, char** argv)
//int main()	
{
	if (argc != 3) {return 1;}
	gama = ArrayXd::Zero(4);
	mu = ArrayXd::Zero(4);
	sigma = ArrayXd::Zero(4);
	p = 2;
	N = 200;
	Nd = (double)N;
	T = 200000; //200000;
	dt = 0.001;
	//mu(p-2) = -2.0;
	int plots = 200;

	int v1 = atoi(argv[1]);
	int v2 = atoi(argv[2]);

	gama(p-2) = (0.02*(double)v1) - 1.0;
	double po = (0.04*(double)v2) - 1.0;
	sigma(p-2) = pow(10.0, po);

	plot(plots);
	
	//cout << cholesky(2) << endl;



	return 0;
}
/*	

	N = 200;
	Nd = (double)N;
	T = 200000;
	dt = 0.001;

	int order = 0; // this is the order we are looking at, 0 means order 2

	for (int ord = 0; ord < 4; ord ++) // sets everything to zero except for "order" sigma and gamma
	{
		mu(ord) = 0.0;
		if (ord != order) {sigma(ord) = 0.0; gama(ord) = 0.0;}
	}


	int v1 = atoi(argv[1]);
	int v2 = atoi(argv[2]);


	string filename = "hoA_" + to_string(v1) + "_" + to_string(v2) + ".txt"; // letter represents order
	ofstream file; file.open(filename);

	//file << v1 << ", " << v2 << endl;

	gama(order) = (0.02*(double)v1) - 1.0;
	double po = (0.04*(double)v2) - 1.0;
	sigma(order) = pow(10.0, po);


	int count0 = 0;
	int count1 = 0;
	int count2 = 0;

	for (int r = 0; r < 200; r ++)
	{
		cout << r << endl;
		simulation sim;
		sim.run();
		if (sim.diverge == false)
		{
			if (sim.fixed) {if (sim.unique) count0 ++; else count1 ++;}
			else count2 ++;
		}
		//cout << sim.fixed << sim.unique << endl;
	}

	//file << count0 << ", " << count1 << ", " << count2 << endl;

	file << (double)count0/200.0 << "," << (double)count1/200.0 << "," << (double)count2/200.0;

	file.close();


	return 0;
}
*/





