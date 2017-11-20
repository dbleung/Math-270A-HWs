#include <vector>
#include <cstdio>
#include <cmath>

using namespace std;

#ifndef _Residual_
#define _Residual_

//
//This class calculates the residual for the target system Au = f for 
//
class Residual 
{
public:
	//
	//Constructors
	//
	Residual() {
		initialize();
	}

	//Copy Constructor
	Residual(const Residual& R) {
		initialize(R);
	}

	//Initialize internal data 
	Residual(vector<double> u, vector<double>& fapprox, double alpha, double M, double hx) {
		initialize(u,fapprox, alpha, M, hx);
	}

	//Initialize Functions
	void initialize() {
		M = 0;
		u.clear();
		alpha=0;
		fapprox.clear();
		hx = 0;
	}

	void initialize(const Residual& R) {
		u = R.u;
		alpha = R.alpha;
		M = R.M;
		fapprox = R.fapprox;
		hx = R.hx;
	}

	void initialize(vector<double> u, vector<double>& fapprox, double alpha, double M, double hx) {
		this->u = u; //"this" is a special pointer
		this->fapprox = fapprox;
		this->alpha = alpha;
		this->M = M;
		this->hx = hx;
	}

	void apply(vector<double> u, vector<double>& fapprox,double alpha, double M, double hx) {
		for (int i = 0; i <= M; i++) {
	
			if (i == 0 || i == M) {
				fapprox[i] = 1;
			}
			else {


				fapprox[i] = alpha*(u[i + 1] - 2 * u[i] + u[i - 1]) / pow(hx, 2);
			}
		

			
		}
	}

//Internal Data
	vector<double> fapprox;
	vector<double> u;
	double alpha;
	double M;
	double hx;
	
};


#endif