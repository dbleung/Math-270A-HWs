#include <vector>
#include <cstdio>
#include <cmath>

using namespace std;

#ifndef _TriSolver_
#define _TriSolver_

class TriSolver
{
public:

	/// Null Constructor -- called when you declare an instance

	TriSolver()
	{
		initialize();
	};

	/// Copy constructor (creates duplicate of T)

	TriSolver(const TriSolver& T) // called when you declare an instance with an existing instance, & is pass be ref
	{
		initialize(T);
	};


	/**
	This constructor initializes the internal class data with entries
	of the n X n tridiagonal matrix T associated with the entries in
	the lower diagonal (loDiag), diagonal (diag), and upper diagonal (upDiag).
	*/

	TriSolver(long systemSize, vector<double>& loDiag, vector<double>& diag,
		vector<double>& upDiag)
	{
		initialize(systemSize, loDiag, diag, upDiag);
	}

	// Default destructor

	virtual ~TriSolver() {};

	/// Null initializer, resizes the internal arrays to zero length

	void initialize()
	{
		systemSize = 0;
		loDiag.clear();
		diag.clear();
		upDiag.clear();

	};

	///  Copy initializer. Duplicates the entries of T

	void initialize(const TriSolver& T) // note: & is pass by reference, changes original T made by constructor
	{
		systemSize = T.systemSize;
		loDiag = T.loDiag;
		diag = T.diag;
		upDiag = T.upDiag;

	};


	/**
	This initializer initializes the internal class data with entries
	of the n X n tridiagonal matrix T associated with the entries in
	the lower diagonal (loDiag), diagonal (diag), and upper diagonal (upDiag).
	*/

	void initialize(long systemSize, vector<double>& loDiag,
		vector<double>& diag, vector<double>& upDiag)
	{
		this->systemSize = systemSize; //"this" is a special pointer
		this->loDiag = loDiag;
		this->diag = diag;
		this->upDiag = upDiag;

	}

	//LU Factors Creates coefficients of LU factors of tridiagonal matrix given 
	//Values of the diagonals

	void LUFactors(vector<double>& c, vector<double>& gamma) {

		c[0] = diag[0];
		gamma[0] = upDiag[0] / c[0];

		for (long i = 1; i <= (systemSize-1); i++) {
			c[i] = diag[i] - loDiag[i - 1] * gamma[i - 1];
			if (i == systemSize-1) {
				break;
			}
			gamma[i] = upDiag[i] / c[i];
		}
	}

	//LUSolve Take LU factors from LUFactors function and solves tridiagonal system

	void LUSolve(vector<double>& f, vector<double>& u, vector<double>& c, vector<double>& gamma, vector<double>& v) {
		// Lv=f -> solve for v
		v[0] = f[0] / c[0];
		for (long j = 1; j <= (systemSize-1); j++) {
			v[j] = (f[j] - loDiag[j - 1] * v[j - 1]) / c[j];
		}

		// Uu=v -> solve for u
		u[systemSize-1] = v[systemSize-1];

		for (long k = systemSize - 2; k >= 0; k--) {
			u[k] = v[k] - gamma[k] * u[k + 1];

		}

	}



	/**
	The apply operator runs the tridiagonal solver for Tx=b to compute x using tridiagonal matrix T
	stored at vectors of the diagonals and vector b

	This uses two above functions: LUFactors and LUSolve

	Currently no bounds checking is performed.

	*/
	void apply(vector<double>& f, vector<double>& u)
	{
		vector<double> c(systemSize, 0.0);
		vector<double> gamma(systemSize-1, 0.0);
		vector<double> v(systemSize, 0.0);

		//LU Factorization (step 1) -> L stored as c on diag, loDiag below that, U stored as 1 on diag, gamma on upper diagonal
		LUFactors(c, gamma);
		//Take LU Factorization to solve for Lv=f and Uu=v (Step 2)
		LUSolve(f, u, c, gamma, v);



	}



	// Internal data
	long        systemSize;
	vector<double>  loDiag;
	vector<double>    diag;
	vector<double>  upDiag;

};


#endif