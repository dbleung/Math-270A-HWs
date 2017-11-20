#include <vector>
#include <cstdio>
#include <cmath>
#include "TriSolver.h"
#include "GridFun1D.h"
using namespace std;

#ifndef _RelaxOp1D_
#define _RelaxOp1D_



class RelaxOp1D
{
public:


	void print_array(vector<double>& a, long M) {
		int i;
		for (i = 0; i < (M); i++) {
			cout << a[i] << ",";
		}
		
		cout << endl;
	}

	//Computes the "A" Matrix for system, stored as loDiag, diag, and upDiag. A doesn't 
	// change throughout apply member function


	void setUpSystem(double h, double alpha, long M,
		vector<double>& loDiag, vector<double>& diag1, vector<double>& diag2, vector<double>& upDiag, double dt)
	{
		
		loDiag.resize(M);
		upDiag.resize(M);
		diag1.resize(M + 1);
		diag2.resize(M + 1);

		long i;

		// First equation is the identity

		i = 0;
		diag1[0] = 0.0;
		diag2[0] = 0.0;
		upDiag[0] = 0.0;

		// Interior equation coefficients associated with a standard
		// second order finite difference approximation

		for (long i = 1; i < M; i++)
		{
			loDiag[i - 1] = 1 / (h*h);
			upDiag[i] = 1/ (h*h);
			diag1[i] = -2.0 / (h*h);
			diag2[i] = -2.0 / (h*h);
		}

		// Last equation is the identity

		i = M;
		diag1[M] = 0.0;
		diag2[M] = 0.0;
		loDiag[M - 1] = 0.0;

	}



	void initialize(double dt, double alpha, const GridFun1D& f) {
		//Form  Astar = (I  + dt/2*A) stored as diagonals
		//f*=dt //ex: doesnt work

	
		vector<double> f1 = f.values;
		long M = f.xPanel;
		vector<double> loDiag(M, 0.0); //local copies
		vector<double> upDiag(M, 0.0);
		vector<double> diag1(M+1, 0.0);
		vector<double> diag2(M+1, 0.0);

		// Create Astar stored as lodiag, diag1, diag1 updiag (for LHS and RHS of equation)
		setUpSystem(f.hx, alpha, f.xPanel, loDiag, diag1, diag2, upDiag,dt); 

		for (long i = 0; i < M; i++) {
			loDiag[i] *= alpha*dt/2;
			upDiag[i] *= alpha*dt/2;
		}
		
		for (long i = 0; i < M+1; i++) { //is there a way to use an overloaded operator here?
			diag1[i] *= alpha*dt / 2;
			diag2[i] *= alpha*dt / 2;
			diag1[i] =  1- diag1[i];
			diag2[i] = 1 + diag2[i];
		}
		cout << "FINAL Astar VALUES" << endl;
		cout << "lodiag:" << endl;
		print_array(loDiag, M);
		cout << "updiag:" << endl;
		print_array(upDiag, M);
		cout << "diag1:" << endl;
		print_array(diag1, M+1);
		cout << "diag2:" << endl;
		print_array(diag2, M+1);

		

		for (long i = 0; i < M; i++) {
			f1[i] *= dt;	
		}

		vector<double> loDiag_n(M,0.0); // NEGATIVE VERSIONS OF loDIag and upDiag for LHS of System of eq
		vector<double> upDiag_n(M, 0.0);

		for (long i = 0; i < M; i++) {
			loDiag_n[i] = -loDiag[i]; //for this case when LHS has I-alpha*dt*A
			upDiag_n[i] = -upDiag[i];
		}


		this->loDiag = loDiag;
		this->upDiag = upDiag;
		this->loDiag_n = loDiag_n;
		this->upDiag_n = upDiag_n;
		this->diag1 = diag1;
		this->diag2 = diag2;
		this->f1 = f1;
		this->M = M;
		
	}

	void applyTriOp(vector<double>& vIn, vector<double>& vOut, long M, vector<double> diag)
	{

		long i;

		if (M == 1)
		{
			vOut[0] = diag[0] * vIn[0]; return;
		}

		vOut[0] = diag[0] * vIn[0] + upDiag[0] * vIn[1];

		for (i = 1; i < M- 1; i++)
		{
			vOut[i] = loDiag[i - 1] * vIn[i - 1] + diag[i] * vIn[i] + upDiag[i] * vIn[i + 1];
		}

		i = M - 1;
		vOut[i] = loDiag[i - 1] * vIn[i - 1] + diag[i] * vIn[i];
	}

	
	//Apply Tridiagonal Solver
	void apply(const GridFun1D& uIn, GridFun1D& uOut) {

		//Form Final right hand side, (I  + dt/2*A)* u_In - dt*f
	
		GridFun1D uIn1 = uIn;
		//cout << "uIn1:" << endl;
		//print_array(uIn1.values, M+1);

		vector<double> vOut(M+1, 0.0);

		// Apply forward operation (I + dt/2 *A) * u:
		applyTriOp(uIn1.values, vOut, M+1,diag2);
	    //	cout << "vOut after TriOP:" << endl;
	    //	print_array(vOut, M+1);

		//Subtract off dt*f (i.e. f1)
		for (long i = 0; i <= M; i++) {
			vOut[i] -= f1[i];
		}

		//cout << "Final RHS:" << endl;
		//print_array(vOut, M + 1);

	//	for (long i = 0; i <= M; i++) {
		//	uOut.values[i] = 0;
	//	}
		//cout << "Before TriSolver, uOut=" << endl;
		//print_array(uOut.values, M + 1);

		//Solve implicit equations 
		//applyTriSolver(vOut, uOut.values,loDiag,upDiag); 

		tSolver.initialize(M+1, loDiag_n, diag1,  upDiag_n);
		
		tSolver.apply(vOut,uOut.values);
	//	cout << "After TriSolver, uOut=" << endl;
	//	print_array(uOut.values, M + 1);

		//cout << "Check TriSolver, T*uOut = " << endl;
		//applyTriOp(uOut.values, vOut, M + 1,diag1);
		//print_array(vOut, M + 1);

	}


	// Internal Data:
	TriSolver     tSolver;  // instantiate a  inverse operator
	//TriOperator   tOperator; //forward operator
	vector<double>  loDiag;
	vector<double>   diag1;
	vector<double>   diag2;
	vector<double>  upDiag;
	vector<double> loDiag_n;
	vector<double> upDiag_n;
	vector<double>  f1;
	long M;
	};

#endif