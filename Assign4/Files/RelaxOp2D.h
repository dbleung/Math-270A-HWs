#include <vector>
#include <cstdio>
#include <cmath>
#include "TriSolver.h"
#include "GridFun2D.h"
#include "GridFun1D.h"

using namespace std;

#ifndef _RelaxOp2D_
#define _RelaxOp2D_

class RelaxOp2D {
public:
	//
	//Functions to print arrays/vectors of data for debugging purposes
	//
	void print_array(const GridFun2D& G, double xPanel, double yPanel) {

		for (long j = yPanel; j >= 0; j--) {
			for (long i = 0; i <=xPanel; i++) {
				cout << G.values(i, j) << ",";
			}
			cout << endl;

		}
	}
	void print_vec(vector<double>& a, long M) {
		int i;
		for (i = 0; i < (M); i++) {
			cout << a[i] << ",";
		}
		i = M;
		cout << a[i] << endl;
	}

//
//Functions to compute second order difference in x or y direction for Douglas ADI
//

//
//Compute delta_y *dt *alphaY
//
void D2y(double dt, double alphaY, const GridFun2D& uIn, GridFun2D& d2yU) {
	double hy = (uIn.hy); 
	//cout << "hy: " << uIn.hy << endl;
	//cout << "index1: " << d2yU.values.getIndex1Size() << endl;
	//cout << "index2: " << d2yU.values.getIndex2Size() << endl;

	for (long i = 1; i < d2yU.values.getIndex1Size()-1; i++){
		for (long j = 1; j < d2yU.values.getIndex2Size()-1; j++) { //only for interior points, -1 because starts at 0, not 1
			//cout << "i=" << i << " j=" << j << endl;
			d2yU.values(i,j)= dt*alphaY*(uIn.values(i, j + 1) - 2 * uIn.values(i, j) + uIn.values(i, j - 1)) / pow(hy, 2);

		}
	}
}

//
//Compute delta_x*dt*alphaX
//
void D2x(double dt, double alphaX, const GridFun2D& uIn,  GridFun2D& d2xU) {
	double hx = (uIn.hx);
	//cout << "hx: " << uIn.hx << endl;
	//cout << "index1: " << d2xU.values.getIndex1Size() << endl;
	//cout << "index2: " << d2xU.values.getIndex2Size() << endl;
	for (long i = 1; i < d2xU.values.getIndex1Size()-1; i++) {
		for (long j = 1; j < d2xU.values.getIndex2Size()-1; j++) { //only for interior points
			//cout << "i=" << i << " j=" << j << endl;
			 d2xU.values(i,j) = dt*alphaX*(uIn.values(i + 1, j) - 2 * uIn.values(i, j) + uIn.values(i - 1, j)) / pow(hx, 2);
		}
	}
}


void initialize(double dt, double alphaX, double alphaY, const GridFun2D& a) {

	this->dt = dt;
	this->alphaX = alphaX;
	this->alphaY = alphaY;
	this->xPanel = a.xPanel;
	this->yPanel = a.yPanel;

    this->f.values = a.values; 
	this->f.xPanel = a.xPanel;
	this->f.xMax = a.xMax;
	this->f.xMin = a.xMin;
	this->f.yPanel = a.yPanel;
	this->f.yMax = a.yMax;
	this->f.yMin = a.yMin;
	this->f.hx = a.hx;
	this->f.hy = a.hy;

	GridFun2D d2xU(f.xPanel,f.xMin,f.xMax,f.yPanel,f.yMin,f.yMax);
	GridFun2D d2yU(f.xPanel, f.xMin, f.xMax, f.yPanel, f.yMin, f.yMax);

	this->d2yU.values = d2yU.values;
	this->d2yU.xPanel = d2yU.xPanel;
	this->d2yU.yPanel = d2yU.yPanel;

	this->d2xU.values = d2xU.values;
	this->d2xU.xPanel = d2xU.xPanel;
	this->d2xU.yPanel = d2xU.yPanel;


	GridFun2D Frhs(f.xPanel, f.xMin, f.xMax, f.yPanel, f.yMin, f.yMax);
	this->Frhs.values = Frhs.values;
	this->Frhs.xPanel = Frhs.xPanel;
	this->Frhs.yPanel = Frhs.yPanel;

	GridFun2D Fstar(f.xPanel, f.xMin, f.xMax, f.yPanel, f.yMin, f.yMax);
	Fstar = f;
	Fstar *= -dt;  //-f*dt
	
	this->Fstar.values = Fstar.values;
	this->Fstar.xPanel = Fstar.xPanel;
	this->Fstar.yPanel = Fstar.yPanel;
	
	GridFun1D f1DX(f.xPanel, f.xMin, f.xMax);
	GridFun1D u1DX(f.xPanel, f.xMin, f.xMax);
	GridFun1D f1DY(f.yPanel, f.yMin, f.yMax);
	GridFun1D u1DY(f.yPanel, f.yMin, f.yMax);

	this->f1DX.values = f1DX.values;
	this->f1DX.xPanel = f1DX.xPanel;
	this->u1DX.values = u1DX.values;


	this->f1DY.values = f1DY.values;
	this->f1DY.xPanel = f1DY.xPanel; //really should be yPanel...
	this->u1DY.values = u1DY.values;

	GridFun2D uOut(f.xPanel, f.xMin, f.xMax, f.yPanel, f.yMin, f.yMax);
	this-> uOut.values = uOut.values;

	//
	//GET LHS SIDE EQUATIONS for Trisolve (I-dt/2*alpha*delta) based on alpha, dt and implicit direction
	//

	//FOR first implicit direction
	vector<double> diag(f.xPanel + 1, 0.0);
	vector<double> upDiag(f.xPanel, 0.0);
	vector<double> loDiag(f.xPanel, 0.0);

	cout << "GET DIAGS for a step" << endl;
	setUpSystem(f.hx, xPanel, loDiag, diag, upDiag); //second order difference for loDiag, diag, upDiag
	
	for (long i = 0; i <= diag.size() - 1; i++) {
		diag[i] *= dt*alphaX / 2;
		diag[i] = 1 - diag[i];
	}
	
	diag[0] = 1;  //identity on boundaries
	diag[diag.size() - 1] = 1;  //identity on boundaries

	for (long i = 0; i <= upDiag.size() - 1; i++) {
		upDiag[i] *= -dt*alphaX / 2;
		loDiag[i] *= -dt*alphaX / 2;
	}

	upDiag[0] = 0;  //identity on boundaries
	loDiag[upDiag.size() - 1] = 0; //identity on boundaries

	cout << "Diag:" << endl;
	print_vec(diag, xPanel);
	cout << "loDiag:" << endl;
	print_vec(loDiag, xPanel - 1);
	cout << "upDiag:" << endl;
	print_vec(upDiag, xPanel - 1);

	this->diag = diag;
	this->loDiag = loDiag;
	this->upDiag = upDiag;

	//
	//Create new LHS for step b stored as diags (for 2nd implicit direction)
	//
	vector<double> diag2(f.yPanel + 1, 0.0);
	vector<double> upDiag2(f.yPanel, 0.0);
	vector<double> loDiag2(f.yPanel, 0.0);

	cout << "GET DIAGS for b step" << endl;
	setUpSystem(f.hy, yPanel, loDiag2, diag2, upDiag2); //second order difference for loDiag, diag, upDiag
	
	for (long i = 0; i <= diag2.size() - 1; i++) {
		diag2[i] *= dt*alphaY / 2;
		diag2[i] = 1 - diag2[i];
	}

	diag2[0] = 1; //identity on boundaries
	diag2[diag2.size() - 1] = 1;  //identity on boundaries

	for (long i = 0; i <= upDiag2.size() - 1; i++) {
		upDiag2[i] *= -dt*alphaY / 2;
		loDiag2[i] *= -dt*alphaY / 2;
	}

	upDiag2[0] = 0;  //identity on boundaries
	loDiag2[upDiag2.size() - 1] = 0;  //identity on boundaries

	cout << "Diag2:" << endl;
	print_vec(diag2, yPanel);
	cout << "loDiag2:" << endl;
	print_vec(loDiag2, yPanel - 1);
	cout << "upDiag2:" << endl;
	print_vec(upDiag2, yPanel - 1);

	this->diag2 = diag2;
	this->loDiag2 = loDiag2;
	this->upDiag2 = upDiag2;
}
//
//Computes the "A" Matrix (second order difference eq) for system, stored as loDiag, diag, and upDiag. A doesn't 
// change throughout apply member function

void setUpSystem(double h, long M,
	vector<double>& loDiag, vector<double>& diag, vector<double>& upDiag)
{
	loDiag.resize(M);
	upDiag.resize(M);
	diag.resize(M + 1);

	long i;

	// First equation 
	i = 0;
	diag[0] = -2.0 / (h*h);
	upDiag[0] = 1 / (h*h);

	
	// Interior equation coefficients associated with a standard
	// second order finite difference approximation
	for (long i = 1; i < M; i++)
	{
		loDiag[i - 1] = 1 / (h*h);
		upDiag[i] = 1 / (h*h);
		diag[i] = -2.0 / (h*h);
	}

	// Last equation
	i = M;
	diag[M] = -2.0 / (h*h);
	loDiag[M - 1] = 1 / (h*h);

}
//
//Apply ADI in 2D using modified form from Douglas's paper
//
void apply(const GridFun2D& uIn, GridFun2D& uOut) {
	/*
	cout << "uOut" << endl;
	print_array(uOut, f.xPanel, f.yPanel);
	cout << "f:" << endl;
	print_array(f, f.xPanel, f.yPanel);
	cout << "uIn: "  << endl;
	print_array(uIn, uIn.xPanel, uIn.yPanel);

	cout << "dt: " << dt << endl;
	cout << "alphaY: " << alphaY << endl;
	cout << "alphaX: " << alphaX << endl;
	*/
	D2y(dt, alphaY, uIn, d2yU); // d2yU = (dt*alphaY * Delta_y) uIn

	//cout << "d2yU: " <<  endl; //print out entire array
	//print_array(d2yU, d2yU.xPanel, d2yU.yPanel);


	D2x(dt*.5, alphaX, uIn, d2xU); // d2xU = (dt/2 * alphaX * Delta_x) uIn

	//cout << "d2xU: " << endl; //print out entire array
	//print_array(d2xU, d2xU.xPanel, d2xU.yPanel);


	//RHS of the first step (part a)
	Frhs = uIn; 
	Frhs += d2yU;
	Frhs += d2xU;
	Frhs += Fstar;              // Fstar = -dt*F

	//cout << "Fstar" << endl;
	//print_array(Fstar, Frhs.xPanel, Frhs.yPanel);

	//cout << "Frhs" << endl;
	//print_array(Frhs, Frhs.xPanel, Frhs.yPanel);

//
// Carry out x-coordinate implicit step for interior y-coordinate indices
// and the identity operator for j = 0 and j = yPanel indices.
//
// This loop overwrites the values of Frhs with the results, Frhs will hold u*


	triSolverX.initialize(xPanel+1, loDiag,diag,upDiag); //systemsize is same length as XPanel+1

	for (long j = 1; j < yPanel; j++)
	{
		Frhs.extractXslice(j, f1DX); //this stores jth row of Frhs as f1DX
		//cout << "f1DX" << endl;
		//print_vec(f1DX.values, f1DX.xPanel);
		//cout << "u1DX" << endl;
		//print_vec(u1DX.values, f1DX.xPanel);

		triSolverX.apply(f1DX.values, u1DX.values); //need to initialize loDiag, Diag, upDiag, etc for apply to work
		//cout << "u1DX" << endl;
		//print_vec(u1DX.values, f1DX.xPanel);
		Frhs.insertXslice(j, u1DX); //inserts calculated u1DX ub Frhs

	}

	// Identity operator for j = 0 and j = yPanel

	uIn.extractXslice(0, u1DX);
	Frhs.insertXslice(0, u1DX);

	uIn.extractXslice(yPanel, u1DX);
	Frhs.insertXslice(yPanel, u1DX);

	//cout << "u* values" << endl;
	//print_array(Frhs, Frhs.xPanel, Frhs.yPanel);

	//
	// Create right hand size for y-coordinate implicit step
	//
	d2yU *= -0.5;    // d2yU now  = -((dt/2)*alphaY * Delta_y) uIn
	Frhs += d2yU;    // Frhs now  = u*_(k+1)  - ((dt/2)*alphaY * Delta_y) uIn

	//cout << "Frhs for part b values" << endl;
	//print_array(Frhs, Frhs.xPanel, Frhs.yPanel);
		
	triSolverY.initialize(yPanel + 1, loDiag2, diag2, upDiag2);

 //
 // Carry out y-coordinate implicit step for interior x-coordinate indices.
 // and identity operator for i = 0 and i = xPanel
 //
 // This loop overwrites the values of uOut with the result
	
	
	for (long i = 1; i < xPanel; i++)
	{
		Frhs.extractYslice(i, f1DY);
		//cout << "f1DY" << endl;
		//print_vec(f1DY.values, f1DY.xPanel); //should be ypanel	

		triSolverY.apply(f1DY.values, u1DY.values);
		//cout << "u1DY" << endl;
		//print_vec(u1DY.values, f1DY.xPanel);// should be ypanel
		uOut.insertYslice(i, u1DY);
	

	}

	// Identity operator for i = 0 and i = xPanel

	Frhs.extractYslice(0, u1DY);
	uOut.insertYslice(0, u1DY);

	Frhs.extractYslice(xPanel, u1DY);
	uOut.insertYslice(xPanel, u1DY);

	//cout << "Final uOut:" << endl;
	//print_array(uOut, Frhs.xPanel, Frhs.yPanel);

	
}

// Internal Class Data Members
TriSolver triSolverX;
TriSolver triSolverY;

double dt; // still need to define these in initialize function
double alphaY;
double alphaX;

GridFun2D d2xU;
GridFun2D d2yU;
GridFun2D f;

GridFun2D Frhs; 
GridFun2D Fstar; 
GridFun2D uIn;
GridFun2D uOut;

GridFun1D u1DX;
GridFun1D u1DY;
GridFun1D f1DX;
GridFun1D f1DY;

long yPanel;
long xPanel;

vector<double> diag;
vector<double> upDiag;
vector<double> loDiag;

vector<double> diag2;
vector<double> upDiag2;
vector<double> loDiag2;
};
#endif

