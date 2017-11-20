#include <cstdio>
#include <cmath>

#include "Math270A_DoubleArray2D.h"
using namespace std;
using namespace Math270A;


#include "GridFun2D.h"    // 2D grid function class
#include "RelaxOp2D.h"    // 2D Douglas relaxation operator class

#include "GNUplotUtility.h"

// cpp file to test the DoubleArray2D class created by prof anderson

void main() {
	double m = 2;
	double n = 2;
	DoubleArray2D values(m,n);
	
	double* dValues = values.getDataPointer();
	for (long i = 0; i < values.getDataSize(); i++)
	{
		dValues[i] = 1;
		cout << dValues[i] << endl;
		
	}

	

	// storage by rows but able to access by (i,j)?
	
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cout << values(i, j) << endl; //changing the pointer values changes values (pointers stored as dValues)

		}
		
	}

	cout << "GRID FUNCTION" << endl;
	long xPanel = 5;
	double xMin = 0;
	double xMax = 1;
	long yPanel = 2;
	double yMin = 0;
	double yMax = 1;

	GridFun2D grid(xPanel, xMin, xMax, yPanel, yMin, yMax);

	grid.setToValue(1.0);

	double* gridValues = grid.values.getDataPointer();
	for (long i = 0; i < grid.values.getDataSize(); i++)
	{
		//dValues[i] = 1;
		cout << gridValues[i] << endl;

	}
	cout << "HI" << endl;

	for (int i = 0; i < yPanel+1; i++) {
		for (int j = 0; j < xPanel+1; j++) {
			cout << grid.values(i, j) << endl; //changing the pointer values changes values (pointers stored as dValues)

		}

	}
	


}
