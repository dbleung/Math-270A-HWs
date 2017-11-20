#include <vector>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "Math270A_DoubleArray2D.h"
#include "GridFun1D.h"

using namespace std;
using namespace Math270A;

#ifndef _GridFun2D_
#define _GridFun2D_

class GridFun2D
{
public:

	//Constructors
	GridFun2D() {
		initialize();
	}

	GridFun2D(const GridFun2D& G) {
		initialize(G);
	}

	GridFun2D(long xPanel, double xMin, double xMax, long yPanel, double yMin, double yMax) {
		initialize(xPanel, xMin, xMax, yPanel, yMin, yMax);
	}

	//Deconstructor
	virtual ~GridFun2D() {};

	//Initialize Functions
	void initialize() {
		xPanel = 0;
		xMin = 0;
		xMax = 0;

		yPanel = 0;
		yMin = 0; 
		yMax = 0;

		values.initialize();
	}

	void initialize(const GridFun2D& G) {
		xPanel = G.xPanel;
		xMin = G.xMin;
		xMax = G.xMax;
		hx = G.hx;

		yPanel = G.yPanel; 
		yMin = G.yMin;
		yMax = G.yMax;
		hy = G.hy;

		values = G.values;
	
	}

	void initialize(long xPanel, double xMin, double xMax, long yPanel, double yMin, double yMax) {
		this->xPanel = xPanel;
		this->xMin = xMin;
		this->xMax = xMax;
		this->hx = (xMax - xMin) / xPanel;

		this->yPanel = yPanel;
		this->yMin = yMin;
		this->yMax = yMax;
		this->hy = (yMax - yMin) / yPanel;
		DoubleArray2D values(xPanel + 1, yPanel + 1); // y goes down (m), x goes accross (n)
		this->values = values;

		GridFun1D u1Dx(xPanel,xMin,xMax);
		GridFun1D u1Dy(yPanel,yMin,yMax);

		this->u1Dx = u1Dx;
		this->u1Dy = u1Dy;

	}


	//Operator Overloading Functions
	void operator=(const GridFun2D& G)
	{
		double* dValues = values.getDataPointer();
		double* gValues = G.values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] = gValues[i];
		}
	}

	void operator+=(const GridFun2D& G)
	{
		double* dValues = values.getDataPointer();
		double* gValues = G.values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] += gValues[i];
		}
	}

	void operator-=(const GridFun2D& G)
	{
		double* dValues = values.getDataPointer();
		double* gValues = G.values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] -= gValues[i];
		}
	}

	void operator*=(double alpha)
	{
		double* dValues = values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] *= alpha;
		}
	}

	void operator/=(double alpha)
	{
		double* dValues = values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] /= alpha;
		}
	}

	//Other functions used...
	void setToValue(double d) {
		//sets values matrix to value d
		double* dValues = values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] =d;
		}
		
	}

	double normInf() {
		//Returns infinity norm of values matrix

		double maxVal = 0;
		double* dValues = values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			if (maxVal < dValues[i])
			maxVal = dValues[i];
		}
		return maxVal;
	}


	//Another operator overload function for << (used for cout)
	friend ostream& operator<<(ostream& outStream, const GridFun2D& G)
	{

		double* gValues = G.values.getDataPointer();
		for (long i = 0; i < G.values.getDataSize(); i++)
		{
			outStream << setw(5) << gValues[i] << " ";
			outStream << endl;
		}
	
		return outStream;
	}

	//Member Functions

	// Extracts 1D values from 2D values with constant y index
	//
	// Assumptions: u1Dx has been dimensioned commensurate
	// with the x-coordinate data of *this
	//
	void extractXslice(long yIndex, GridFun1D& u1Dx) const
	{
		for (long i = 0; i <= xPanel; i++)
		{
			u1Dx.values[i] = this->values(i, yIndex);
		}
	}
	// Inserts 1D values into 2D values with constant y index
	//
	// Assumptions: u1Dx has been dimensioned commensurate
	// with the x-coordinate data of *this
	//
	void insertXslice(long yIndex, const GridFun1D& u1Dx)
	{
		for (long i = 0; i <= xPanel; i++)
		{
			this->values(i, yIndex) = u1Dx.values[i];
		}
	}
	
	void extractYslice(long xIndex, GridFun1D& u1Dy) const
	{
		for (long j = 0; j <= yPanel; j++)
		{
			u1Dy.values[j] = this->values(xIndex,j);
		}
	}
	// Inserts 1D values into 2D values with constant y index
	//
	// Assumptions: u1Dx has been dimensioned commensurate
	// with the x-coordinate data of *this
	//
	void insertYslice(long xIndex, const GridFun1D& u1Dy)
	{
		for (long j = 0; j <= yPanel; j++)
		{
			this->values(xIndex, j) = u1Dy.values[j];
		}
	}


	//Internal Data
	double hx;
	double xMin;
	double xMax;
	long   xPanel;
	double hy;
	double yMin;
	double yMax;
	long yPanel;

	DoubleArray2D values;
	GridFun1D u1Dx;
	GridFun1D u1Dy;


};

#endif
