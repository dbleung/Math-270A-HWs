#include <vector>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

#ifndef _GridFun1D_
#define _GridFun1D_

class GridFun1D
{
public:

	//Constructors

	//Null Constructor
	GridFun1D(){
		initialize();
	};

	//Copy Constructor
	GridFun1D(const GridFun1D& G){
		initialize(G);
	};

	//Initialize internal class data with data given by object in main
	GridFun1D(long xPanel, double xMin, double xMax){
		initialize(xPanel, xMin, xMax);
	};

	//Not sure what the constructors should initialize

	
	void initialize() {
		xPanel = 0;
		xMin = 0;
		xMax = 0;
		double hx = 0;
		vector<double> values(xPanel + 1, 0.0);
	}

	
	void initialize(const GridFun1D& G) {
		xPanel = G.xPanel;
		xMin = G.xMin;
		xMax = G.xMax;
		hx = G.hx;
		values = G.values;
		
	}


	void initialize(long xPanel, double xMin, double xMax) {
		this->xPanel = xPanel;
		this->xMin = xMin;
		this->xMax = xMax;
		this -> hx = (xMax-xMin)/xPanel;
		vector<double> values(xPanel+1,0.0);
		this -> values = values; 

	}

	// Operator Overloading Functions

	void operator=(const GridFun1D& G)
	{
		for (unsigned long i = 0; i < values.size(); i++)
		{
			values[i] = G.values[i];
		}
	}

	void operator+=(const GridFun1D& G)
	{
		for (unsigned long i = 0; i < values.size(); i++)
		{
			values[i] += G.values[i];
		}
	}

	void operator-=(const GridFun1D& G)
	{
		for (unsigned long i = 0; i < values.size(); i++)
		{
			values[i] -= G.values[i];
		}
	}

	void operator*=(double alpha)
	{
		for (unsigned long i = 0; i < values.size(); i++)
		{
			values[i] *= alpha;
		}
	}

	void operator/=(double alpha)
	{
		for (unsigned long i = 0; i < values.size(); i++)
		{
			values[i] /= alpha;
		}
	}

	// Other functions

	void setToValue(double d) {
		for (unsigned long i = 0; i < values.size(); i++)
		{
			values[i] = d;
		}
	}

	double normInf() {
		double maxval = 0;
		for (unsigned long i = 0; i < values.size(); i++)
		{
			if (maxval < values[i]) {
				maxval = values[i];
			}
		}
		return maxval;
	}

	//Another operator overload function for <<
	friend ostream& operator<<(ostream& outStream, const GridFun1D& G)
	{
		for (unsigned long i = 0; i < G.values.size(); i++)
		{
			outStream << setw(5) << G.values[i] << " ";
			outStream << endl;
		}
		return outStream;
	}



	// Internal Variables
	double             hx;
	double           xMin;
	double           xMax;
	long           xPanel;
	vector<double> values;

};

#endif
