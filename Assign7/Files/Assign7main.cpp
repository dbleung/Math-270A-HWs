#include <cstdio>
#include <cmath>
#include <string>

using namespace std;

//
// Math 270A : Assignment 7
//
// Test program for variable time step relaxation to a solution of
//
// alphaX * d^2 U/dx^2  alphaY * d^2 U/dy^2 = f(x,y)  with Dirichlet boundary conditions
//
// for (x,y) in [xMin, xMax]X[yMin,yMax] using the 2D version of Douglas ADI
//
//
// Date: Oct. 24, 2017
//

#include "GridFun2D.h"    // 2D grid function class
#include "RelaxOp2D.h"    // 2D Douglas relaxation operator class
#include "GNUplotUtility.h" //used to output files readable by GNUplot
#include "ClockIt.h" // used to time functions

#ifdef _OPENMP
#include <omp.h>              // required for OpenMP
#endif


class TestProblem2D
{
public:

	TestProblem2D(double alphaX, double waveNumberX,double xMin, double xMax,
			      double alphaY, double waveNumberY,double yMin, double yMax)
	{
		this->alphaX      = alphaX;
		this->alphaY      = alphaY;
		this->waveNumberX = waveNumberX;
		this->waveNumberY = waveNumberY;
		this->pi          = 3.1415926535897932;
	}

	double operator()(double x, double y)
	{
		double d2Xfact = (2.0*pi*waveNumberX)*(2.0*pi*waveNumberX);
		double d2Yfact = (2.0*pi*waveNumberY)*(2.0*pi*waveNumberY);
		return -(cos(2.0*pi*x*waveNumberX)*cos(2.0*pi*y*waveNumberY))/(alphaX*d2Xfact + alphaY*d2Yfact);
	}

	double normInf()
	{
		double d2Xfact = (2.0*pi*waveNumberX)*(2.0*pi*waveNumberX);
		double d2Yfact = (2.0*pi*waveNumberY)*(2.0*pi*waveNumberY);
		return 1.0/(alphaX*d2Xfact + alphaY*d2Yfact);
	}

	double pi;
	double alphaX;
	double alphaY;
	double waveNumberX;
	double waveNumberY;
};

//
// This routine computes the residual of the five-point discretization of Poisson's
// equation with Dirichlet boundary conditions at interior points of a rectangular domain.
//
// Input : Approximate solution u

//         Right hand side f which specifies the forcing function at interior grid points
//         (boundary values of f are ignored)
//
// Output : The residual obtained by evaluating
//
// r = f - (alpha_x*Delta_x + alpha_y*Delta_y) u   at interior points
//
// r = 0                                           at boundary points
//

void evaluatePoissonResidual(double alphaX, double alphaY, GridFun2D& u, GridFun2D& f, GridFun2D& residual)
{
	residual.initialize(f);   // Initialize residual with f
	residual.setToValue(0.0); // Set values to 0.0

	// Residual at interior points

	double hy = u.hy;
	double hx = u.hx;
	for(long i = 1; i <  u.xPanel; i++)
	{
	for(long j = 1; j <  u.yPanel; j++)
	{
	    residual.values(i,j)  = f.values(i,j) - (alphaX*((u.values(i+1,j) - 2.0*u.values(i,j) + u.values(i-1,j))/(hx*hx))
				                        + alphaY*((u.values(i,j-1) - 2.0*u.values(i,j) + u.values(i,j+1))/(hy*hy)));
	}
    }
}


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

// Function to find the variable timestep dt based on current iteration value of j and initial grid size ho=[b - a]/M
double find_dt(int j, double ho, double alpha, double pi) {

	return 4 * pow(ho*pow(2, j), 2) / (alpha*pi*pi);
}

int main()
{

	//
	//Initialize number of threads and number of panels (user input)
	//

	long               xPanel=40;
	long               yPanel=40;
	long repetitions = 1;
#ifdef _OPENMP
	int threadCount = -1;
	cout << "Enter in number of threads: ";
	cin >> threadCount;

	if (threadCount > omp_get_max_threads()) { threadCount = omp_get_max_threads(); }
	if (threadCount <= 0) { threadCount = omp_get_max_threads(); }
	omp_set_num_threads(threadCount);

	printf("\n");
	printf("#############\n");
	printf("############# Using OpenMP With %d Threads\n", omp_get_max_threads());
	printf("#############\n");
	printf("\n");



	cout << "Enter number of panel in x direction:  ";
	cin >> xPanel;

	cout << "Enter number of panel in y direction:  ";
	cin >> yPanel;

	//cout << "Enter in number of repetitions: ";
	//cin >> repetitions;

#endif
	ClockIt clock_1;
	clock_1.start();
	// Set up test problem

	double alphaX      =  2.0;  // Laplace operator-x prefactor
	double alphaY      =  2.0;  // Laplace operator-y prefactor
	double waveNumberX =  1.0;  // test problem x-coordinate wave number
	double waveNumberY =  1.0;  // test problem y-coordinate wave number

	double dt     =     .1;  // Relaxation timestep

	//long   xPanel     =   40;  // X panel count
	//long   yPanel     =   40;  // Y panel count
	double tol    =  1.0e-6;  // Stopping tolerance

	double xMin   =     0.0;
	double xMax   =     1.0;
	double hx     = (xMax-xMin)/(double)xPanel;

	double yMin   =     0.0;
	double yMax   =     1.0;
	double hy     = (yMax-yMin)/(double)yPanel;
	double pi = 3.1415926535897932;

	int M;
	double ho;

	if (xPanel > yPanel) {
		M = xPanel;
		ho = hx;
	}
	else {
		M = yPanel;
		ho = hy;
	}

	int jstar = int(log2(M)) + 1; // step count for variable timestep (nearest integer of log2(M)), M should be greater number of xPanel or yPanel


	// List of files output

	vector<string> outputFileList;

	// Echo input parameters

	cout << "XXXX   Douglas ADI Program Start      XXXX " << endl;
	cout << "alpha_x : " << alphaX << endl;
	cout << "alpha_y : " << alphaY << endl;
	cout << "X Panel Count : " << xPanel << endl;
	cout << "Y Panel Count : " << yPanel << endl;
	cout << "Timestep      : " << dt << endl;

	//
	// Instantiate the test problem solution
	//

	TestProblem2D     testSoln(alphaX,waveNumberX,xMin,xMax,
                               alphaY,waveNumberY,yMin,yMax);

	double solnNorm = testSoln.normInf();

	//
	// Instantiate 2D grid functions
	//

	GridFun2D f(xPanel,xMin,xMax,yPanel,yMin,yMax);
	GridFun2D uk(xPanel,xMin,xMax,yPanel,yMin,yMax);
	GridFun2D ukp1(xPanel,xMin,xMax,yPanel,yMin,yMax);
	GridFun2D ukf(xPanel, xMin, xMax, yPanel, yMin, yMax);
	GridFun2D ukfp1(xPanel, xMin, xMax, yPanel, yMin, yMax);

	GridFun2D vTmp(xPanel, xMin, xMax,yPanel,yMin,yMax); // for variable time stepping
	GridFun2D vTmpp1(xPanel, xMin, xMax, yPanel, yMin, yMax); // for variable time stepping, plus 1 term
	GridFun2D vBarTmp(xPanel, xMin, xMax, yPanel, yMin, yMax);// for variable time stepping
	GridFun2D vBarTmpm1(xPanel, xMin, xMax, yPanel, yMin, yMax);// for variable time stepping, minus 1 term

	GridFun2D uExact(xPanel,xMin,xMax,yPanel,yMin,yMax);
	GridFun2D uSolnErr(xPanel,xMin,xMax,yPanel,yMin,yMax);
	GridFun2D ufSolnErr(xPanel, xMin, xMax, yPanel, yMin, yMax); // for fixed time step
	GridFun2D residual(xPanel,xMin,xMax,yPanel,yMin,yMax);
	GridFun2D residualf(xPanel, xMin, xMax, yPanel, yMin, yMax); // for fixed time step
	vector<double> values1d;
	GridFun1D dtTemp(jstar, xMin, xMax,values1d);//for variable time stepping, keeps values of dt for forward and backward sweep

	vector<RelaxOp2D> relaxOps(jstar+1); //check jstar sweep number?


   
	//
	// Construct the right hand side and the exact solution
	//

	double x; double y;
	
	for(long i = 0; i <= xPanel; i++)
	{
		x = xMin + i*hx;
		for(long j = 0; j <= yPanel; j++)
		{
			y =  yMin + j*hy;
			f.values(i,j)      = cos(2.0*pi*x*waveNumberX)*cos(2.0*pi*y*waveNumberY);
			uExact.values(i,j) = testSoln(x,y);
	    }
	}

	//
	// Set initial iterate with interior values equal to zero
	// and perimeter (boundary) values set the exact solution.
	//
	// Set boundary values of f to zero. (when j =0, and when i =0)

	uk.setToValue(0.0);

	long i; long j;
	for(i = 0; i <= xPanel; i++)
	{
		j   = 0;
		x   =   xMin + i*hx;
	    y   =   yMin + j*hy;
	    uk.values(i,j) = testSoln(x,y);
	    f.values(i,j)  = 0.0;

	    j   = yPanel;
		x   = xMin + i*hx;
	    y   = yMin + j*hy;
	    uk.values(i,j) = testSoln(x,y);
	    f.values(i,j)  = 0.0;
	}

	for(long j = 0; j <= yPanel; j++)
	{
		i = 0;
        x = xMin + i*hx;
		y = yMin + j*hy;

		uk.values(i,j) = testSoln(x,y);
		f.values(i,j)  = 0.0;

		i = xPanel;
		x = xMin + i*hx;
		y = yMin + j*hy;

		uk.values(i,j) = testSoln(x,y);
		f.values(i,j)  = 0.0;
	}

	// Output right hand side, initial data and exact solution

	GNUplotUtility::output(f,"f.dat");
	outputFileList.push_back("f.dat      : Forcing function");

	GNUplotUtility::output(uk,"u0.dat");
	outputFileList.push_back("u0.dat     : Starting iterate");

	GNUplotUtility::output(uExact,"uExact.dat");
	outputFileList.push_back("uExact.dat : Exact solution");

	//
	// Instantiate and initialize the relaxation operator
	//

	//
	//Get dt values and initialize relaxops for variable time stepping
	//
	for (int j = 0; j <= jstar; j++) {
		dtTemp.values[j] = find_dt(j, ho, alphaX, pi);//alphaX or alphay is fine
		relaxOps[j].initialize(dtTemp.values[j], alphaX, alphaY, f);
	}

    //
	//Initialize RelaxOperator for fixed time step
	//
	RelaxOp2D relaxOp;
	relaxOp.initialize(dt,alphaX,alphaY,f);

	//Fixed Time step Stuff
	ukf = uk;
	ukfp1=ukp1;

	clock_1.stop(); //stop timing

	double initializeTime = clock_1.getMilliSecElapsedTime();


	//
	// Relaxation loop
	//

	ClockIt clock_2;
	clock_2.start(); //start timing

	double diffNorm = 2*tol;
	long   iterMax  = 10000;
	long   iter     = 0;
	double uNorm;
	double relaxsteps = 0;



	while ((diffNorm > tol) && (iter < iterMax))
	{
		//relaxOp.initialize(dt, alphaX, alphaY, f);
		cout << "******ITERATION " << iter + 1 << "**************" << endl;

		//initialize forward sweep
		vTmp = uk;

		// Go through forward sweep 

		for (int j = 0; j <= (jstar); j++) {
			//relaxOp.initialize(dtTemp.values[j], alphaX,alphaY, f);
			relaxOps[j].apply(vTmp, vTmpp1);
			vTmp = vTmpp1;
			relaxsteps++;
		} // end for loop

		// Go through backward sweep
		vBarTmp = vTmp;
		for (int j = (jstar); j >= 0; j--) {
			//relaxOp.initialize(dtTemp.values[j], alphaX, alphaY, f);
			relaxOps[j].apply(vBarTmp, vBarTmpm1);
			vBarTmp = vBarTmpm1;
			relaxsteps++;
		} // end for loop


		ukp1 = vBarTmp;


		// Check scaled relative difference between iterates

		uNorm = ukp1.normInf();

		uk -= ukp1;
		diffNorm = uk.normInf();
		diffNorm /= (dt*uNorm);

		// Update iterate

		uk = ukp1;
		iter++;

		if (iter % 100 == 0) cout << "Iteration : " << iter << endl;
	}// end while loop

	clock_2.stop(); // stop timing

	double relaxTime = clock_2.getMilliSecElapsedTime();

	// Evaluate the residual and the realtive norm of the residual

	double residualNormInf;


	evaluatePoissonResidual(alphaX,  alphaY, uk, f, residual);

	residualNormInf = residual.normInf();
	residualNormInf /= f.normInf();

	GNUplotUtility::output(residual,"residual.dat");
	outputFileList.push_back("Residual   : residual.dat");


	// Evaluate the error with respect to the continuous solution


	uSolnErr =         uk;
	uSolnErr -=    uExact;

	GNUplotUtility::output(uk,"uK.dat");
	outputFileList.push_back("uK.dat     : Computed solution");

	GNUplotUtility::output(uSolnErr,"uSolnErr.dat");
	outputFileList.push_back("uSolnErr.dat : Computed solution error");

	double absSolnError    = uSolnErr.normInf();
	double relSolnError    = absSolnError/solnNorm;

	//
	// Print out the results using standard C I/0
	//
	relaxsteps = int(relaxsteps);
    printf("\nXXXX   Variable Timestep Relaxation Output XXXX \n\n");
    printf("Number of iterations               : %ld \n",iter);
	printf("Number of relaxation steps         : %5e \n", relaxsteps);
    printf("Scaled Difference between iterates : %10.5e \n",diffNorm);
    printf("Relative residual Size (InfNorm)   : %10.5e \n",residualNormInf);
    printf("Absolute solution error (InfNorm)  : %10.5e \n",absSolnError);
    printf("Relative solution error (InfNorm)  : %10.5e \n",relSolnError);
	printf("Initialization Time (ms)           : %10.5e \n", initializeTime);
	printf("Relaxation Sweep Time (ms)         : %10.5e \n", relaxTime);

    // Print out list of files output

    printf("\nXXXX   Files Output  XXXX \n");

    for(size_t k = 0; k < outputFileList.size(); k++)
    {
    	cout << outputFileList[k] << endl;
    }

	
	//
	// Run Fixed Time Step Solution
	//
	ClockIt clock_3;
	clock_3.start();
	iter = 0;
	diffNorm = 2 * tol;
	relaxOp.initialize(dt, alphaX, alphaY, f);
	while ((diffNorm > tol) && (iter < iterMax))
	{
		relaxOp.apply(ukf, ukfp1);

		// Check scaled relative difference between iterates

		uNorm = ukfp1.normInf();

		ukf -= ukfp1;
		diffNorm = ukf.normInf();
		diffNorm /= (dt*uNorm);

		// Update iterate

		ukf = ukfp1;
		iter++;

		if (iter % 100 == 0) cout << "Iteration : " << iter << endl;
	}// end while loop

	clock_3.stop();
	double relaxFixedTime = clock_3.getMilliSecElapsedTime();

	evaluatePoissonResidual(alphaX, alphaY, ukf, f, residualf);
	double residualfNormInf;
	residualfNormInf = residualf.normInf();
	residualfNormInf /= f.normInf();

	// Evaluate the error with respect to the continuous solution
	
	ufSolnErr = ukf;
	ufSolnErr -= uExact;

	double absSolnErrorf = ufSolnErr.normInf();
	double relSolnErrorf = absSolnError / solnNorm;

	//
	// Print out the results using standard C I/0
	//

	printf("\nXXXX   Fixed Timestep Relaxation Output XXXX \n\n");
	printf("Number of iterations               : %ld \n", iter);
	printf("Scaled Difference between iterates : %10.5e \n", diffNorm);
	printf("Relative residual Size (InfNorm)   : %10.5e \n", residualfNormInf);
	printf("Absolute solution error (InfNorm)  : %10.5e \n", absSolnErrorf);
	printf("Relative solution error (InfNorm)  : %10.5e \n", relSolnErrorf);
	printf("Time for Relaxation Sweep (ms)     : %10.5e \n", relaxFixedTime);

	return 0;
}
