#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

//
// Math 270A : Assignment 3
//
// Test program for varaible time step relaxation to a solution of
//
// alpha * d^2 U/dx^2  = f  with Dirichlet boundary conditions
//
// for x in [xMin, xMax]
//
//
// Version : Tues. Oct. 21, 2017 
//
#include "GridFun1D.h"  // 1D grid function class
#include "RelaxOp1D.h"  // 1D Crank-Nicholson relaxation operator class
#include "Residual.h"

class TestProblem1D
{
public:

    TestProblem1D(double alpha)
    {
        this->alpha = alpha;
        this->pi    = 3.1415926535897932;
    }

    double operator()(double x)
    {
        return -cos(2.0*pi*x)/(4.0*pi*pi*alpha);
    }

    double pi;
    double alpha;

};

//Function to print arrays to the screen
void print_array(vector<double>& a, long M) {
	int i;
	for (i = 0; i < (M); i++) {
		cout << a[i] << ",";
	}
	i = M;
	cout << a[i] << endl;
}

// Function to find the variable timestep dt based on current iteration value of j and initial grid size ho=[b - a]/M
double find_dt(int j, double ho, double alpha, double pi) {

	return 4*pow(ho*pow(2, j), 2) / (alpha*pi*pi);
}

int main()
{
	//
	// Set up test problem
	//

	double alpha = 2.0;  // Laplace operator prefactor
	double dt = 0.1;  // Relaxation timestep
	//double dtTemp; // For variable time stepping
	long   M = 10;   // Panel count
	int jstar = int(log2(M)) + 1; // step count for variable timestep (nearest integer of log2(M))
 	double tol = 1.0e-6;  // Stopping tolerance


	double xMin = 0.0;
	double xMax = 1.0;
	double h = (xMax - xMin) / (double)M;
	double ho = (xMax - xMin) / (double)M;

	//
	// Instantiate the test problem solution
	//

	TestProblem1D      testSoln(alpha);


	// Instantiate grid functions for a discretization 
	// of [xMin,xMax] with M panels.
	//

	GridFun1D f(M, xMin, xMax); // RHS vector for target system of odes (for Trisolver calculation of du/dt = Bu -f
	GridFun1D fExactResidual(M, xMin, xMax); // RHS vector for target discrete equation (used for residual calculation of Au=f)
	GridFun1D uk(M, xMin, xMax); //u_k
	GridFun1D ukp1(M, xMin, xMax);//u_k+1
	GridFun1D vTmp(M, xMin, xMax); // for variable time stepping
	GridFun1D vTmpp1(M, xMin, xMax); // for variable time stepping, plus 1 term
	GridFun1D vBarTmp(M, xMin, xMax);// for variable time stepping
	GridFun1D vBarTmpm1(M, xMin, xMax);// for variable time stepping, minus 1 term
	GridFun1D dtTemp(jstar-1, xMin, xMax);//for variable time stepping, keeps values of dt for forward and backward sweep
	GridFun1D fApprox(M, xMin, xMax); // For residual, approximate value of RHS f
	//
	// Initialize right hand side
	//

	double x;
	double pi = 3.1415926535897932;
	for (long i = 0; i <= M; i++)
	{
		x = xMin + i*h;
		f.values[i] = cos(2.0*pi*x);
		fExactResidual.values[i] = cos(2.0*pi*x);
	}

	// Set edge values of f to zero

	f.values[0] = 0.0;
	f.values[M] = 0.0;
	cout << "f is:" << endl;
	print_array(f.values, M);

	//
	// Instantiate and initialize relaxation operator
	//

	RelaxOp1D relaxOp;
	
	//cout << "INITIALIZE" << endl;
	//relaxOp.initialize(dt,alpha,f); 

	//
	// Set up initial iterate: zero in the interior, boundary
	// values at the edges.
	//

	uk.setToValue(0.0); //uses function in GridFun1D
	uk.values[0] = testSoln(xMin); // BC, uses overloaded operator
	uk.values[M] = testSoln(xMax); // BC, uses overloaded operator
	cout << "uk is:" << endl;
	print_array(uk.values, M);
	
	//
	//Get dt vector for multi scale forward and backward sweep
	//
	
	for (int j = 0; j <= jstar-1; j++) {
		dtTemp.values[j] = find_dt(j, ho, alpha, pi);
	}

	cout << "dtTemp" << endl;
	print_array(dtTemp.values, jstar - 1);

    //
    // Initialize relaxation loop
    //


    //double diffNorm = 2*tol;
    long   iterMax  = 4000;
    long   iter     = 0;
	int relaxsteps = 0;
	vector<double> fapprox(M + 1, 0.0);
	double twoNorm = 2*tol;
	double MaxNorm = 2 * tol;

	//Initialize residual function
	Residual residualOp;

    while((MaxNorm > tol)&&(iter < iterMax))
    {

	//initialize forward sweep
		vTmp = uk;

	// Go through forward sweep
	
		for (int j = 0; j <= (jstar-1); j++) {
			relaxOp.initialize(dtTemp.values[j], alpha, f);
			relaxOp.apply(vTmp, vTmpp1); 
			vTmp = vTmpp1;
			relaxsteps++;
		}
	

	// Go through backward sweep
		vBarTmp = vTmp;
		for (int j = (jstar - 1); j >= 0; j--) {
			relaxOp.initialize(dtTemp.values[j], alpha, f);
			relaxOp.apply(vBarTmp, vBarTmpm1);
			vBarTmp = vBarTmpm1;
			relaxsteps++;
		}
		

	ukp1 = vBarTmp;
	cout << "**********ITERATION " << iter +1 << "***********" << endl;
    //relaxOp.apply(uk,ukp1); //Solve for next iteration given current iteration


    // Check relative difference between iterates
	residualOp.apply(uk.values,fApprox.values,alpha,M,h); //Gets the approximate f value in Au = f

	//Calculate 2 norm of difference bewteen f and fapprox
	double sumDiff = 0;
	double absDiff = 0;
	MaxNorm = 0;
		for(int i = 0; i <= M; i++) {
			sumDiff += pow((fApprox.values[i] - fExactResidual.values[i]),2);
			absDiff = fabs(fApprox.values[i] - fExactResidual.values[i]);
			//cout << "abdDiff:" << absDiff << endl;
			if (absDiff > MaxNorm) {
				MaxNorm = absDiff;
			}
		}

	twoNorm = sqrt(sumDiff);
	cout << "Two Norm of Residual:" << endl;
	cout << twoNorm << endl;
	cout << "Inf Norm of Residual:" << endl;
	cout << MaxNorm << endl;
			
//    uk       -= ukp1;

 //   diffNorm  = uk.normInf();   // Norm of time derivative
 //   diffNorm /= dt;             

    // Update iterate

    uk  = ukp1;
    iter++;
    }

    //
    // Evaluate the error
    //

	cout << "Variable Time Step Final Soulution Is:" << endl;
	print_array(uk.values, M);


	//
	// Calculate Fixed Time Step Results
	//

	cout << "************FIXED TIME STEP**************" << endl;
	double diffNorm = 2 * tol;
	long iterF = 0;
	GridFun1D ukf(M, xMin, xMax); //u_k fixed step
	GridFun1D ukfp1(M, xMin, xMax);//u_k+1 fixed time step
//
// Set up initial iterate: zero in the interior, boundary
// values at the edges.
//

	ukf.setToValue(0.0); //uses function in GridFun1D
	ukf.values[0] = testSoln(xMin); // BC, uses overloaded operator
	ukf.values[M] = testSoln(xMax); // BC, uses overloaded operator

	relaxOp.initialize(dt, alpha, f);

	while ((diffNorm > tol) && (iterF < iterMax))
	{
		//cout << "**********ITERATION " << iter + 1 << "***********" << endl;
		relaxOp.apply(ukf, ukfp1); //Solve for next iteration given current iteration

		 // Check relative difference between iterates

		ukf -= ukfp1;
		diffNorm = ukf.normInf();   // Norm of difference
		diffNorm /= dt;

		// Update iterate

		ukf = ukfp1;
		iterF++;
	}

	cout << "Fixed Time Step Final Soulution Is:" << endl;
	print_array(ukf.values, M);

	double errorMax = 0.0;
	double errVal = 0.0;
	for (long i = 0; i <= M; i++)
	{
		x = xMin + i*h;
		errVal = abs(uk.values[i] - testSoln(x)); //, uses overloaded operator
		if (errorMax < errVal) { errorMax = errVal; }
	}

    //
    // Print out the results (using standard C I/O because I find it
    // easier to create nice output).
    //

    printf("XXXX   Output XXXX \n\n");
    printf("Panel Count                                               : %ld  \n",M);
    //printf("Timestep    : %10.5e \n",dt);

    printf("Number of iterations for Variable Time Step               : %ld \n",iter);
	printf("Number of relaxation steps for Variable Time Step         : %1d \n", relaxsteps);
	printf("Number of iterations for Fixed Time Step                  : %ld \n", iterF);
	printf("Max Difference between iterates, Fixed Time Step          : %10.5e \n", diffNorm);
    //printf("Iterative solution error    : %10.5e \n",errorMax);
	printf("Residual using Inf Norm, Variable Time Step               : %10.5e \n", twoNorm);
	printf("Iterative solution error using Inf Norm, Fixed Time Step  : %10.5e \n", errorMax);


    return 0;
}

