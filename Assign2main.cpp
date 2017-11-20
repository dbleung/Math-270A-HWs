#include <cstdio>
#include <cmath>
#include <cstdlib>
using namespace std;

//
// Math 270A : Assignment 2
//
// Test program for fixed step relaxation to a solution of
//
// alpha * d^2 U/dx^2  = f  with Dirichlet boundary conditions
//
// for x in [xMin, xMax]
//
//
// Version : Tues. Oct. 10, 2017 06:48:20 PM PST 
//
#include "GridFun1D.h"  // 1D grid function class
#include "RelaxOp1D.h"  // 1D Crank-Nicholson relaxation operator class

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

int main()
{
    //
    // Set up test problem
    //
	cout << "test" << endl;

    double alpha  =     2.0;  // Laplace operator prefactor
    double dt     =     0.1;  // Relaxation timestep
    long   M      =     10;   // Panel count
    double tol    =  1.0e-6;  // Stopping tolerance


    double xMin   =     0.0;
    double xMax   =     1.0;
    double h      = (xMax-xMin)/(double)M;
	vector<double> values;
    //
    // Instantiate the test problem solution
    //

    TestProblem1D      testSoln(alpha);


    // Instantiate grid functions for a discretization 
    // of [xMin,xMax] with M panels.
    //

    GridFun1D f(M,xMin,xMax,values);          
    GridFun1D uk(M,xMin,xMax,values);
    GridFun1D ukp1(M,xMin,xMax,values);
    GridFun1D uTmp(M,xMin,xMax,values); //not used?

    //
    // Initialize right hand side
    //

    double x;
    double pi =  3.1415926535897932;
    for(long i = 0; i <= M; i++)
    {
        x          =   xMin + i*h;
       f.values[i] = cos(2.0*pi*x);
    }

    // Set edge values of f to zero

    f.values[0] = 0.0;
    f.values[M] = 0.0;

    //
    // Instantiate and initialize relaxation operator
    //

    RelaxOp1D relaxOp;
    relaxOp.initialize(dt,alpha,f); //what is this function supposed to do?

    //
    // Set up initial iterate: zero in the interior, boundary
    // values at the edges.
    //

    uk.setToValue(0.0); //uses function in GridFun1D
    uk.values[0] = testSoln(xMin); // BC, uses overloaded operator
    uk.values[M] = testSoln(xMax); // BC, uses overloaded operator

    //
    // Initialize relaxation loop
    //


    double diffNorm = 2*tol;
    long   iterMax  = 4000;
    long   iter     = 0;

    while((diffNorm > tol)&&(iter < iterMax))
    {
    relaxOp.apply(uk,ukp1); //Solve for next iteration given current iteration


    // Check relative difference between iterates

    uk       -= ukp1;
    diffNorm  = uk.normInf();   // Norm of time derivative
    diffNorm /= dt;             

    // Update iterate

    uk  = ukp1;
    iter++;
    }

    //
    // Evaluate the error
    //

    double errorMax = 0.0;
    double errVal   = 0.0;
    for(long i = 0; i <= M; i++)
    {
            x = xMin + i*h;
            errVal = abs(uk.values[i] -testSoln(x)); //, uses overloaded operator
            if(errorMax < errVal) {errorMax = errVal;}
    }

    //
    // Print out the results (using standard C I/O because I find it
    // easier to create nice output).
    //

    printf("XXXX   Fixed Timstep Relaxation Output XXXX \n\n");
    printf("Panel Count : %ld  \n",M);
    printf("Timestep    : %10.5e \n",dt);
    printf("Difference between iterates : %10.5e \n",diffNorm);
    printf("Number of iterations        : %ld \n",iter);
    printf("Iterative solution error    : %10.5e \n",errorMax);

    return 0;
}
