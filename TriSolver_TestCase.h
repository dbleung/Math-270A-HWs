/*
 * TriSolver_TestCase.h
 *
 *  Created on: Oct. 14, 2017
 *      Author: anderson
 */

#ifndef _TriSolver_TestCase_
#define _TriSolver_TestCase_


#include <vector>
#include <cmath>
using namespace std;

/**
Class  TriSolver_TestCase provides data for test problems used to verify
the functionality of the TriSolver class.

Currently  three test problems are supported.

problemA: alpha [D+D-       ] u = 0              u(0) = 1                         u(1) = 2
problemB: alpha [D+D-       ] u = cos(2 pi x)    u(0) = -1/(4 alpha pi^2)         u(1) = -1/(4 alpha pi^2)
problemC: alpha [D+D- + beta] u = cos(2 pi x)    u(0) =  1/(beta - 4 alpha pi^2)  u(1) = 1/(beta - 4 alpha pi^2)

    
D+D- is the standard second order three-point discretization of the Laplacian.
*/

class TriSolver_TestCase
{
public:

     enum {problemA = 0, problemB = 1, problemC = 2}; // Specifically specifying value (equivalent to default)
                                                      // to emphasize the problem can also be referenced as 0, 1 and 2
     TriSolver_TestCase()
     {initialize();}

     /**
      Initializes the test case to the specified problem type (int index)
      with the coefficients alpha and beta
      */

     TriSolver_TestCase(int problemType, double alpha, double beta)
     {initialize(alpha,beta,problemType);}

     void initialize()
     {
            this->alpha       = 0.0;
            this->beta        = 0.0;
            this->problemType = -1;
            this->pi    = 3.14159265358979323846264338327950288419;
     }

     void initialize(int problemType, double alpha, double beta)
     {
            this->alpha       = alpha;
            this->beta        = beta;
            this->problemType = problemType;
            this->pi          = 3.14159265358979323846264338327950288419;
     }

     /// Evaluates the solution of the test problem

     double solution(double x)
     {
        switch(problemType)
        {
        case  problemA : return 1.0 + x;
        case  problemB : return -cos(2.0*pi*x)/(4.0*pi*pi*alpha);
        case  problemC : return  cos(2.0*pi*x)/(beta - 4.0*pi*pi*alpha);
        default:
            cerr << "TriSolver_UnitTest class error : Improper problem specification " << endl;
            cerr << "XXX  Execution halted XXXX ";
            exit(1);
        }
     }

     /// Evaluates the right hand side of the test problem

     double rightHandSide(double x)
     {
        switch(problemType)
        {
        case  problemA : return 0.0;
        case  problemB : return cos(2.0*pi*x);
        case  problemC : return cos(2.0*pi*x);
        default:
            cerr << "TriSolver_UnitTest class error : Improper problem specification " << endl;
            cerr << "XXX  Execution halted XXXX ";
            exit(1);
        }
     }

     /// Returns the name of the test problem

     string getName()
     {
      switch(problemType)
      {
        case  problemA : return "Problem A";
        case  problemB : return "Problem B";
        case  problemC : return "Problem C";
     }
     }


     /** This routine sets up the tri-diagonal matrix corresponding to the discrete operator

               = u_0                        // identity at i = 0
        L u    = (alpha*D+D- + beta) u_i    // second order discretization for interior points
               = u_M                        // identity at i = M

       The associated grid is one of M panels and M+1 grid points.

       M = the number of panels

     */

     void setUpSystem(double h, double alpha, double beta, long M,
     vector<double>& loDiag, vector<double>& diag, vector<double>& upDiag)
     {
         loDiag.resize(M);
         upDiag.resize(M);
         diag.resize(M+1);


         // identity at left edge

         diag[0]   =  1.0;
         upDiag[0] =  0.0;


         for(long i = 1; i < M; i++)
         {
         loDiag[i-1] =  alpha/(h*h);
         upDiag[i]   =  alpha/(h*h);
         diag[i]     = -2.0*alpha/(h*h) + beta;
         }

         // identity at right edge

         diag[M]    =  1.0;
         loDiag[M-1] = 0.0;
     }

     double    alpha;
     double     beta;
     int problemType;
     double       pi;
};


#endif
