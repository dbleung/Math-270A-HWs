#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <string>
using namespace std;

#include "TriOperator.h"
#include "TriSolver.h"

#include "TriSolver_TestCase.h"

// Detect if using a compiler that uses the C++11 compiler.
// To enable C++11 in g++ specify the option -std=c++11


int main()
{
    
    TriSolver            tSolver;
    TriOperator        tOperator;
    TriSolver_TestCase  testCase;

    // Declaring an output filestream and associating it with the file TriSolverTest.dat

    ofstream outputFile;
    outputFile.open("TriSolverTest.dat");

    // Set up test problems

    //g++ -std=c++11 TriSolver_UnitTest.cpp -o TriSolver_UnitTest.exe

    vector<double> alpha   = {2.0, 2.0,  2.0};
    vector<double>  beta   = {0.0, 0.0, -1.0};
    int problemCount       = alpha.size();

    // Set up sequence of grid approximations

    vector<long> mValues   = {10,  20,   40};

    double xMin   =  0.0;
    double xMax   =  1.0;
    double h;

    // Vectors to hold discretization coefficients

    vector<double> loDiag;
    vector<double> upDiag;
    vector<double> diag;

    // Vectors to hold right hand side and temporary for error checking

    vector<double> f;
    vector<double> u;
    vector<double> fStar;

    long   M; double x;
    double    errorMax; double residualMax;
    double     errorL2; double  residualL2;
    double        rate; double      diffVal;

    vector<double>  errorM_L2;
    vector<double> errorM_Inf;

    for(int problemIndex = 0; problemIndex < problemCount; problemIndex++)
    {
        testCase.initialize(problemIndex, alpha[problemIndex],beta[problemIndex]);
        errorM_L2.clear();

        cout       << endl << endl <<  "XXXX " << testCase.getName()  << " XXXX " << endl << endl;
        outputFile << endl << endl <<  "XXXX " << testCase.getName()  << " XXXX " << endl << endl;

        for(int mCount   = 0; mCount < mValues.size(); mCount++)
        {
        M = mValues[mCount];
        h = (xMax-xMin)/(double)M;

        f.resize(M+1,0.0);
        u.resize(M+1,0.0);
        fStar.resize(M+1,0.0);

        // Pack boundary values in f

        f[0] =  testCase.solution(xMin);
        f[M] =  testCase.solution(xMax);

		
        for(long i = 1; i < M; i++)
        {
            x    = xMin + i*h;
            f[i] = testCase.rightHandSide(x);
        }

        testCase.setUpSystem(h, alpha[problemIndex], beta[problemIndex], M, loDiag, diag, upDiag);

        tOperator.initialize(M+1,loDiag,diag,upDiag);
        tSolver.initialize(M+1,loDiag,diag,upDiag);

        // Apply the inverse operator (solver)

        tSolver.apply(f,u);

        // Apply the forward operator
		
        tOperator.apply(u,fStar);

        errorMax    = 0.0;
        residualMax = 0.0;
        errorL2     = 0.0;
        residualL2  = 0.0;

        // Compute the max and discrete L2 norm of the residual
        // and the error

        for(long i = 0; i <= M; i++)
        {
            x = xMin + i*h;
            diffVal = abs(u[i] - testCase.solution(x));
            if(errorMax < diffVal) {errorMax = diffVal;}
            errorL2 += diffVal*diffVal*h;

            diffVal = abs(fStar[i] - f[i]);
            if(residualMax < diffVal) {residualMax = diffVal;}
            residualL2 += diffVal*diffVal*h;
        }

        errorL2    = sqrt(errorL2);
        residualL2 = sqrt(residualL2);

        // Using C++ I/O for both the screeen and file output

        // Print out the results to the screen


        cout << "M :  " << M << endl;

        cout.precision(12);
        cout << "Residual error (Inf Norm) : " << residualMax << endl;
        cout << "Solution error (Inf Norm) : " << errorMax << endl << endl;

        cout << "Residual error (L2  Norm) : " << residualL2<< endl;
        cout << "Solution error (L2  Norm) : " << errorL2 << endl << endl << endl;


        // Print out the results to the file associated with outputFile fstream

        outputFile << "M :  " << M << endl;

        outputFile.precision(12);
        outputFile << "Residual error (Inf Norm) : " << residualMax << endl;
        outputFile << "Solution error (Inf Norm) : " << errorMax << endl << endl;

        outputFile << "Residual error (L2  Norm) : " << residualL2<< endl;
        outputFile << "Solution error (L2  Norm) : " << errorL2 << endl << endl << endl;


        // Cache the L2 error for rate convergence estimation

        errorM_L2.push_back(errorL2);

    }

    // Estimate the rates of convergence

    outputFile << endl;

    for(long i = 1; i < errorM_L2.size(); i++)
    {
    rate = log2(errorM_L2[i-1]/errorM_L2[i]);
    cout       << "Rate of convergence ["<< mValues[i-1] << ", " << mValues[i] << "]  : " << rate << endl;
    outputFile << "Rate of convergence ["<< mValues[i-1] << ", " << mValues[i] << "]  : " << rate << endl;
    }

    }

    // Closing file stream

    outputFile.close();
    return 0;
}
