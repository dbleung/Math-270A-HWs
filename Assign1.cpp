#define _USE_MATH_DEFINES // For pi constant
#include <iostream>   // For C++ input/output
#include <vector>     // For STL vectors
#include <cstdio>     // For C input/output
#include <cmath>      // For math functions
using namespace std; // Using the "standard" (std) standard template library components

/**
  Instances of class TriOperator are linear operators whose apply member function
  multiplies a tri-diagonal matrix with an input vector. 

  The tri-diagonal matrix associated with the class instance is 
  specified by it's lower diagonal, diagonal and upper diagonal entries, e.g.

   n X n tridiagonal matrix T where n = systemSize


       T =  |  diag[0]  upDiag[0]                                  |
            |loDiag[0]    diag[1] upDiag[1]                        |
            |           loDiag[1]                                  |
            |                *       *       *                     |
            |                   *         *        *               |
            |                loDiag[n-3]   diag[n-2]  upDiag[n-2]  |
            |                             loDiag[n-2]   diag[n-1]  |



  Currently, STL (Standard Template Library) vectors are used to store
  the diagonals that define the non-zero values of T

  Version: Friday Sept 29 2017 01:50:45 PM PST 
 */
class TriOperator
{
public:

    /// Null Constructor -- called when you declare an instance

    TriOperator()
    {initialize();};

    /// Copy constructor (creates duplicate of T)

    TriOperator(const TriOperator& T) // called when you declare an instance with an existing instance, & is pass be ref, so that changes to 
		// variables changes argument passed rather than value of parameter in called function
    {initialize(T);};


    /**
    This constructor initializes the internal class data with entries
    of the n X n tridiagonal matrix T associated with the entries in
    the lower diagonal (loDiag), diagonal (diag), and upper diagonal (upDiag).
    */

    TriOperator(long systemSize, vector<double>& loDiag, vector<double>& diag,
               vector<double>& upDiag)
    {
        initialize(systemSize,loDiag,diag,upDiag);
    }

    // Default destructor

    virtual ~TriOperator(){};

    /// Null initializer, resizes the internal arrays (listed on bottom of class) to zero length

    void initialize()
    {
        systemSize = 0;
        loDiag.clear();
        diag.clear();
        upDiag.clear();

    };

    ///  Copy initializer. Duplicates the entries of T

    void initialize(const TriOperator& T) // note: & is pass by reference to modify object outside function and not just a copt
		// putting a class in a variable declaration because calling object
    {
        systemSize = T.systemSize;
        loDiag     = T.loDiag;
        diag       = T.diag;
        upDiag     = T.upDiag;
    };


    /**
    This initializer initializes the internal class data with entries
    of the n X n tridiagonal matrix T associated with the entries in
    the lower diagonal (loDiag), diagonal (diag), and upper diagonal (upDiag).
    */

    void initialize(long systemSize, vector<double>& loDiag,
              vector<double>& diag, vector<double>& upDiag)
    {
        this->systemSize = systemSize;
        this->loDiag     = loDiag;
        this->diag       = diag;
        this->upDiag     = upDiag;
    }

    /**
     The apply operator evaluates

     vOut = T*vIn

     where T is the tri-diagonal operator specified by the class constructor
     or initializer.

     Currently no bounds checking is performed.

    */


    void apply(vector<double>& vIn, vector<double>& vOut) //don't you need to input diag, etc b/c they are internal
		//(because it's in"internal data below"
    {
    long i;
    
    if(systemSize == 1)
    {
    vOut[0] = diag[0]*vIn[0]; return;
    }

    vOut[0] = diag[0]*vIn[0] + upDiag[0]*vIn[1];

    for(i = 1; i < systemSize-1; i++)
    {
    vOut[i] = loDiag[i-1]*vIn[i-1] + diag[i]*vIn[i] + upDiag[i]*vIn[i+1];
    }

    i = systemSize-1;
    vOut[i] = loDiag[i-1]*vIn[i-1] + diag[i]*vIn[i];
    }


    // Internal data
    long        systemSize;
    vector<double>  loDiag;
    vector<double>    diag;
    vector<double>  upDiag;

};


class TriSolver
{
public:

	/// Null Constructor -- called when you declare an instance

	TriSolver()
	{
		initialize();
	};

	/// Copy constructor (creates duplicate of T)

	TriSolver(const TriSolver& T) // called when you declare an instance with an existing instance, & is pass be ref
	{
		initialize(T);
	};


	/**
	This constructor initializes the internal class data with entries
	of the n X n tridiagonal matrix T associated with the entries in
	the lower diagonal (loDiag), diagonal (diag), and upper diagonal (upDiag).
	*/

	TriSolver(long systemSize, vector<double>& loDiag, vector<double>& diag,
		vector<double>& upDiag)
	{
		initialize(systemSize, loDiag, diag, upDiag);
	}

	// Default destructor

	virtual ~TriSolver() {};

	/// Null initializer, resizes the internal arrays to zero length

	void initialize()
	{
		systemSize = 0;
		loDiag.clear();
		diag.clear();
		upDiag.clear();

	};

	///  Copy initializer. Duplicates the entries of T

	void initialize(const TriSolver& T) // note: & is pass by reference, changes original T made by constructor
	{
		systemSize = T.systemSize;
		loDiag = T.loDiag;
		diag = T.diag;
		upDiag = T.upDiag;

	};


	/**
	This initializer initializes the internal class data with entries
	of the n X n tridiagonal matrix T associated with the entries in
	the lower diagonal (loDiag), diagonal (diag), and upper diagonal (upDiag).
	*/

	void initialize(long systemSize, vector<double>& loDiag,
		vector<double>& diag, vector<double>& upDiag)
	{
		this->systemSize = systemSize; //"this" is a special pointer
		this->loDiag = loDiag;
		this->diag = diag;
		this->upDiag = upDiag;

	}

	//LU Factors Creates coefficients of LU factors of tridiagonal matrix given 
	//Values of the diagonals

	void LUFactors(vector<double>& c, vector<double>& gamma) {

		c[0] = diag[0];
		gamma[0] = upDiag[0] / c[0];

		for (long i = 1; i <= (systemSize-1); i++) {
			c[i] = diag[i] - loDiag[i - 1] * gamma[i - 1];
			if (i == systemSize-1) {
				break;
			}
			gamma[i] = upDiag[i] / c[i];
		}
	}

	//LUSolve Take LU factors from LUFactors function and solves tridiagonal system

	void LUSolve(vector<double>& f, vector<double>& u, vector<double>& c, vector<double>& gamma, vector<double>& v) {
		// Lv=f -> solve for v
		v[0] = f[0] / c[0];
		for (long j = 1; j <= (systemSize-1); j++) {
			v[j] = (f[j] - loDiag[j - 1] * v[j - 1]) / c[j];
		}

		// Uu=v -> solve for u
		u[systemSize-1] = v[systemSize-1];

		for (long k = systemSize - 2; k >= 0; k--) {
			u[k] = v[k] - gamma[k] * u[k + 1];

		}

	}

	/**
	The apply operator runs the tridiagonal solver for Tx=b to compute x using tridiagonal matrix T
	stored at vectors of the diagonals and vector b

	This uses two above functions: LUFactors and LUSolve

	Currently no bounds checking is performed.

	*/
	void apply(vector<double>& f, vector<double>& u, vector<double>& c, vector<double>& gamma, vector<double>& v)
	{
		
		cout << "*******RUN TRIDIAGONAL SOLVER*******" << endl;
		//LU Factorization (step 1) -> L stored as c on diag, loDiag below that, U stored as 1 on diag, gamma on upper diagonal
		LUFactors(c,gamma);
		//Take LU Factorization to solve for Lv=f and Uu=v (Step 2)
		LUSolve(f, u, c, gamma, v);
		cout << "*******TRIDIAGONAL SOLVER COMPLETE*******" << endl;
	}

	// Internal data
	long        systemSize;
	vector<double>  loDiag;
	vector<double>    diag;
	vector<double>  upDiag;

};

//
// Local utility function to set up the coefficients for the
// discrete system.
//
void setUpSystem(double h, double alpha, double beta, long M,
vector<double>& loDiag, vector<double>& diag, vector<double>& upDiag)
{
    loDiag.resize(M);
    upDiag.resize(M);
    diag.resize(M+1);

    long i;

    // First equation is the identity

    i = 0;
    diag[0]   =  1.0;
    upDiag[0] =  0.0;

    // Interior equation coefficients associated with a standard
    // second order finite difference approximation

    for(long i = 1; i < M; i++)
    {
    loDiag[i-1] =  alpha/(h*h);
    upDiag[i]   =  alpha/(h*h);
    diag[i]     = -2.0*alpha/(h*h) + beta;
    }

    // Last equation is the identity

    i          =    M;
    diag[M]    =  1.0;
    loDiag[M-1] = 0.0;
}

// Exact solutions

double Exact_a(double x)
{
    return 1.0 + x;
}

double Exact_b(double x, double alpha)
{
	return -cos(2*M_PI*x)/(4*pow(M_PI,2)*alpha);
}

double Exact_c(double x, double alpha, double beta)
{
	return cos(2 * M_PI*x) / (beta-(4 * pow(M_PI, 2)*alpha));
}

//Print Function

void print_array(vector<double>& a, long M) {
	int i;
	for (i = 0; i <= (M); i++) {
		cout << a[i] << ",";
	}
	i = M;
	cout << a[i] << endl;
}

//Main Function
int main()
{
	TriOperator tOperator;  // instantiate a forward operator ---  null constructor initializes
	TriSolver     tSolver;  // instantiate a forward operator inverse operator

	long   M = 10; // panels
	double alpha = 2.0; // coefficients
	//double beta = 0; //part A
	double beta = -1; //part c

	double xMin = 0.0; // limits of the domain
	double xMax = 1.0;
	double h = (xMax - xMin) / (double)M;


	vector<double> loDiag;   // instantiating arrays for coefficients
	vector<double> upDiag;
	vector<double> diag;

	vector<double> c(M+1, 0.0);   
	vector<double> gamma(M , 0.0);
	vector<double> v(M + 1, 0.0);

	vector<double> f(M + 1, 0.0);      // declaring solution and right hand side vectors
	vector<double> u(M + 1, 0.0);      // and initializing to zero.
	vector<double> fStar(M + 1, 0.0);

		// Assign boundary values to first and last element of right hand side

	//f[0] = 1.0; //Change this based on BC
	//f[M] = 2.0;

	//f[0] = -1 / (4 * (pow(M_PI,2)) * alpha); //Part b
	//f[M] = -1 / (4 * (pow(M_PI, 2)) * alpha);

	f[0] = 1 / (beta-(4 * (pow(M_PI,2)) * alpha)); // Part c
	f[M] = 1 / (beta - (4 * (pow(M_PI, 2)) * alpha));

	for (long i = 1; i < M; i++)
	{
		double x = xMin + i*h;
		//f[i] = 0.0; // change this based on RHS of equation
		f[i] = cos(2 * M_PI*x); // part b &c
		

	}

	// Setting up system. Outputs loDiag,diag,upDiag based on coefficients and number of panels
	setUpSystem(h, alpha, beta, M, loDiag, diag, upDiag); // finds loDiag, diag, and upDiag based of coefficient

	// Initializing the forward operator 
	tOperator.initialize(M + 1, loDiag, diag, upDiag);

	// Initializing the inverse operator
	tSolver.initialize(M + 1, loDiag, diag, upDiag);

	//print initial parameters
	cout << "********INITIAL PARAMETERS********" << endl;
	cout << "M=" << M << endl;

	cout << "diag is: " << endl;
	print_array(diag, M);

	cout << "loDiag is: " << endl;
	print_array(loDiag, M-1);

	cout << "updiag is: " << endl;
	print_array(upDiag, M-1);

	cout << "u is: " << endl;
	print_array(u, M);

	cout << "f is: " << endl;
	print_array(f, M);

	// Applying inverse
	tSolver.apply(f, u, c, gamma, v);
	
	// Applying forward operator to the result (to evaluate the residual): fstar is value of T*u from calculated u
	tOperator.apply(u, fStar); // outputs new fstar based on forward operator based on solution u just solved for

	cout << "u after solver is: " << endl;
	print_array(u, M);

	// Determining the size of the residual and the size of the error

	double aErrorMax = 0.0; //solution error infinity norm
	double fErrorMax = 0.0; //residual error infinity norm

	double aErrorTwoNorm = 0.0; // two norm of solution error
	double fErrorTwoNorm = 0.0; // two norm of the residual error
	double aErrortmp = 0.0;
	double fErrortmp = 0.0;

	double errVal;
	double x;

	for (long i = 0; i <= M; i++)
	{
		//find max error with u (solution error)
		x = xMin + i*h;
		//errVal = fabs(u[i] - Exact_a(x));
		//errVal = fabs(u[i] - Exact_b(x,alpha));
		errVal = fabs(u[i] - Exact_c(x, alpha,beta));

		if (aErrorMax < errVal) {
			aErrorMax = errVal; // inf norm
		}
		
		aErrortmp += pow(errVal, 2)*h; //stored for 2 norm


		//Now do the same thing with f (residual error)
		errVal = fabs(fStar[i] - f[i]);

		if (fErrorMax < errVal) {
			fErrorMax = errVal;
		}

		fErrortmp += pow(errVal, 2)*h;
	}

	aErrorTwoNorm = sqrt(aErrortmp);
	fErrorTwoNorm = sqrt(fErrortmp);

	
	// C style output using printf
	cout << "********ERROR RESULTS********" << endl;
	printf("(1) Inf Norm Residual error %20.15e \n", fErrorMax);
	printf("(2) Inf Norm Solution error %20.15e \n", aErrorMax);
	printf("(3) Two Norm Residual error %20.15e \n", fErrorTwoNorm);
	printf("(4) Two Norm Solution error %20.15e \n", aErrorTwoNorm);


	// C++ style output using the cout iostream
	/**
	cout.setf(ios::scientific); // Scientific notation
	cout.precision(15);         // 15 digits

	cout << endl;
	cout << "(a) Residual error " << fErrorMax << endl;
	cout << "(a) Solution error " << aErrorMax << endl;
	*/

	return 0;
}