/*
 * TriOperator.h
 *
 *  Created on: Oct. 14, 2017
 *      Author: anderson
 */

#include <vector>
#include <cstdio>
#include <cmath>

using namespace std;

#ifndef _TriOperator_
#define _TriOperator_

class TriOperator
{
public:
    
    TriOperator()
    {initialize();};

    TriOperator(const TriOperator& T)
    {initialize(T);};

    TriOperator(long systemSize, vector<double>& loDiag,
    vector<double>& diag, vector<double>& upDiag)
    {
    initialize(systemSize,loDiag,diag,upDiag);
    }

    ~TriOperator(){};

    void initialize()
    {
        systemSize = 0;
        loDiag.clear();
        diag.clear();
        upDiag.clear();

    };

    void initialize(const TriOperator& T)
    {
        systemSize = T.systemSize;
        loDiag     = T.loDiag;
        diag       = T.diag;
        upDiag     = T.upDiag;
    };

    void initialize(long systemSize, vector<double>& loDiag, vector<double>& diag, vector<double>& upDiag)
    {
        this->systemSize = systemSize;
        this->loDiag     = loDiag;
        this->diag       = diag;
        this->upDiag     = upDiag;
    }


    void apply(vector<double>& vIn, vector<double>& vOut)
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

    long        systemSize;

    vector<double>  loDiag;
    vector<double>    diag;
    vector<double>  upDiag;

};




#endif /* TRIOPERATOR_H_ */
