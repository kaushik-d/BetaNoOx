#include "stdafx.h"

#include "largeMatrix.hpp"
#include "utility/excepts.hpp"
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <cmath>

//=============================================================
CreateErrorHandler(LargeMatrix);
//=============================================================
LargeMatrix::LargeMatrix(const int numEqns)
{
	numberOfEquations = numEqns; 
	zeroDiagonalReplacementValue=0;
	verbose=Off;

	UsePreviousFillInOrdering=false;
	UsePreviousFactorizedMatrix=false;
	OrderingCompleted=false;
	FactorizeCompleted=false;
	num_threads_for_solver=1;
}
//=============================================================
LargeMatrix::LargeMatrix(const LargeMatrix& arg)
{
	numberOfEquations = arg.numberOfEquations; 
	zeroDiagonalReplacementValue=arg.zeroDiagonalReplacementValue;
	verbose=arg.verbose;

	UsePreviousFillInOrdering=arg.UsePreviousFillInOrdering;
	UsePreviousFactorizedMatrix=arg.UsePreviousFactorizedMatrix;
	OrderingCompleted=arg.OrderingCompleted;
	FactorizeCompleted=arg.FactorizeCompleted;
	num_threads_for_solver=arg.num_threads_for_solver;
}
//=============================================================
LargeMatrix::~LargeMatrix(void)
{
}
//=============================================================
void LargeMatrix::CheckAndFixZeroDiagonal()
{
	bool zeroOnDiagonal=false;
	const int L=numberOfEquations;
	int diag_index=-1;
	int i;
	// Check for positive definite
	for(i=0;i<L;i++){  // Quick check not a true check.
		if(operator()(i,i) == 0) {zeroOnDiagonal = true; diag_index = i; break;}
		if(operator()(i,i) < 0) {
			cerr<<"The diagonal term(s) are less than zero - here is one of them : (" << i <<"," << i << ") = " << operator()(i,i);
			exit(1);
		}
	}

	if(zeroOnDiagonal) 
	{
		if(zeroDiagonalReplacementValue==0)
		{
			cerr <<"The diagonal term(s) are zero - here is one of them : (" << i <<"," << i << ")\n";
			cerr <<"If you want to replace the zero diagonals with another value, use the command 'ReplaceZeroDiagonal_Option' to set the replacement value";
			exit(1);
		}
		if(verbose >= Basic){
			cerr << "Atleast one zero on diagonal. Eq. number: " << diag_index << endl;
			Warning("Zero(s) on diagonal.  Will attempt to fix.");
		}

		for(i=0;i<L;i++)
			if(operator()(i,i)==0) 
			{
//				if(verbose >= Basic)	
//					cerr << "Changing [" << i << ',' << i << "] to "<< zeroDiagonalReplacementValue<<".\n";
				operator()(i,i) = zeroDiagonalReplacementValue;
			}
	}
}
