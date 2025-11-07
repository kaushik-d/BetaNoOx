#include "stdafx.h"

#ifdef MKL
#include "mklgbmat.hpp"
#include "utility/excepts.hpp"
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <iomanip>

//#include "mkl.h"         // Math Kernel Library 
extern "C" {
void    dgbtrf(int *m,int *n,int *kl,int *ku,double *ab,int *ldab,int *ipiv,int *info);
void    dgbtrs(char *trans,int *n,int *kl,int *ku,int *nrhs,double *ab,int *ldab,int *ipiv,double *b,int *ldb,int *info);
}

#define FORMAT setw(30) << setprecision(17)

#define STATIC_DATA static

CreateErrorHandler(MKLGeneralMatrix); 

MKLGeneralMatrix::MKLGeneralMatrix(const int numEqns) : LargeMatrix(numEqns)
{
	WhoClass("MKLGeneralMatrix");
	WhoMethod("MKLGeneralMatrix(const int)");
	int i;
	Assert(bandwidth = new int [numberOfEquations]);
	for(i=0;i<numberOfEquations;i++) bandwidth[i]=1;
	mat = 0;
	maxBandwidth = 0;
	trans = 'N';			// A x = B, then trans = 'N'
	nrhs = 1;				// # of rhs in the eqn
	NumSubDiagonals=1; //FIX THIS !!!
	NumSuperDiagonals=1; //FIX THIS !!!
}

MKLGeneralMatrix::~MKLGeneralMatrix(void)
{
	WhoMethod("~MKLGeneralMatrix()");
	if(mat!=0)	delete [] mat;	
	if(ipiv!=0) delete [] ipiv;	
	if(bandwidth !=0) delete [] bandwidth ;	
}

double& MKLGeneralMatrix::operator() (const int i, const int j)
{
	int rank;	// max rank of the Ke. gives a measure of 
				//the number of sub-diagonals and super-diagonals

	//assumption is that number of sub-diagonals=number of super-diagonals

	rank=( (compact_length-1)/3) +1;

	if( i<0 || i >numberOfEquations || j<0 || j >numberOfEquations || 
		abs(i-j) > (rank-1) )  {
			WhoMethod("operator()");
			cerr << "i,j= " << i << ", " <<j << endl;
			FatalError("Subscript out of range.\n");
	}
	
	return mat[(rank-1)*2 + i-j + (j*compact_length) ];

}

void MKLGeneralMatrix::zeroYourself(void)
{
	int i;
	WhoMethod("zeroYourself");
	if(mat==0) FatalError("Matrix not allocated yet. Can't zero.\n");
	for(i=0; i<compact_length * numberOfEquations ;i++)
			mat[i]=0;
	factored = false;
}

void MKLGeneralMatrix::solve(double *vector,double *dummy)
{
	WhoMethod("Solve");
//______________________________________________________________________
//
//           LU Factorization Using Math Kernel Library
//_______________________________________________________________________
	cout << "(6) Factorizing Global K " << endl;
    dgbtrf( &numberOfEquations, &numberOfEquations, &NumSubDiagonals, &NumSuperDiagonals, 
		    mat, &compact_length, ipiv, &info);
	cout << "MKL info code - " << info << endl;
	if(info !=0) FatalError("MKL cannot solve matrix. singular matrix!!!");
    //printVectorAsMatrix(A_compact,compact_length,totaldof);
//______________________________________________________________________
//
//           Linear Solver Using Math Kernel Library
//______________________________________________________________________	 
	cout << "(7) Solving [ K ] { q } = { F } for Linear FEA " << endl;
    dgbtrs( &trans, &numberOfEquations, &NumSubDiagonals, &NumSuperDiagonals, &nrhs,
		    mat, &compact_length, ipiv, vector, &numberOfEquations, &info);
    //Printing potentials at nodal points
//	printVector(Fg,totaldof);


}

inline double MKLGeneralMatrix::dotProduct(int length, double* v1, double *v2)
{
	STATIC_DATA int i;
	register double val;
	val=0;
	for(i=0;i<length;i++)
		val+=v1[i]*v2[i];
	return val;			
}

void MKLGeneralMatrix::printMatrix(char *s)
{
	int i,j;
	WhoMethod("printMatrix");
	cout << s<< '\n';

	for(i=0;i<compact_length;i++){
		cout << "| ";
		for(j=0;j<numberOfEquations;j++){
			cout << FORMAT << mat[j*compact_length+i] << " ";
		}
		cout << " |";
		cout << endl;
	}
	cout << endl;	
}

void MKLGeneralMatrix::specifyNonZeroLocation(const int i, const int j)
{
	// used to find out the number of sub and super diagonals.

	// For row i check to see if required bandwidth is large enough
	//		if not increase it to required size.
	if(i-j+1>bandwidth[i]) bandwidth[i]=i-j+1;
}

void MKLGeneralMatrix::allocate(void)
{
	WhoMethod("allocate(void)");
	int i,totalMemory=0;
	maxBandwidth = 0;

	for(i=0;i<numberOfEquations;i++) {
		if(bandwidth[i] > maxBandwidth) maxBandwidth = bandwidth[i];
	}
	NumSubDiagonals=NumSuperDiagonals=(maxBandwidth-1); // assuming NumSubDiagonals=NumSuperDiagonals

	compact_length = 2*NumSubDiagonals+NumSuperDiagonals+1;


	cout << "MKLGeneralMatrix::totalMemory = " << compact_length * numberOfEquations << " words." << endl;
	Assert(mat = new double [compact_length * numberOfEquations]);
	Assert(ipiv = new int[numberOfEquations]);
	
}
#ifdef DEBUG_ProfileClass
void main(void)
{ 
	int i;
	double AA[6] = {5,33,47,85,61,78};
	ProfileMatrix A,B; 

	A.allocateStorage(); 

	printf("\n-------------\n");
	A.printMatrix("A =");
	for(i=0;i<6;i++)
		printf("%4lf  ",AA[i]);
	printf("\n");
	A.solve(AA);
	printf("\n{");
	for(i=0;i<5;i++)
		printf("%lf,",AA[i]);
	printf("%lf}\n",AA[5]);
	A.printMatrix("A =");
	printf("End of Test Matrix A::Solve()\n\n");

	B.allocateStorage();
	B.printMatrix("B =");
	B(0,0)=2;
	B.printMatrix("B =");
	B(1,0)=4;
	B.printMatrix("B =");
	B(1,1)=5;
	B.printMatrix("B =");
	B(5,2)=99;
	B.printMatrix("B =");
}
#endif
#undef FORMAT
#endif