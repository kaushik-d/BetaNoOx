#include "stdafx.h"

#if defined(WIN32) 
//currently only for win32. will need to rebuild library for 64bit if you want to use with WIN64

#include "oSparse.hpp"
#include "utility/excepts.hpp"
#include <iostream>
#include <iomanip>
#include <stdlib.h>
using namespace std;

#ifdef OLAF
//#define VSS97
#ifndef VSS97
   #define VSS96
#endif

//whit
#if defined(WIN32) || defined(WIN64)
   #define vss  VSS
   #define vss_ VSS_
   #define vss2  VSS2
   #define vss2_ VSS2_
#endif
     

#ifdef SGI
   #define VSS vss_
   #define VSS2 vss2_
#else
   #define VSS vss
   #define VSS2 vss2
#endif

 
#include <cmath>

// Added "void"  8/19/97

#ifdef VSS96
	#ifdef SGI
		extern "C"  void VSS(int *, int *, double *, double *, double *, int *, int *);
		extern "C"  void VSS2(int *, int *, double *, double *, double *, int *, int *);
	#else
		extern "C"  void __stdcall VSS(int *, int *, double *, double *, double *, int *, int *);
		extern "C"  void __stdcall VSS2(int *, int *, double *, double *, double *, int *, int *);
	#endif
#endif

#ifdef VSS97
extern "C"  void __stdcall VSS(int *, int *, double *, double *, double *, int *, int *,double *);
extern "C"  void __stdcall VSS2(int *, int *, double *, double *, double *, int *, int *,double *);
#endif

CreateErrorHandler(SymmetricOSparseMatrix);
SymmetricOSparseMatrix::SymmetricOSparseMatrix(const int numEqns) : LargeMatrix(numEqns)
{
	int i;
	Assert(col = new SortedList<int>[numberOfEquations]);
	diag = 0; 	// Allocation done in method allocate
	coefs = 0; 	// Allocation done in method allocate
	indxs = 0; 	// Allocation done in method allocate
	ptrs = 0; 	// Allocation done in method allocate
	for(i=0;i<numberOfEquations;i++)
		specifyNonZeroLocation(i,i); // Make sure there is a diagonal component
	factored=false;

	//JV101905
	K_PTRS=0;

}

SymmetricOSparseMatrix::~SymmetricOSparseMatrix()
{
	WhoMethod("~SymmetricOSparseMatrix()");
	if(coefs) delete [] coefs;
	if(diag) delete [] diag;
	if(indxs) delete [] indxs;
	if(ptrs) delete [] ptrs;
	if(col) delete [] col;

	//JV101905
	if(K_PTRS) delete [] K_PTRS;

}
SymmetricOSparseMatrix::SymmetricOSparseMatrix(const SymmetricOSparseMatrix& a) : LargeMatrix(a) 
{
	WhoClass("SymmetricOSparseMatrix");
	WhoMethod("SymmetricOSparseMatrix(const SymmetricOSparseMatrix& a)");
	int i,j;
	int num;

	neq		= a.neq;
	ncoefs	= a.ncoefs;

	K_PTRS=0;

	Assert(col = new SortedList<int>[a.numberOfEquations]);
	for(i=0;i<numberOfEquations;i++) {
		num=a.col[i].getNum();
		for(j=0;j<num;j++) {
			col[i].add(a.col[i][j]);
		}
	}
	
	Assert(coefs = new double [ncoefs]);
	for(i=0;i<ncoefs;i++) coefs[i] = a.coefs[i];

	Assert(diag = new double [numberOfEquations]);
	for(i=0;i<numberOfEquations;i++) diag[i] = a.diag[i];

	Assert(indxs = new int [ncoefs]);
	for(i=0;i<ncoefs;i++) indxs[i] = a.indxs[i];

	Assert(ptrs = new int [numberOfEquations]);
	for(i=0;i<numberOfEquations;i++) ptrs[i] = a.ptrs[i];
}


//=================================================================
void SymmetricOSparseMatrix::CopyDataFromInterface(SparseMatrixInterface *smi)
{
	//olaf stores the lower triangular matrix

	//copying information from the interface
	numberOfEquations=smi->numEquations;
	ncoefs=smi->numNonZeros-numberOfEquations;
	const int L=numberOfEquations;

	int i,j,k,numnz;
	int ctr=0;
	int ii,jj;

	Assert(coefs = new double [ncoefs]);
	Assert(diag = new double [numberOfEquations]);
	Assert(indxs = new int [ncoefs]);
	Assert(ptrs = new int [numberOfEquations]);
	for(i=0,ptrs[0]=0,k=0;i<numberOfEquations-1;i++) {
		ptrs[i+1] = ptrs[i]+smi->col[i].getNum()-1;
		for(j=0;j<smi->col[i].getNum()-1;j++,k++) {
			indxs[k]=smi->col[i][j+1];
		}
	}

	if(smi->type==LowerTriangular){
	//this is copying lower triangular information into this olaf matrix object
		for(i=0;i<L;i++) { //loop thru columns
			numnz = smi->col[i].getNum();
			ii=smi->col[i][0];
			jj=smi->row[ii].getNum();
			diag[i]		= smi->mat[ii][jj-1];
			for(j=1;j<numnz;j++) {
				ii=smi->col[i][j];
				coefs[ctr]		= (*smi)(ii,i);
				ctr++;
			}
		}
	}else if(smi->type==UpperTriangular || smi->type==Unsymmetric){
	//this is converting upper triangle to lower triangular
		cout << "SymmetricOSparseMatrix::CopyDataFromInterface not implemented for upper triangle as yet" << endl;
		/*
		for(i=0; i<L; i++) {
			rowindex[i]=ctr+1;
			numnz = smi->col[i].getNum();
			for(j=0;j<numnz;j++) {
				ii=smi->col[i][j];
				jj=smi->row[ii].find(i);
				coefs[ctr]		= smi->mat[ii][jj];
				colindex[ctr]	= smi->col[i][j]+1;
				ctr++;
			}
		}
		*/
	}
}
//=================================================================

// OSparse matrix is stored in column form
double &SymmetricOSparseMatrix::operator()(const int ii,
		const int jj)
{
	int j;
	int limit;
	if(ii>=numberOfEquations || ii < jj) {
//_operatorError:
		cerr << "numberOfEquations = " << numberOfEquations << endl;
		WhoMethod("operator()(const unsigned Int,const unsigned Int)");
		cerr << "ii,jj = " << ii << " , " << jj << endl;
		FatalError("Bad operand...");
	}
	if(ii==jj) return diag[ii];
	j=ptrs[jj];
	limit=ptrs[jj+1];
	while(indxs[j]!=ii && j<limit) j++;
	if(j==limit) {
		cerr << ii << " not in column " << jj << endl;
		cerr << "Column " << jj << "-" << ptrs[jj+1]-ptrs[jj-1] << '\n';
		for(j=ptrs[jj];j<limit;j++) cerr << ": " << indxs[j];
		cerr << endl;
		FatalError("Bad operand...");
	}
	return coefs[j];
}

void SymmetricOSparseMatrix::setaij(const int i, const int j, double val)
{
	operator()(i,j)=val;
}
double SymmetricOSparseMatrix::getaij(const int i, const int j)
{
	return operator()(i,j);
}

void SymmetricOSparseMatrix::solve(double *b,double *initialGuess) 
{
	WhoMethod("solve(double *)");

	// Check for positive definite
	CheckAndFixZeroDiagonal();

	// Create K_PTRS like vss wants
	int *K_PTRS,i;
	const int neqm1 = numberOfEquations-1;
	const int nc = ncoefs;
	Assert(K_PTRS = new int [numberOfEquations]);
	for(i=0;i<neqm1;i++) K_PTRS[i] = ptrs[i+1]-ptrs[i];
	K_PTRS[neqm1]=0;
	for(i=0;i<nc;i++) indxs[i]++;

//#define PR
#ifdef PR
	for(i=0;i<numberOfEquations;i++) cout << K_PTRS[i] << ':';
	cout << endl;
	for(i=0;i<ncoefs;i++) cout << indxs[i] << ':';
	cout << endl;
	for(i=0;i<ncoefs;i++) cout << coefs[i] << ':';
	cout << endl;
	for(i=0;i<numberOfEquations;i++) cout << diag[i] << ':';
	cout << endl;
	for(i=0;i<numberOfEquations;i++) cout << b[i] << ':';
	cout << endl;
#endif

//Whitcomb ....mod for new version of vss
double * vssAnswer;
vssAnswer = new double [numberOfEquations];


#ifdef VSS96
    VSS(&numberOfEquations,&ncoefs,diag,b,coefs,K_PTRS,indxs);
	VSS2(&numberOfEquations,&ncoefs,diag,b,coefs,K_PTRS,indxs);
#endif

#ifdef VSS97
    VSS(&numberOfEquations,&ncoefs,diag,b,coefs,K_PTRS,indxs,vssAnswer);
	VSS2(&numberOfEquations,&ncoefs,diag,b,coefs,K_PTRS,indxs,vssAnswer);
	for(i=0; i<numberOfEquations; i++)
	{b[i] = vssAnswer[i];}
#endif

  	for(i=0;i<nc;i++) indxs[i]--;
	delete [] K_PTRS;
	K_PTRS=0;


//Whitcomb
	delete [] vssAnswer;

}
#undef VSS
#undef VSS2

void SymmetricOSparseMatrix::vectorProduct(double *result, double *vector)
{
	int i,j,nt,p,ip;
//	double *mRow;
//	int *rowColNum;
	double vector_i;
	const int L=numberOfEquations;
	const int Lm1=L-1;

	for(i=0;i<L;i++) result[i] = 0;
	for(i=0;i<Lm1;i++) {
		p = ptrs[i];
		nt = ptrs[i+1]-p;
		double &result_i=result[i];
		vector_i = vector[i];
		result_i += vector_i * diag[i];
		for(j=0;j<nt;j++) {
			ip = indxs[p+j];
			result[ip] += vector_i * coefs[p+j];
			result_i += vector[ip] * coefs[p+j];
		}
	}
	result[numberOfEquations-1] += diag[numberOfEquations-1] * vector[numberOfEquations-1];
}
	
void SymmetricOSparseMatrix::zeroYourself(void)
{
	int i;
	const int L=numberOfEquations;
	const int nc=ncoefs;
	for(i=0;i<nc;i++) coefs[i]=0;
	for(i=0;i<L;i++) diag[i]=0;

	factored=false;
}

void SymmetricOSparseMatrix::printMatrix(char *message)
{
	int i,j;
	cout << message << "\nSymmetricOSparseMatrix with " << numberOfEquations <<
		" equations.\n";
	for(i=0;i<numberOfEquations;i++) cout << i << ':' << diag[i] << '\t';
	cout << '\n';
	j=1;
	for(i=0;i<ncoefs;i++) {
		if(ptrs[j]==i) {
			cout << '\n';
			j++;
			if(j>=numberOfEquations) cerr << "j >= numberOfEquations : STOP!\n";
		}
		cout << indxs[i] << ':' << coefs[i] << '\t';
	}
	cout << '\n';
}

void SymmetricOSparseMatrix::printMatrixAddressFormat(char * message, ostream* out)
{ 
	//*out << message <<endl; 

	int numNonZeros=ncoefs+numberOfEquations;
	int i,j;
	const int L=numberOfEquations;

	*out << "1" << endl; //shows that it is lower triangular.
	*out << numberOfEquations << " " << numNonZeros << endl;

#define FORMAT setw(30)<<setprecision(20)
	for(i=0;i<L;i++) 
		*out << '\t' << i << " " << i << " " << FORMAT<< diag[i];
	*out << endl;
	j=1;
	//cout << operator()(23,22) << endl;
	for(i=0;i<ncoefs;i++) {
		while (ptrs[j]==i) {
			j++;
		}
		if(j>=L) cerr << "j >= numberOfEquations : STOP!\n";
		*out << endl;
	
		*out << '\t' << indxs[i] << " " << j-1 << " " << FORMAT<< coefs[i];
	}
	*out << endl;
}
#undef FORMAT 

void SymmetricOSparseMatrix::readMatrixAddressFormat(istream* instream)
{ 
	int numNonZeros;
	*instream >> numNonZeros;

	int i;
	int ii,jj;
	double value;
	const int L=numberOfEquations;

	streampos m_posStartData;
	m_posStartData = instream->tellg(); 


	for(i=0;i<numNonZeros;i++) {
		*instream >> ii >> jj >> value;
		specifyNonZeroLocation(ii,jj);
	}
	allocate();

	instream->seekg(m_posStartData);

	for(i=0;i<numNonZeros;i++) {
		*instream >> ii >> jj >> value;
		operator()(ii,jj) = value;
	}

}

void SymmetricOSparseMatrix::specifyNonZeroLocation(const int I, const int J)
{ col[J].add(I); }

void SymmetricOSparseMatrix::allocate(void)
{
	WhoMethod("allocate(void)");
	int i,j,k;
	for(i=0,ncoefs=0;i<numberOfEquations;i++) {
//		cout << i << " - (" << col[i].getNum() << ") ";
		if(col[i].getNum()>0) {
			ncoefs+=col[i].getNum()-1;
//			col[i].print(cout);
		}
//		cout << endl;
	}
	Assert(coefs = new double [ncoefs]);
	Assert(diag = new double [numberOfEquations]);
	Assert(indxs = new int [ncoefs]);
	Assert(ptrs = new int [numberOfEquations]);
	for(i=0,ptrs[0]=0,k=0;i<numberOfEquations-1;i++) {
		ptrs[i+1] = ptrs[i]+col[i].getNum()-1;
		for(j=0;j<col[i].getNum()-1;j++,k++) {
			indxs[k]=col[i][j+1];
		}
	}

	//JV101805 keeping the non-zero locations information so that it can be used later on as well
	// used in checkNonZeroLocation(const int i, const int j)
	//so commenting out the next two lines :
	//delete [] col; 
	//col = NULL;

	zeroYourself();
}
#endif
#endif