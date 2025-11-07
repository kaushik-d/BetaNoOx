#include "stdafx.h"

#ifdef MKL
#include <cmath>
#include <iostream>
#include <iomanip>
#include "MKLPardisoUnsymm.hpp"
#include "utility/excepts.hpp"
#include "utility/progressIndicator.hpp"
#include "sparseMatrixInterface.hpp"

extern "C" {int omp_get_max_threads();}
/* PARDISO prototype. */
#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif
extern "C" {
int PARDISO
	(void *, int *, int *, int *, int *, int *,
	double *, int *, int *, int *, int *, int *,
	int *, double *, double *, int *);
}


CreateErrorHandler(MKLPardisoUnsymmMatrix);
//=============================================================
MKLPardisoUnsymmMatrix::MKLPardisoUnsymmMatrix(const int numEqns) :
		MKLPardisoSymmMatrix(numEqns)
{
	WhoClass("MKLPardisoUnsymmMatrix");
	WhoMethod("MKLPardisoUnsymmMatrix(const int)");
}
//=============================================================
MKLPardisoUnsymmMatrix::~MKLPardisoUnsymmMatrix()
{
	DeleteWhoMethod("~MKLPardisoUnsymmMatrix()");
//	cerr << "Deleting MKLPardisoUnsymmMatrix matrix" << endl;
}
//=============================================================
void MKLPardisoUnsymmMatrix::binaryWrite(FILE *outfile)
{
	WhoMethod("binaryWrite(FILE *)");
	const int L=numberOfEquations;
	/*
	fwrite(&numberOfEquations,sizeof(int),1,outfile);
	for(i=0;i<L;i++) row[i].fwrite(outfile);
	for(i=0;i<L;i++) fwrite(mat[i],sizeof(double),row[i].getNum(),outfile);
	*/
}
//=============================================================
void MKLPardisoUnsymmMatrix::binaryRead(FILE *infile)
{
	WhoMethod("binaryRead(FILE *)");
	const int L=numberOfEquations;
	/*
	if(L>0) {
		for(i=0;i<L;i++)
			delete [] mat[i];
		delete [] mat;
		delete [] row;
	}
	fread(&numberOfEquations,sizeof(int),1,infile);
	Assert(row = new SortedList<int>[numberOfEquations]);
	for(i=0;i<L;i++) row[i].fread(infile);
	Assert(mat = new double *[numberOfEquations]);
	for(i=0;i<L;i++) {
		Assert(mat[i] = new double[row[i].getNum()]);
		fread(mat[i],sizeof(double),row[i].getNum(),infile);
	}
	*/
}
//=============================================================
double& MKLPardisoUnsymmMatrix::operator()(const int ii,
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

	j=rowindex[jj]-1;
	limit=rowindex[jj+1]-1;
	while(colindex[j]-1!=ii && j<limit) j++;
	if(j==limit) {
		cerr << ii << " not in row " << jj << endl;
		cerr << "Row " << jj << "-" << rowindex[ii+1]-rowindex[ii] << '\n';
		for(j=rowindex[ii];j<limit;j++) cerr << ": " << colindex[j];
		cerr << endl;
		FatalError("Bad operand...");
	}
	return coefs[j];
}
//=============================================================
void MKLPardisoUnsymmMatrix::addtoaij(const int i, const int j, double val)
{
	if(i>=numberOfEquations || j>=numberOfEquations ) {
		cerr << "numberOfEquations = " << numberOfEquations << endl;
		WhoMethod("setaij(const int i, const int j, double val)");
		cerr << "ii,jj = " << i << " , " << j << endl;
		FatalError("Bad operand...");
	}
	if(i>j)
		operator()(j,i) += val;
	else
		operator()(i,j) += val;
}
//=============================================================
void MKLPardisoUnsymmMatrix::solve(double *b,double *initialGuess) 
{
	WhoMethod("solve(double *)");
	register int i=0;
//	int j,k;
	const int L=numberOfEquations;
	bool convergedSolution = false,zeroOnDiagonal = false;
	DEBUG3(printMatrix("DEBUG is ON"));
	// Check for positive definite
	for(i=0;i<L;i++){  // Quick check not a true check.
		if(operator()(i,i) == 0) {
			zeroOnDiagonal = true; 
			break;
		}
	}

	if(zeroOnDiagonal) 	{	//+DG_Feb2004
		if(zeroDiagonalReplacementValue==0) { //+DG_Feb2004
			cerr<<"The values are zero on diagonal - (" << i <<"," << i << ")";
			cerr <<"\nReplacement for Zero Diagonal not set.\nUse ReplaceZeroDiagonal_Option\n\n";
			exit(1);
		}
		Warning("Zero on diagonal.  Will attempt to fix.");
		for(i=0;i<L;i++)
			if(operator()(i,i)==0) {
				cerr << "Changing [" << i << ',' << i << "] to "<< zeroDiagonalReplacementValue<<".\n";
				operator()(i,i) = zeroDiagonalReplacementValue;
			}
	}
	//-DG_Feb2004

	int mtype = 11; /* real and unsymmetric matrix */
	int nrhs = 1; /* Number of right hand sides. */
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	int iparm[64];
	int maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	double ddum; /* Double dummy */
	int idum; /* Integer dummy. */
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}
	iparm[0] = 1; /* No solver default */
	iparm[1] = 2; /* Fill-in reordering from METIS */
	/* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = omp_get_max_threads();
	iparm[3] = 0; /* No iterative-direct algorithm */
	iparm[4] = 0; /* No user fill-in reducing permutation */
	iparm[5] = 0; /* Write solution into x */
	iparm[6] = 0; /* Not in use */
	iparm[7] = 2; /* Max numbers of iterative refinement steps */
	iparm[8] = 0; /* Not in use */
	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0; /* Not in use */
	iparm[12] = 0; /* Not in use */
	iparm[13] = 0; /* Output: Number of perturbed pivots */
	iparm[14] = 0; /* Not in use */
	iparm[15] = 0; /* Not in use */
	iparm[16] = 0; /* Not in use */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0; /* Output: Numbers of CG Iterations */
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 1; /* Print statistical information in file */
	error = 0; /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		pt[i] = 0;
	}
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&numberOfEquations, coefs, rowindex, colindex, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&numberOfEquations, coefs, rowindex, colindex, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */

	double *x;
	Assert(x = new double [numberOfEquations]);

	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&numberOfEquations, coefs, rowindex, colindex, &idum, &nrhs,
		iparm, &msglvl, b, x, &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	for(i=0;i<L;i++) b[i]=x[i];
	delete [] x;
	
	phase = -1; /* Release internal memory. */
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&numberOfEquations, &ddum, rowindex, colindex, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);

}
//=============================================================
void MKLPardisoUnsymmMatrix::printMatrix(char *message)
{
	const int L=numberOfEquations;
	cout << message << "\nMKLPardisoUnsymmMatrix with " << numberOfEquations <<
		" equations.\n";
	/*
	for(i=0;i<L;i++) {
		for(j=0,limit=row[i].getNum();j<limit;j++)
			cout << '\t' << row[i][j] << ':' << mat[i][j];
		cout << '\n';
	}
	*/
}
//=============================================================
void MKLPardisoUnsymmMatrix::printMatrixAddressFormat(char * message, ostream* out)
{ 
	int numNonZeros=0;
	int i;
	const int L=numberOfEquations;

	*out << "3" << endl;
	*out << numberOfEquations << " " << ncoefs << endl;

#define FORMAT setw(30)<<setprecision(20)

	int ii=0, jj;

	for(i=0,ii=0; i<ncoefs; i++) {
		if(i==rowindex[ii+1]-1){
			*out << '\n';
			ii++;
		}
		jj=colindex[i]-1;
		*out << '\t' << ii << " " << jj << " " << FORMAT<< coefs[i];
	}
	*out << endl;
}
#undef FORMAT 

//=============================================================
void MKLPardisoUnsymmMatrix::readMatrixAddressFormat(istream* instream)
{ 
}
//=============================================================
void MKLPardisoUnsymmMatrix::specifyNonZeroLocation(const int I, const int J)
{ 
	//check this !!!!!!!!!!!!!!!!
	//row[I].add(J); 
	row[J].add(I); // pardiso stores the upper triangle
	// the element matrices store the lower triangle 
	//so the indices need to be switched.
}
//=============================================================
void MKLPardisoUnsymmMatrix::allocate(void)
{
	WhoMethod("allocate(void)");
	const int L=numberOfEquations;
	int i,j,k;
	for(i=0,ncoefs=0;i<L;i++) {
//		cout << i << " - (" << row[i].getNum() << ") ";
		if(row[i].getNum()>0) {
			ncoefs+=row[i].getNum();
//			row[i].print(cout);
		}
//		cout << endl;
	}
	Assert(coefs = new double [ncoefs]);
	Assert(colindex = new int [ncoefs]);
	Assert(rowindex = new int [L+1]);
	for(i=0,rowindex[0]=1,k=0; i<L; i++) {
		rowindex[i+1] = rowindex[i] + row[i].getNum();
		for(j=0;j<row[i].getNum();j++,k++) {
			colindex[k]=row[i][j]+1;
		}
	}
	//rowindex[L] = ncoefs+1;
	delete [] row; 
	row = NULL;
	zeroYourself();
}
//=============================================================
#ifdef DEBUG_MKLPardisoUnsymmMatrix
ofstream Logfile;
void main(void)
{
	int i;
	ifstream is;
	ofstream os;
	is.open("GlobalStiffness.txt");
	int NumEquations;
	is >> NumEquations;
	MKLPardisoUnsymmMatrix A(NumEquations);
	cout << "Start reading GlobalStiffness" << endl;
	A.readMatrixAddressFormat(&is);
	is.close();
	cout << "Finished reading GlobalStiffness" << endl;

	is.open("LoadVector.txt");
	double *Load;
	Load=new double[NumEquations];
	for(i=0;i<NumEquations;i++){
		is >> Load[i];
	}
	is.close();
	cout << "Finished reading Load Vector" << endl;

	A.solve(Load);
	cout << "Finished Solving" << endl;

	os.open("Solution.test.txt");
	for(i=0;i<NumEquations;i++) os << Load[i] << endl;
	os.close();
	cout << "Finished writing solution" << endl;

	delete [] Load;
}
#endif
#endif
