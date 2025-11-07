#include "stdafx.h"

#ifdef MKL
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
//#include "mkl.h"         // Math Kernel Library 
#include "MKLPardiso.hpp"
#include "utility/utility_defines.h"
#include "utility/excepts.hpp"
#include "sparseMatrixInterface.hpp"
//#include <omp.h>

extern "C" {void MKL_Set_Num_Threads(int nth);}
extern "C" {void  ompc_set_num_threads(int);}
extern "C" {void omp_set_num_threads(int);}
extern "C" {int omp_get_max_threads();}
extern "C" {int omp_get_num_threads();}
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

void writeCSRMatrixInBinary(char* filename, int* &ia, int* &ja, double* &a, int numberOfEquations, int ncoefs);
void readCSRMatrixInBinary(char* filename, int* &ia, int* &ja, double* &a, int &numberOfEquations, int &ncoefs);
void writeCSRMatrixInASCII(char* filename, int* &ia, int* &ja, double* &a, int numberOfEquations, int ncoefs);

CreateErrorHandler(MKLPardisoSymmMatrix);
//=============================================================
MKLPardisoSymmMatrix::MKLPardisoSymmMatrix(const int numEqns) :
		LargeMatrix(numEqns)
{
	WhoClass("MKLPardisoSymmMatrix");
	WhoMethod("MKLPardisoSymmMatrix(const int)");
	Assert(row = new SortedList<int>[numberOfEquations]);
	coefs = 0; 	// Allocation done in method allocate
	rowindex = 0; 	// Allocation done in method allocate
	colindex = 0; 	// Allocation done in method allocate
	pardiso_x = 0;

	int i;
	for(i=0;i<numberOfEquations;i++)
		specifyNonZeroLocation(i,i); // Make sure there is a diagonal component

	DoInitializePardisoParameters=true;
	DoReleasePardisoMemory=true;
	if(pardiso_x==0) //this array is required by pardiso when solving (back substitution)
		Assert(pardiso_x = new double [numberOfEquations]);
}
//=============================================================
MKLPardisoSymmMatrix::~MKLPardisoSymmMatrix()
{
	DeleteWhoMethod("~MKLPardisoSymmMatrix()");
//	cerr << "Deleting MKLPardisoSymmMatrix matrix" << endl;

	//if(DoReleasePardisoMemory) ReleasePardisoMemory();
	if(pardiso_x) {delete [] pardiso_x; pardiso_x=0;}

	int L=numberOfEquations;
	if(row) delete [] row;
	if(coefs) delete [] coefs;
	if(colindex) delete [] colindex;
	if(rowindex) delete [] rowindex;
}
//=============================================================
void MKLPardisoSymmMatrix::ReleaseMemory()
{
	//pardiso call
	phase = -1; // Release internal memory. 
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&numberOfEquations, &ddum, rowindex, colindex, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);

	FactorizeCompleted = false;

	if(!UsePreviousFillInOrdering){
		OrderingCompleted = false; }

}
//=============================================================
void MKLPardisoSymmMatrix::binaryWrite(FILE *outfile)
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
void MKLPardisoSymmMatrix::binaryRead(FILE *infile)
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
double& MKLPardisoSymmMatrix::operator()(const int ii,
		const int jj)
{
	//=============================
	//testing!!!!
//	return coefs[0];
	//=============================
// the input parameters assume that the matrix is storing lower diagonal info
// but the matrix is actually storing upper diagonal info, so the parameter check below
// is to check if the parameters follow a lower diagonal matrix
// but after that, parameters are in effect switched for the algorithm to locate the matrix entry
	int j;
	int limit;
	if(ii>=numberOfEquations || ii < jj) {
//_operatorError:
		cerr << "numberOfEquations = " << numberOfEquations << endl;
		WhoMethod("operator()(const unsigned Int,const unsigned Int)");
		cerr << "ii,jj = " << ii << " , " << jj << endl;
		FatalError("Bad operand...");
	}

	//actually locating Kg(jj,ii): (because Kg(ii,jj) assumes lower diagonal storage)
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
void MKLPardisoSymmMatrix::setaij(const int i, const int j, double val)
{
	if(i>=numberOfEquations || j>=numberOfEquations ) {
		cerr << "numberOfEquations = " << numberOfEquations << endl;
		WhoMethod("setaij(const int i, const int j, double val)");
		cerr << "ii,jj = " << i << " , " << j << endl;
		FatalError("Bad operand...");
	}
	operator()(i,j)=val;
}
//=============================================================
void MKLPardisoSymmMatrix::addtoaij(const int i, const int j, double val)
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
double MKLPardisoSymmMatrix::getaij(const int i, const int j)
{
	return operator()(i,j);
}
//=============================================================
void MKLPardisoSymmMatrix::InitializePardisoParameters()
{
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
		for (int i = 0; i < 64; i++) iparm[i] = 0;

//		ompc_set_num_threads(2);
		MKL_Set_Num_Threads(num_threads_for_solver);

		iparm[0] = 1; /* No solver default */	
		if(reorderingScheme==METIS)
			iparm[1] = 2; /* Fill-in reordering from METIS */
		else if(reorderingScheme==MMD)
			iparm[1] = 0; /* Fill-in reordering from MMD */
		else if(reorderingScheme==PAR)
			iparm[1] = 3; /* Fill-in reordering from parallel version of nested disection algorithm */
			
		/* Numbers of processors, value of OMP_NUM_THREADS */
		//omp_set_num_threads(2);
		//iparm[2] = omp_get_max_threads();
		//iparm[2] = omp_get_num_threads();
		//iparm[2] = 2;

		iparm[3] = 0; /* No iterative-direct algorithm */
		iparm[4] = 0; /* No user fill-in reducing permutation */
		iparm[5] = 0; // if 0, Write solution into x. if 1, Write solution into b
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
		if(verbose >= Basic)
			msglvl = 1; /* Print statistical information in file */
		else
			msglvl = 0; /* Print statistical information in file */
		error = 0; /* Initialize error flag */

		for (int i = 0; i < 64; i++) pt[i] = 0;
}
//=============================================================
void MKLPardisoSymmMatrix::printCSRFormatASCII(char * filename, ostream &out)
{
	writeCSRMatrixInASCII(filename, rowindex, colindex, coefs, numberOfEquations, ncoefs);
}
//=============================================================
void MKLPardisoSymmMatrix::printCSRFormatBinary(char * filename, ostream &out)
{
	writeCSRMatrixInBinary(filename, rowindex, colindex, coefs, numberOfEquations, ncoefs);
}
//=============================================================
void MKLPardisoSymmMatrix::solve(double *b,double *initialGuess) 
{
	// Note: If all the dofs are constrained, i.e. ONLY diagonal is nonzero,
	// then PARDISO does NOT print the 'solver report' to cout, but it does 'solve' 
	// the system of equations, a zero vector is returned.

	// Check for positive definite
	CheckAndFixZeroDiagonal();

//	writeCSRMatrixInBinary("CSR_binary.txt", rowindex, colindex, coefs, numberOfEquations, ncoefs);
//	readCSRMatrixInBinary("CSR_binary.txt", rowindex, colindex, coefs, numberOfEquations, ncoefs);
//	writeCSRMatrixInASCII("CSR_ascii.txt", rowindex, colindex, coefs, numberOfEquations, ncoefs);

	WhoMethod("solve(double *)");
	register int i=0;
	const int L=numberOfEquations;
	DEBUG3(printMatrix("DEBUG is ON"));

	mtype = 2; /* Real PSD symmetric matrix */
	nrhs = 1; /* Number of right hand sides. */

	if(OrderingCompleted == false ){	

		InitializePardisoParameters(); 

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
		OrderingCompleted=true;
		if(verbose >= Basic){
			printf("\nReordering completed ... ");
			printf("\nNumber of nonzeros in factors = %d", iparm[17]);
			printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
		}
	}
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
	if(FactorizeCompleted ==false || UsePreviousFactorizedMatrix == false){
		phase = 22;
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&numberOfEquations, coefs, rowindex, colindex, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			printf("\nERROR during numerical factorization: %d", error);
			exit(2);
		}
		if(verbose >= Basic){
			printf("\nFactorization completed ... ");
		}
		FactorizeCompleted=true;
	}
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */
	if(pardiso_x==0) //this array is required by pardiso when solving (back substitution)
		Assert(pardiso_x = new double [numberOfEquations]);

	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&numberOfEquations, coefs, rowindex, colindex, &idum, &nrhs,
		iparm, &msglvl, b, pardiso_x, &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
//test
	for(i=0;i<L;i++) b[i]=pardiso_x[i];
	delete [] pardiso_x;pardiso_x=0;
//

	if(UsePreviousFactorizedMatrix == false && UsePreviousFillInOrdering == false){
		ReleaseMemory();
	}

	if(UsePreviousFactorizedMatrix == false && UsePreviousFillInOrdering == true){
		FactorizeCompleted = false;
	}
}
//=============================================================
void MKLPardisoSymmMatrix::zeroYourself(void)
{
	int i; for(i=0;i<ncoefs;i++) coefs[i]=0;
    FactorizeCompleted = false;
}
//=============================================================
void MKLPardisoSymmMatrix::printMatrix(char *message)
{
	const int L=numberOfEquations;
	cout << message << "\nMKLPardisoSymmMatrix with " << numberOfEquations <<
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
void MKLPardisoSymmMatrix::printMatrixAddressFormat(char * message, ostream* out)
{ 
	int numNonZeros=0;
	int i;
	const int L=numberOfEquations;
	//*out << "\nSymmetricSparseMatrix\n";
	*out << "2" << endl; //shows that it is upper triangular.
	*out << numberOfEquations << " " << ncoefs << endl;

	int ii=0, jj;

	for(i=0,ii=0; i<ncoefs; i++) {
		if(i==rowindex[ii+1]-1){
			*out << '\n';
			ii++;
		}
		jj=colindex[i]-1;
		*out << '\t' << ii << " " << jj << " " << DOUBLE_FORMAT<< coefs[i];
	}
	*out << endl;
}
//=================================================================
void MKLPardisoSymmMatrix::printCompressedRowFormat(char * message, ostream* out)
{ 
	//prints row major upper triangular storage format
	cout << "\nPrinting row major Upper triangular storage format" << endl;
	*out << "2-Upper triangular storage format" <<endl;  //indicates Upper triangular storage format

	int numNonZeros=0;
	int i;
	const int L=numberOfEquations;
	*out << numberOfEquations << endl;

	int *numnz=0;
	numnz=new int[L];

	int ii=0, jj;
	//calculate num non zeros in each row
	int nzctr=0,rowctr=0;
	for(i=0; i<ncoefs; i++) {
		if(i==rowindex[ii+1]-1){
			ii++;
			numnz[rowctr]=nzctr;
			nzctr=0;
			rowctr++;
		}
		nzctr++;
	}
	numnz[L-1]=nzctr;

	//writing to file
	*out << numnz[0];
	ii=0, jj=0;
	for(i=0; i<ncoefs; i++) {
		if(i==rowindex[ii+1]-1){
			ii++;
			*out << endl << numnz[ii];
		}
		jj=colindex[i]-1;
		*out << ' ' << DOUBLE_FORMAT << coefs[i] << ' ' << jj ;
	}
	*out << endl;
	delete [] numnz;
}
//=================================================================
void MKLPardisoSymmMatrix::CopyDataFromInterface(SparseMatrixInterface *smi)
{
	//copying information from the interface
	numberOfEquations=smi->numEquations;
	ncoefs=smi->numNonZeros;
	const int L=numberOfEquations;
	Assert(coefs = new double [ncoefs]);
	Assert(rowindex	= new int[numberOfEquations+1]);
	Assert(colindex	= new int[ncoefs]);

	rowindex[numberOfEquations]	=ncoefs + 1;
	rowindex[0]		=1;

	int i,j,numnz;
	int ctr=0;
	int ii,jj;

	if(smi->type==LowerTriangular){
	//this is converting lower triangle to upper triangular
		for(i=0;i<L;i++) {
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
	}else if(smi->type==UpperTriangular || smi->type==Unsymmetric){
	//this is copying upper triangular information into this pardiso matrix object
		for(i=0;i<L;i++) {
			rowindex[i]=ctr+1;
			numnz = smi->row[i].getNum();
			for(j=0;j<numnz;j++) {
				coefs[ctr]		= smi->mat[i][j];
				colindex[ctr]	= smi->row[i][j]+1;
				ctr++;
			}
		}
	}
}
//=================================================================
void MKLPardisoSymmMatrix::readMatrixAddressFormat(istream* instream)
{ 
	/*
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
		//need to fix this make another class for unsymmetric matrix
		if(ii != jj)
			specifyNonZeroLocation(jj,ii);
	}
	allocate();

	instream->seekg(m_posStartData);

	for(i=0;i<numNonZeros;i++) {
		*instream >> ii >> jj >> value;
		//this block if symmetric matrix
		if(ii == jj)
			MatSetValue(mat, ii, jj, value, ADD_VALUES);
		if(ii > jj)
			MatSetValue(mat, jj, ii, value, ADD_VALUES);//stores upper triangle
		//=========
	}
*/

}
//=============================================================
void MKLPardisoSymmMatrix::vectorProduct(double *result, double *vector)
{
	/*
	int i,j,nt;
	double *mRow;
	double vector_i,tmpResult_i;
	const int L=numberOfEquations;
	SortedList<int> *Row = row;
	double **Mat = mat;
	
	for(i=0;i<L;i++) result[i] = 0;
	for(i=0;i<L;i++) {
		nt = Row[i].getNum()-1;
		if(nt>=0) { 
			SortedList<int> &rowColNum = Row[i];
			vector_i = vector[i];
			mRow = Mat[i];
			tmpResult_i = 0;
			for(j=0;j<nt;j++) {
				tmpResult_i += mRow[j]*vector[rowColNum[j]];
				result[rowColNum[j]] += mRow[j]*vector_i;
			}
			result[i] += mRow[nt]*vector[rowColNum[nt]] + tmpResult_i;
			if(rowColNum[nt]!=i) result[rowColNum[nt]] += mRow[nt]*vector_i;
		}
	}
	*/
}
//=============================================================
void MKLPardisoSymmMatrix::specifyNonZeroLocation(const int I, const int J)
{ 
	//row[I].add(J); 
	row[J].add(I); // pardiso stores the upper triangle info
	//but the equations class assumes that the LargeMatrix stores the lower diagonal info,
	//so the indices need to be switched.
}
//=============================================================
void MKLPardisoSymmMatrix::allocate(void)
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
	delete [] row; 
	row = NULL;
	zeroYourself();
}
//=============================================================
void MKLPardisoSymmMatrix::copyFrom(LargeMatrix *a)
{
	//phh 07/18/08 make a deep copy of the coefs, colindex
	//and rowindex
	ncoefs = ((MKLPardisoSymmMatrix*)a)->ncoefs;
	Assert(coefs = new double [ncoefs]);
	Assert(colindex = new int [ncoefs]);
	Assert(rowindex = new int [numberOfEquations+1]);

	for(int i = 0; i < ncoefs; i++){
		colindex[i] = ((MKLPardisoSymmMatrix*)a)->colindex[i];
	}

	for(int i = 0; i < numberOfEquations+1; i++){
		rowindex[i] = ((MKLPardisoSymmMatrix*)a)->rowindex[i];
	}
}
//=============================================================
void MKLPardisoSymmMatrix::add(LargeMatrix *a, 	bool *restraint)
{
	for(int i = 0; i < ncoefs; i++)
	{
		coefs[i] += ((MKLPardisoSymmMatrix*)a)->coefs[i];
	}
}
//=============================================================
#ifdef DEBUG_MKLPardisoSymmMatrix
ofstream Logfile;
void main(void)
{
	int i;
	ifstream is;
	ofstream os;
	is.open("GlobalStiffness.txt");
	int NumEquations;
	is >> NumEquations;
	MKLPardisoSymmMatrix A(NumEquations);
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
