#include "stdafx.h"
//#define WSMP
#ifdef WSMP

#include <cmath>
#include <iostream>
#include <iomanip>
#include "WSMPSymm.hpp"
#include "utility/utility_defines.h"
#include "utility/excepts.hpp"
#include "sparseMatrixInterface.hpp"
//#include <omp.h>

void readMatrix(char* filename, int* &ia, int* &ja, double* &a, int &n);

/* Change this, if your Fortran compiler does not append underscores. */
/* e.g. the AIX compiler:  #define F77_FUNC(func) func                */
#ifdef AIX
#define F77_FUNC(func)  func 
#else
#define F77_FUNC(func)  func ## _
#endif

extern "C" {
	double F77_FUNC(wsmprtc)();
	double F77_FUNC(wsmp_initialize)();
	void   F77_FUNC(wsetmaxthrds)(int*);
	void F77_FUNC(wssmp)
        (int*, int*, int*, double*, double*, int*, int*, 
		double*, int*, int*, double*, int*, int*, int*, double*);
}


CreateErrorHandler(WSMPSymmMatrix);
//=============================================================
WSMPSymmMatrix::WSMPSymmMatrix(const int numEqns) :
		LargeMatrix(numEqns)
{
	WhoClass("WSMPSymmMatrix");
	WhoMethod("WSMPSymmMatrix(const int)");
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
WSMPSymmMatrix::~WSMPSymmMatrix()
{
	DeleteWhoMethod("~WSMPSymmMatrix()");
//	cerr << "Deleting WSMPSymmMatrix matrix" << endl;

	//if(DoReleasePardisoMemory) ReleasePardisoMemory();
	if(pardiso_x) {delete [] pardiso_x; pardiso_x=0;}

	int L=numberOfEquations;
	if(numberOfEquations>0) {
		delete [] row;
	}
	if(coefs) delete [] coefs;
	if(colindex) delete [] colindex;
	if(rowindex) delete [] rowindex;
}
//=============================================================
void WSMPSymmMatrix::binaryWrite(FILE *outfile)
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
void WSMPSymmMatrix::binaryRead(FILE *infile)
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
double& WSMPSymmMatrix::operator()(const int ii,
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
void WSMPSymmMatrix::setaij(const int i, const int j, double val)
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
void WSMPSymmMatrix::addtoaij(const int i, const int j, double val)
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
double WSMPSymmMatrix::getaij(const int i, const int j)
{
	return operator()(i,j);
}
//=============================================================
void WSMPSymmMatrix::solve(double *b,double *initialGuess) 
{
	WhoMethod("solve(double *)");

	// Check for positive definite
	CheckAndFixZeroDiagonal();

	register int i=0;
	const int L=numberOfEquations;
	DEBUG3(printMatrix("DEBUG is ON"));

    int      nrhs = 1;          /* Number of right hand sides. */
    int      iparm[64];
    double	 dparm[64];
    double waltime;
	double diag,aux;
	int ldb=L;
	int naux=0;
	int mrp;
	int *perm, *invp;
	perm=new int[L];
	invp=new int[L];

/* Executable part of the program starts here.*/
        waltime = F77_FUNC(wsmprtc)();
/*
// Fill 'iparm' and 'dparm' arrays with default values.

//  As an alternative to this step, the values in 'iparm' and 'dparm' can be
//   filled with values suitable for the application either manually or by 
//   using a "data" statement according to the description in the User's guide.
*/
	//set num threads
	int nthds=4;
	F77_FUNC(wsetmaxthrds)(&nthds);
	
	F77_FUNC(wsmp_initialize)();
        iparm[0] = 0;
        iparm[1] = 0;
        iparm[2] = 0;
        F77_FUNC(wssmp) (&numberOfEquations, rowindex, colindex, coefs, &diag, perm, invp, b, &ldb, &nrhs, &aux, &naux, &mrp, iparm, dparm);

        printf("\nNumber of CPUs used - %d",iparm[32]);
        printf("\nInitialization complete in time - %f",F77_FUNC(wsmprtc)()-waltime);
		if (iparm[63] != 0) {
          printf("The following ERROR was detected: %d",iparm[63]);
          exit(1);
		}

/* Ordering.*/

        waltime = F77_FUNC(wsmprtc)();
        iparm[1] = 1;
        iparm[2] = 1;
        F77_FUNC(wssmp) (&numberOfEquations, rowindex, colindex, coefs, &diag, perm, invp, b, &ldb, &nrhs, &aux, &naux, &mrp, iparm, dparm);
        printf("\nNumber of CPUs used - %d",iparm[32]);
        printf("\nOrdering complete in time - %f",F77_FUNC(wsmprtc)()-waltime);

		if (iparm[63] != 0) {
          printf("The following ERROR was detected: %d",iparm[63]);
          exit(1);
		}

/* Symbolic Factorization.*/

        waltime = F77_FUNC(wsmprtc)();
        iparm[1] = 2;
        iparm[2] = 2;
        F77_FUNC(wssmp) (&numberOfEquations, rowindex, colindex, coefs, &diag, perm, invp, b, &ldb, &nrhs, &aux, &naux, &mrp, iparm, dparm);

        printf("\nNumber of CPUs used - %d",iparm[32]);
        printf("\nSymbolic complete in time - %f",F77_FUNC(wsmprtc)()-waltime);
		if (iparm[63] != 0) {
          printf("The following ERROR was detected: %d",iparm[63]);
          exit(1);
		}

        printf("\nNumber of nonzeros in factor L = %d",iparm[23]);
        printf("\nNumber of FLOPS in factorization = %f",dparm[22]);
        printf("\nDouble words needed to factor = %d",iparm[22]);

/* Cholesky Factorization.*/

        waltime = F77_FUNC(wsmprtc)();
        iparm[1] = 3;
        iparm[2] = 3;
        F77_FUNC(wssmp) (&numberOfEquations, rowindex, colindex, coefs, &diag, perm, invp, b, &ldb, &nrhs, &aux, &naux, &mrp, iparm, dparm);

        waltime = F77_FUNC(wsmprtc)() - waltime;
        printf("\nNumber of CPUs used - %d",iparm[32]);
        printf("\nCholesky complete in time - %f",waltime);
        printf("\nFactorization MegaFlops = %f",(dparm[22] * 1.e-6) / waltime);
		if (iparm[63] != 0) {
          printf("The following ERROR was detected: %d",iparm[63]);
		  if (iparm[63] > 0) printf("Probably a non spd matrix");
          exit(1);
		}

/* Back substitution.*/
        waltime = F77_FUNC(wsmprtc)();
        iparm[1] = 4;
        iparm[2] = 4;
        F77_FUNC(wssmp) (&numberOfEquations, rowindex, colindex, coefs, &diag, perm, invp, b, &ldb, &nrhs, &aux, &naux, &mrp, iparm, dparm);

        printf("\nNumber of CPUs used - %d",iparm[32]);
        printf("\nBack substitution done in time - %f",F77_FUNC(wsmprtc)()-waltime);
		if (iparm[63] != 0) {
          printf("The following ERROR was detected: %d",iparm[63]);
          exit(1);
		}

/* Itertative refinement.*/

        waltime = F77_FUNC(wsmprtc)();
        iparm[1] = 5;
        iparm[2] = 5;
        F77_FUNC(wssmp) (&numberOfEquations, rowindex, colindex, coefs, &diag, perm, invp, b, &ldb, &nrhs, &aux, &naux, &mrp, iparm, dparm);

        printf("\nNumber of CPUs used - %d",iparm[32]);
        printf("\nItertative ref. done in time - %f",F77_FUNC(wsmprtc)()-waltime);
		if (iparm[63] != 0) {
          printf("The following ERROR was detected: %d",iparm[63]);
          exit(1);
		}

	//memory release
	delete [] perm;
	delete [] invp;
}
//=============================================================
void WSMPSymmMatrix::zeroYourself(void)
{
	int i; for(i=0;i<ncoefs;i++) coefs[i]=0;
}
//=============================================================
void WSMPSymmMatrix::printMatrix(char *message)
{
	const int L=numberOfEquations;
	cout << message << "\nWSMPSymmMatrix with " << numberOfEquations <<
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
void WSMPSymmMatrix::printMatrixAddressFormat(char * message, ostream* out)
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
void WSMPSymmMatrix::printCompressedRowFormat(char * message, ostream* out)
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
void WSMPSymmMatrix::CopyDataFromInterface(SparseMatrixInterface *smi)
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
void WSMPSymmMatrix::readMatrixAddressFormat(istream* instream)
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
void WSMPSymmMatrix::vectorProduct(double *result, double *vector)
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
void WSMPSymmMatrix::specifyNonZeroLocation(const int I, const int J)
{ 
	//row[I].add(J); 
	row[J].add(I); // pardiso stores the upper triangle
	// the element matrices store the lower triangle 
	//so the indices need to be switched.
}
//=============================================================
void WSMPSymmMatrix::allocate(void)
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
#ifdef DEBUG_WSMPSymmMatrix
ofstream Logfile;
void main(void)
{
	int i;
	ifstream is;
	ofstream os;
	is.open("GlobalStiffness.txt");
	int NumEquations;
	is >> NumEquations;
	WSMPSymmMatrix A(NumEquations);
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
