#include "stdafx.h"

//#define DEBUG_Level3
//#define PR
#include <iostream>
#include <iomanip>
#include "sparse.hpp"
#include "utility/utility_defines.h"
#include "utility/excepts.hpp"
#include "utility/progressIndicator.hpp"
#include "sparseMatrixInterface.hpp"

CreateErrorHandler(SymmetricSparseMatrix);

SymmetricSparseMatrix::SymmetricSparseMatrix(const int numEqns) :
		LargeMatrix(numEqns)
{
	WhoClass("SymmetricSparseMatrix");
	WhoMethod("SymmetricSparseMatrix(const int)");
	Assert(row = new SortedList<int>[numberOfEquations]);
	mat = 0;
	setMAX_NUMBER_ITERATIONS(10000);
	setSOLUTION_EPSILON(1e-4);
}

SymmetricSparseMatrix::~SymmetricSparseMatrix()
{
	DeleteWhoMethod("~SymmetricSparseMatrix()");
//	cerr << "Deleting sparse matrix" << endl;
	int i,L=numberOfEquations;
	if(numberOfEquations>0) {
		if(mat){
			for(i=0;i<L;i++) 
				delete [] mat[i];
			delete [] mat;
		}
		if(row) delete [] row;
	}
}

void SymmetricSparseMatrix::setMAX_NUMBER_ITERATIONS(int max_iter)
{
	MAX_NUMBER_ITERATIONS=max_iter;
}
void SymmetricSparseMatrix::setSOLUTION_EPSILON(double sol_epsilon)
{
	SOLUTION_EPSILON=sol_epsilon;
}

void SymmetricSparseMatrix::binaryWrite(FILE *outfile)
{
	WhoMethod("binaryWrite(FILE *)");
	int i;
	const int L=numberOfEquations;
	fwrite(&numberOfEquations,sizeof(int),1,outfile);
	for(i=0;i<L;i++) row[i].fwrite(outfile);
	for(i=0;i<L;i++) fwrite(mat[i],sizeof(double),row[i].getNum(),outfile);
}
void SymmetricSparseMatrix::binaryRead(FILE *infile)
{
	WhoMethod("binaryRead(FILE *)");
	int i;
	const int L=numberOfEquations;
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
}

double &SymmetricSparseMatrix::operator()(const int ii,
		const int jj)
{
	int j;
	if(ii>=numberOfEquations || jj > ii) {
//_operatorError:
		cerr << "numberOfEquations = " << numberOfEquations << endl;
		WhoMethod("operator()(const unsigned Int,const unsigned Int)");
		cerr << "ii,jj = " << ii << " , " << jj << endl;
		FatalError("Bad operand...");
	}
	j = row[ii].find(jj);
	if(j==row[ii].getNum()) {
		cerr << jj << " not on row " << ii << endl;
		cerr << "Row " << ii << "-" << row[ii].getNum();
		row[ii].print(cerr) << endl;
		FatalError("Bad operand...");
	}
	return mat[ii][j];
}

void SymmetricSparseMatrix::setaij(const int i, const int j, double val)
{
	operator()(i,j)=val;
}
double SymmetricSparseMatrix::getaij(const int i, const int j)
{
	return operator()(i,j);
}


#include <cmath>


void SymmetricSparseMatrix::solve(double *b,double *initialGuess) 
{
	WhoMethod("solve(double *)");
	ProgressIndicator pi;
	register int i;

	const int L=numberOfEquations;
	double *x,*p,*r,*q,*temp,*diag;
	double pAp,rq,a_k,b_k,newrq,rnorm,xnorm;
	bool convergedSolution = false;
	DEBUG3(printMatrix("DEBUG is ON"));

	// Check for positive definite
	CheckAndFixZeroDiagonal();

	if(verbose >= Basic)	{
		cout << "MAX_NUMBER_ITERATIONS " << MAX_NUMBER_ITERATIONS << endl;
		cout << "SOLUTION_EPSILON " << SOLUTION_EPSILON << endl;
	}


// 	cout << "Initial Guess:" << endl;
//	for(i=0;i<L;i++) cout << i << " : " << initialGuess[i] << endl;

	// Create solution and residual vector
	Assert(x = new double [numberOfEquations]);
	Assert(r = new double [numberOfEquations]);

	// Create vector to hold diagonal terms
	Assert(diag = new double [numberOfEquations]);
	for(i=0;i<L;i++) diag[i] = mat[i][row[i].getNum()-1];

	// Set initial guess = {0}
	if(initialGuess==0) for(i=0;i<L;x[i++]=0);
	else {
		if(verbose >= Basic)	cerr << "I have the initial Guess." << endl;
		for(i=0;i<L;i++) x[i] = initialGuess[i];
	}
	
	// p0 = r = b - A x0
	vectorProduct(r,x);
	delete [] x;

	Assert(p = new double [numberOfEquations]);
	Assert(q = new double [numberOfEquations]);
	Assert(temp = new double [numberOfEquations]);

	for(i=0;i<L;i++) r[i] = b[i] - r[i];

	// q[i] = r[i]/mat[i][i]
	// Diagonal Scaling Preconditioner
	for(i=0;i<L;i++) p[i] = q[i] = r[i]/diag[i];
	// Add other preconditioners here

	// Start Iteration
	{
		double firstnorm=0,lastnorm;
		int numIterations = 0;
		if(verbose >= Basic)	
			pi.start("Solve Status:");
		// Set initial guess = {0}
		if(initialGuess==0) for(i=0;i<L;i++) b[i] = 0;
		else for(i=0;i<L;i++) b[i] = initialGuess[i];
		while(!convergedSolution && numIterations < MAX_NUMBER_ITERATIONS) {
			// Find pAp = p(traspose) A p    and  rq = r(transpose) q
			vectorProduct(temp,p);
			for(i=0,pAp=0,rq=0;i<L;i++) {
				pAp += p[i]*temp[i];
				rq += r[i]*q[i];
			}
			a_k = rq/pAp;
			// x_k+1[i] = x_k[i] + a_k*p[i]
			for(i=0;i<L;i++) b[i] += a_k*p[i];
	
			// Calculate new residual
			for(i=0;i<L;i++) r[i] -= a_k*temp[i];
		
			// Calculate new q
			for(i=0;i<L;i++) q[i] = r[i]/diag[i];
	
			// Calculate new rq
			for(i=0,newrq=0;i<L;i++) newrq += r[i]*q[i];
			b_k = newrq/rq;
			// Calculate new p
			for(i=0;i<L;i++) p[i] = q[i] + b_k*p[i];
	
			// Calculate Error
			// rnorm = sqrt(r.L2norm());
			for(i=0,rnorm=0;i<L;i++) rnorm += r[i]*r[i];
			// xnorm = sqrt(b.L2norm());
			for(i=0,xnorm=0;i<L;i++) xnorm += b[i]*b[i];
			numIterations++;
			if(numIterations%100==0) {
				if(verbose >= Basic)
					cout << numIterations << " - Residual Norm = " << sqrt(rnorm/xnorm) << endl;
				//getchar();
				/* //JV090507 commenting output for now. when needed, uncomment it
				FILE *fp; fp = fopen("sparse.bin","wb");
				fwrite(b,sizeof(double),L,fp);
				fclose(fp);
				*/
			}
					
			if(firstnorm==0) firstnorm = sqrt(rnorm/xnorm);
			lastnorm = sqrt(rnorm/xnorm);
//cout<< "Residual Norm = " << lastnorm << endl;
		if(verbose >= Basic)
			pi.notifyPercentDone((float) (100.0-100.0*(log(lastnorm)-log(SOLUTION_EPSILON))/
			      (log(firstnorm)-log(SOLUTION_EPSILON)+1)) );
			if(sqrt(rnorm/xnorm)<SOLUTION_EPSILON) convergedSolution = true;
		}
		if(numIterations==MAX_NUMBER_ITERATIONS) {
			cout << "Solution not converged...\n";
			cerr << "Solution not converged..." << endl;
			exit(1);
		}
		if(verbose >= Basic){
			cout <<"Iterations for convergence="<< numIterations << " - Residual Norm = " << sqrt(rnorm/xnorm) << endl;
			pi.done();
		}
	}
//	cout << "Solution:" << endl;
//	for(i=0;i<L;i++) cout << i << " : " << b[i] << endl;
	delete [] q;
	delete [] p;
	delete [] temp;
	delete [] r;
	delete [] diag; // xt 062901
}
void SymmetricSparseMatrix::zeroYourself(void)
{
	int i,j,limit;
	const int L=numberOfEquations;
	for(i=0;i<L;i++) for(j=0,limit=row[i].getNum();j<limit;j++) mat[i][j]=0;
}

void SymmetricSparseMatrix::printMatrix(char *message)
{
	int i,j,limit;
	const int L=numberOfEquations;
	cout << message << "\nSymmetricSparseMatrix with " << numberOfEquations <<
		" equations.\n";
	for(i=0;i<L;i++) {
		for(j=0,limit=row[i].getNum();j<limit;j++)
			cout << '\t' << row[i][j] << ':' << mat[i][j];
		cout << '\n';
	}
}

void SymmetricSparseMatrix::printMatrixAddressFormatUpper(char * message, ostream* out)
{ 
	//prints row major UPPER triangular storage format
	//*out << message <<endl; 
	cout << "\nPrinting row major Upper triangular storage format" << endl;
	*out << "2" <<endl; //indicates upper triangular storage format

	int numNonZeros=0;
	int i,limit;
	const int L=numberOfEquations;
	//*out << "\nSymmetricSparseMatrix\n";
	for(i=0;i<L;i++) {
		limit=row[i].getNum();
		numNonZeros+=limit;
	}

	*out << numberOfEquations << " " << numNonZeros << endl;

	int ii, jj;

	for(i=0;i<L;i++) {
		//cout << i << endl;
		for(ii=i;ii<L;ii++) {
			//cout << '\t' << ii << endl;
			jj = row[ii].find(i);
			if(jj==row[ii].getNum()) continue;
			*out << '\t' << i << " " << ii << " " << DOUBLE_FORMAT<< mat[ii][jj];
		}
		*out << '\n';
	}
}

void SymmetricSparseMatrix::printMatrixAddressFormat(char * message, ostream* out)
{ 
	//temporary call to write upper triangular format
	//return printMatrixAddressFormatUpper(message, out);

	//prints row major lower triangular storage format
	cout << "\nPrinting row major Lower triangular storage format" << endl;
	//*out << message <<endl; 
	*out << "1" <<endl;  //indicates Lower triangular storage format

	int numNonZeros=0;
	int i,j,limit;
	const int L=numberOfEquations;
	//*out << "\nSymmetricSparseMatrix\n";
	for(i=0;i<L;i++) {
		limit=row[i].getNum();
		numNonZeros+=limit;
	}

	*out << numberOfEquations << " " << numNonZeros << endl;

	for(i=0;i<L;i++) {
		for(j=0,limit=row[i].getNum();j<limit;j++)
			//if(mat[i][j] != 0){
				*out << '\t' << i << " " << row[i][j] << " " << DOUBLE_FORMAT<< mat[i][j];
			//}else{
			//	*out << '\t' << "\n\nZERROOOOOOOOOO " << i << " " << row[i][j] << " " << FORMAT<< mat[i][j];
			//}
		*out << '\n';
	}
}
//=================================================================
void SymmetricSparseMatrix::printCompressedRowFormat(char * message, ostream* out)
{ 
	//temporary call to write upper triangular format
	//return printMatrixAddressFormatUpper(message, out);

	//prints row major lower triangular storage format
	cout << "\nPrinting row major Lower triangular storage format" << endl;
	*out << "1-Lower triangular storage format" <<endl;  //indicates Lower triangular storage format

	int numNonZeros=0;
	int i,j,limit;
	const int L=numberOfEquations;
	//*out << "\nSymmetricSparseMatrix\n";
	for(i=0;i<L;i++) {
		limit=row[i].getNum();
		numNonZeros+=limit;
	}

	*out << numberOfEquations << endl;

	int numdiag=0;
	for(i=0;i<L;i++) {
		*out << row[i].getNum();
		for(j=0,limit=row[i].getNum();j<limit;j++){
			*out << " " << DOUBLE_FORMAT<< mat[i][j] << " " << row[i][j];
			if(i==row[i][j]) numdiag++;
		}
		*out << endl;
	}
	cout << "numdiag = " << numdiag << endl;
}
//=================================================================
void SymmetricSparseMatrix::printCompressedRowFormatUpper(char * message, ostream* out)
{ 
	//prints row major UPPER triangular storage format
	cout << "\nPrinting row major Upper triangular storage format" << endl;
	*out << "2-Upper triangular storage format" <<endl;  //indicates Lower triangular storage format

	int numNonZeros=0;
	int i,limit;
	const int L=numberOfEquations;
	for(i=0;i<L;i++) {
		limit=row[i].getNum();
		numNonZeros+=limit;
	}

	*out << numberOfEquations << endl;

/*
	int ii, jj;
	for(i=0;i<L;i++) {
		for(ii=i;ii<L;ii++) {
			jj = row[ii].find(i);
			if(jj==row[ii].getNum()) continue;
			*out << '\t' << i << " " << ii << " " << DOUBLE_FORMAT<< mat[ii][jj];
		}
		*out << '\n';
	}
*/
}
//=================================================================
void SymmetricSparseMatrix::CopyDataFromInterface(SparseMatrixInterface *smi)
{
	//copying information from the interface
	numberOfEquations=smi->numEquations;
	const int L=numberOfEquations;
	Assert(mat = new double *[L]);
	Assert(row = new SortedList<int>[L]);
	int i,j,numnz;

	if(smi->type==LowerTriangular){
		for(i=0;i<L;i++) {
			numnz = smi->row[i].getNum();
			//row[i]=smi.row[i]; //no valid copy constructor for this class
			Assert(mat[i] = new double [numnz]);
			for(j=0;j<numnz;j++) {
				row[i].add(smi->row[i][j]); //this is added cos there is not valid copy constructor for this class
				mat[i][j]=smi->mat[i][j];
			}
		}	
	}else if(smi->type==UpperTriangular){
		for(i=0;i<L;i++) {
			numnz = smi->col[i].getNum();
			Assert(mat[i] = new double [numnz]);
			for(j=0;j<numnz;j++) {
				row[i].add(smi->col[i][j]); 
				mat[i][j]=smi->mat[ smi->col[i][j] ][ smi->row[smi->col[i][j]].find(i) ];
			}
		}	
	}else {
		cerr << "Unacceptable format when reading SymmetricSparseMatrix!!" << endl;
		cerr << "void SymmetricSparseMatrix::CopyDataFromInterface(SparseMatrixInterface *smi)" << endl;
		exit(1);
	}
}
//=================================================================
void SymmetricSparseMatrix::readMatrixAddressFormat(istream* instream)
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

void SymmetricSparseMatrix::vectorProduct(double *result, double *vector)
{
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
			if(rowColNum[nt]!=i) 
				result[rowColNum[nt]] += mRow[nt]*vector_i;
		}
	}
}
void SymmetricSparseMatrix::specifyNonZeroLocation(const int I, const int J)
{row[I].add(J);}

void SymmetricSparseMatrix::allocate(void)
{
	int i;
	const int L=numberOfEquations;
	Assert(mat = new double *[numberOfEquations]);
	for(i=0;i<L;i++) {
		if(row[i].getNum()==0) row[i].add(i);
#ifdef PR
		cout << "Row " << i << " - (" << row[i].getNum() << ") - ";
		row[i].print(cout) << endl;
#endif
		Assert(mat[i] = new double [row[i].getNum()]);
	}
}
void SymmetricSparseMatrix::copyFrom(LargeMatrix *a)
{
	int i;
	const int L = numberOfEquations;
	Assert(mat = new double *[numberOfEquations]);
	Assert(row = new SortedList<int>[numberOfEquations]);
	for(i=0;i<L;i++)
	{
		row[i] = ((SymmetricSparseMatrix*)a)->row[i];
		Assert(mat[i] = new double[row[i].getNum()]);
	}
}

void SymmetricSparseMatrix::add(LargeMatrix* a)
{
	for(int i=0;i<numberOfEquations;i++)
	{
		for(int j=0;j<row[i].getNum();j++)
		{
			mat[i][j] += ((SymmetricSparseMatrix*)a)->mat[i][j];
		}
	}
}


#ifdef DEBUG_SymmetricSparseMatrix
ofstream Logfile;
void main(void)
{
	int i;
	ifstream is;
	ofstream os;
	is.open("GlobalStiffness.txt");
	int NumEquations;
	is >> NumEquations;
	SymmetricSparseMatrix A(NumEquations);
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
/*
void main(void)
{
	int i;
	SymmetricSparseMatrix A;
	A.allocateStorage();
	A.printMatrix("Test Mat...");

	double result[4],vector[4]={1,1,1,1};
	A.vectorProduct(result,vector);
	for(i=0;i<4;i++) cout << '\t' << result[i];
	cout << endl;
	A(0,0) = 10;
	A(1,1) = 10;
	A.printMatrix("Test Mat...");
	A.vectorProduct(result,vector);
	for(i=0;i<4;i++) cout << '\t' << result[i];
	cout << endl;

	vector[0] = 18;
	vector[1] = 20;
	vector[2] = 9;
	vector[3] = 13;
	A.solve(vector);
	for(i=0;i<4;i++) cout << '\t' << vector[i];
	cout << endl;
	cout << "Answer should be all 1's" << endl;
}
*/
#endif
