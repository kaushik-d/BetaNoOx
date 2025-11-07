#include "stdafx.h"
#ifdef PARCG

//#define DEBUG_Level3
//#define PR
#include <iostream>
#include <iomanip>
#include "ParCG.hpp"
#include "utility/excepts.hpp"
#include "utility/progressIndicator.hpp"
#include "sparseMatrixInterface.hpp"
#include <cmath>
#include <omp.h>

CreateErrorHandler(ParCGMatrix);

ParCGMatrix::ParCGMatrix(const int numEqns) :
		SymmetricSparseMatrix(numEqns)
{
	WhoClass("ParCGMatrix");
	WhoMethod("ParCGMatrix(const int)");
}

ParCGMatrix::~ParCGMatrix()
{
	DeleteWhoMethod("~ParCGMatrix()");
//	cerr << "Deleting ParCGMatrix" << endl;
}

void ParCGMatrix::binaryWrite(FILE *outfile)
{
	WhoMethod("binaryWrite(FILE *)");
	int i;
	const int L=numberOfEquations;
	fwrite(&numberOfEquations,sizeof(int),1,outfile);
	for(i=0;i<L;i++) row[i].fwrite(outfile);
	for(i=0;i<L;i++) fwrite(mat[i],sizeof(double),row[i].getNum(),outfile);
}
void ParCGMatrix::binaryRead(FILE *infile)
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

double &ParCGMatrix::operator()(const int ii,
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

void ParCGMatrix::setaij(const int i, const int j, double val)
{
	operator()(i,j)=val;
}
double ParCGMatrix::getaij(const int i, const int j)
{
	return operator()(i,j);
}

#define FORMAT setw(30)<<setprecision(20)


void ParCGMatrix::solve(double *b,double *initialGuess) 
{
	WhoMethod("solve(double *)");
//	ProgressIndicator pi;
	register int i;

	const int L=numberOfEquations;
	double *x,*p,*r,*q,*temp,*diag;
	double pAp,rq,a_k,b_k,newrq,rnorm,xnorm;
	bool convergedSolution = false,zeroOnDiagonal = false;
	DEBUG3(printMatrix("DEBUG is ON"));
	// Check for positive definite
	for(i=0;i<L;i++)  // Quick check not a true check.
		if(operator()(i,i) == 0) zeroOnDiagonal = true; 
	if(zeroOnDiagonal) 
	{
		//+DG_Feb2004
		if(zeroDiagonalReplacementValue==0) //+DG_Feb2004
		{
			cerr<<"The values are zero on diagonal - (" << i <<"," << i << ")";
			cerr <<"\nReplacement for Zero Diagonal not set.\nUse ReplaceZeroDiagonal_Option\n\n";
			exit(1);
		}
		Warning("Zero on diagonal.  Will attempt to fix.");
		for(i=0;i<L;i++)
			if(operator()(i,i)==0) 
			{
				cerr << "Changing [" << i << ',' << i << "] to "<< zeroDiagonalReplacementValue<<".\n";
				operator()(i,i) = zeroDiagonalReplacementValue;
			}
	}
	//-DG_Feb2004
	cout << "MAX_NUMBER_ITERATIONS " << MAX_NUMBER_ITERATIONS << endl;
	cout << "SOLUTION_EPSILON " << SOLUTION_EPSILON << endl;


// 	cout << "Initial Guess:" << endl;
//	for(i=0;i<L;i++) cout << i << " : " << initialGuess[i] << endl;

	// Create solution and residual vector
	Assert(x = new double [numberOfEquations]);
	Assert(r = new double [numberOfEquations]);

	// Create vector to hold diagonal terms
	Assert(diag = new double [numberOfEquations]);
	
// TEMPORARY !!!!
cout << "Setting num threads." << endl;
omp_set_num_threads(2);
cout << "finished setting num threads." << endl;

#pragma omp parallel
{
#pragma omp critical
printf("I am thread # %d of %d threads.\n", omp_get_thread_num(), omp_get_num_threads() );
}

int j;
rowcount = new int [L];
colindex = new int * [L];
for(i=0;i<L;i++){
	rowcount[i]=row[i].getNum();
	colindex[i]=new int [rowcount[i]];
	for(j=0;j<rowcount[i];j++){
		colindex[i][j]=row[i][j];
	}
	
}
	
double **K=mat;
int * Rowcount=rowcount;
//#pragma omp parallel for shared(diag, K, Rowcount) private(i)
	for(i=0;i<L;i++) diag[i] = K[i][Rowcount[i]-1];


	// Set initial guess = {0}
	if(initialGuess==0) 
//#pragma omp parallel for shared(x) private(i) 
		for(i=0;i<L;i++) x[i]=0;
	else {
		cerr << "I have the initial Guess." << endl;
//#pragma omp parallel for shared(x, initialGuess) private(i) 
		for(i=0;i<L;i++) x[i] = initialGuess[i];
	}
	
	// p0 = r = b - A x0
	vectorProduct(r,x);
	delete [] x;

	Assert(p = new double [numberOfEquations]);
	Assert(q = new double [numberOfEquations]);
	Assert(temp = new double [numberOfEquations]);


//#pragma omp parallel for shared(p, q, r, b, diag) private(i) 
	for(i=0;i<L;i++) {
		r[i] = b[i] - r[i];

	// q[i] = r[i]/mat[i][i]
	// Diagonal Scaling Preconditioner
		p[i] = q[i] = r[i]/diag[i];
	}
	
	// Add other preconditioners here

	// Start Iteration
	{
		double firstnorm=0,lastnorm;
		int numIterations = 0;
		
//		pi.start("Solve Status:");
		// Set initial guess = {0}
		double start=omp_get_wtime();				
		if(initialGuess==0) 
//#pragma omp parallel for shared(b) private(i) 
			for(i=0;i<L;i++) b[i] = 0;
		else 
//#pragma omp parallel for shared(b, initialGuess) private(i) 
			for(i=0;i<L;i++) b[i] = initialGuess[i];

		//double end=omp_get_wtime();
		//printf("time = %.16g\n",end - start);
			
int L1=L/2;

//		MAX_NUMBER_ITERATIONS=5;
		while(!convergedSolution && numIterations < MAX_NUMBER_ITERATIONS) {
			// Find pAp = p(traspose) A p    and  rq = r(transpose) q
//			double start3=omp_get_wtime();				
			vectorProduct(temp,p);
//			double start2=omp_get_wtime();				
//			printf("time matvec= %.16g\n",start2 - start3);
			
			pAp=0;rq=0;
//			double pap1=0,pap2=0,rq1=0,rq2=0;
#pragma omp parallel shared(p,temp,r,q) \
		     private(i) firstprivate(L, L1) reduction(+:pAp,rq)
{
if(omp_get_thread_num() == 0){
			for(i=0;i<L1;i++) {
				pAp = pAp + p[i] * temp[i];
				rq  = rq  + r[i] * q[i];
			}
}
if(omp_get_thread_num() == 1){
			for(i=L1;i<L;i++) {
				pAp = pAp + p[i] * temp[i];
				rq  = rq + r[i] * q[i];
			}
}
}
//pAp=pap1 + pap2;
//rq =rq1 + rq2;

			
//			cout << "pAp " << pAp << endl;
//			cout << "a_k " << a_k << endl;
			a_k = rq/pAp;
			// x_k+1[i] = x_k[i] + a_k*p[i]
			newrq=0;
//			double newrq1=0,newrq2=0;
#pragma omp parallel shared(b, p, temp, r, q, diag) \
		     private(i) firstprivate(a_k, L, L1) reduction(+:newrq)
{
if(omp_get_thread_num() == 0){
			for(i=0;i<L1;i++) {
				b[i] += a_k*p[i];
				r[i] -= a_k*temp[i];
				q[i] = r[i]/diag[i];
				newrq = newrq + r[i]*q[i];
			}
}
if(omp_get_thread_num() == 1){
			for(i=L1;i<L;i++) {
				b[i] += a_k*p[i];
				r[i] -= a_k*temp[i];
				q[i] = r[i]/diag[i];
				newrq = newrq + r[i]*q[i];
			}
}
}
//newrq = newrq1 + newrq2;


			b_k = newrq/rq;
	
			// Calculate Error
			rnorm=0;
			xnorm=0;
//			double xnorm1=0,xnorm2=0,rnorm1=0,rnorm2=0;
#pragma omp parallel shared(r,b,p,q) \
		     private(i) firstprivate(b_k, L, L1) reduction(+:xnorm,rnorm)
{
if(omp_get_thread_num() == 0){
			for(i=0;i<L1;i++) {
				p[i] = q[i] + b_k*p[i];
				rnorm = rnorm + r[i]*r[i];
				xnorm = xnorm + b[i]*b[i];
			}
}
if(omp_get_thread_num() == 1){
			for(i=L1;i<L;i++) {
				p[i] = q[i] + b_k*p[i];
				rnorm = rnorm + r[i]*r[i];
				xnorm = xnorm + b[i]*b[i];
			}
}
}
//xnorm=xnorm1+xnorm2;
//rnorm=rnorm1+rnorm2;


//			double end2=omp_get_wtime();
//			printf("time = %.16g\n",end2 - start2);
			

//			cout << "norm " << FORMAT <<  rnorm << " " << xnorm << endl;

			numIterations++;
			if(numIterations%100==0) {
				//cout << numIterations << " - Residual Norm = " << rnorm << " " << xnorm << endl;
				cout << numIterations << " - Residual Norm = " << sqrt(rnorm/xnorm) << endl;
				//getchar();
//				FILE *fp; fp = fopen("sparse.bin","wb");
//				fwrite(b,sizeof(double),L,fp);
//				fclose(fp);
			}
					
			if(firstnorm==0) firstnorm = sqrt(rnorm/xnorm);
			lastnorm = sqrt(rnorm/xnorm);
//cout<< "Residual Norm = " << lastnorm << endl;
//			pi.notifyPercentDone((float) (100.0-100.0*(log(lastnorm)-log(SOLUTION_EPSILON))/
//			      (log(firstnorm)-log(SOLUTION_EPSILON)+1)) );
			if(sqrt(rnorm/xnorm)<SOLUTION_EPSILON) convergedSolution = true;
		}
		if(numIterations==MAX_NUMBER_ITERATIONS) {
			cout << "Solution not converged...\n";
			cerr << "Solution not converged..." << endl;
			exit(1);
		}
cout <<"Iterations for convergence="<< numIterations << " - Residual Norm = " << rnorm << " " << xnorm << endl;
//cout <<"Iterations for convergence="<< numIterations << " - Residual Norm = " << sqrt(rnorm/xnorm) << endl;
//		pi.done();
		double end=omp_get_wtime();
		printf("time = %.16g\n",end - start);

	}
//	cout << "Solution:" << endl;
//	for(i=0;i<L;i++) cout << i << " : " << b[i] << endl;
	delete [] q;
	delete [] p;
	delete [] temp;
	delete [] r;
	delete [] diag; // xt 062901
}


void ParCGMatrix::zeroYourself(void)
{
	int i,j,limit;
	const int L=numberOfEquations;
	for(i=0;i<L;i++) for(j=0,limit=row[i].getNum();j<limit;j++) mat[i][j]=0;
}

void ParCGMatrix::printMatrix(char *message)
{
	int i,j,limit;
	const int L=numberOfEquations;
	cout << message << "\nParCGMatrix with " << numberOfEquations <<
		" equations.\n";
	for(i=0;i<L;i++) {
		for(j=0,limit=row[i].getNum();j<limit;j++)
			cout << '\t' << row[i][j] << ':' << mat[i][j];
		cout << '\n';
	}
}

void ParCGMatrix::printMatrixAddressFormatUpper(char * message, ostream* out)
{ 
	//prints row major UPPER triangular storage format
	//*out << message <<endl; 
	cout << "\nPrinting row major Upper triangular storage format" << endl;
	*out << "2" <<endl; //indicates upper triangular storage format

	int numNonZeros=0;
	int i,limit;
	const int L=numberOfEquations;
	//*out << "\nParCGMatrix\n";
	for(i=0;i<L;i++) {
		limit=row[i].getNum();
		numNonZeros+=limit;
	}

	*out << numberOfEquations << " " << numNonZeros << endl;

#define FORMAT setw(30)<<setprecision(20)

	int ii, jj;

	for(i=0;i<L;i++) {
		//cout << i << endl;
		for(ii=i;ii<L;ii++) {
			//cout << '\t' << ii << endl;
			jj = row[ii].find(i);
			if(jj==row[ii].getNum()) continue;
			*out << '\t' << i << " " << ii << " " << FORMAT<< mat[ii][jj];
		}
		*out << '\n';
	}
}
#undef FORMAT 

void ParCGMatrix::printMatrixAddressFormat(char * message, ostream* out)
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
	//*out << "\nParCGMatrix\n";
	for(i=0;i<L;i++) {
		limit=row[i].getNum();
		numNonZeros+=limit;
	}

	*out << numberOfEquations << " " << numNonZeros << endl;

#define FORMAT setw(30)<<setprecision(20)

	for(i=0;i<L;i++) {
		for(j=0,limit=row[i].getNum();j<limit;j++)
			//if(mat[i][j] != 0){
				*out << '\t' << i << " " << row[i][j] << " " << FORMAT<< mat[i][j];
			//}else{
			//	*out << '\t' << "\n\nZERROOOOOOOOOO " << i << " " << row[i][j] << " " << FORMAT<< mat[i][j];
			//}
		*out << '\n';
	}
}
#undef FORMAT 
//=================================================================
void ParCGMatrix::CopyDataFromInterface(SparseMatrixInterface *smi)
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
		cerr << "Unacceptable format when reading ParCGMatrix!!" << endl;
		cerr << "void ParCGMatrix::CopyDataFromInterface(SparseMatrixInterface *smi)" << endl;
		exit(1);
	}
}
//=================================================================
void ParCGMatrix::readMatrixAddressFormat(istream* instream)
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

void ParCGMatrix::vectorProduct(double *result, double *vector)
{
	int i,j,nt;
	double *mRow;
	double vector_i,tmpResult_i,tmp_Result;
	const int L=numberOfEquations;
//	SortedList<int> *Row = row;
	double **Mat = mat;
	int tmp_index;
	
	int * Rowcount=rowcount;
	int ** Colindex=colindex;

#pragma omp parallel for shared(result) private(i)
	for(i=0;i<L;i++) result[i] = 0;
	
int L1=L/2;
int L2=L-L1;
//cout << "before first section" << endl;

#pragma omp parallel shared(result,Rowcount,Mat,vector,Colindex, L1, L2, L) \
		     private(i, j, nt, vector_i, mRow, tmpResult_i, tmp_Result, tmp_index)
{
if(omp_get_thread_num() == 0){
	for(i=0;i<L1;i++) {
		nt = Rowcount[i]-1;
		if(nt>=0) { 
			vector_i = vector[i];
			mRow = Mat[i];
			tmpResult_i = 0;
			for(j=0;j<nt;j++) {
				tmpResult_i += Mat[i][j]*vector[Colindex[i][j]];
				tmp_Result=mRow[j]*vector_i;
				tmp_index=Colindex[i][j];
//#pragma omp critical
//				cout << "first " << tmp_index << endl;
				result[tmp_index] += tmp_Result;
			}
			tmp_Result = mRow[nt]*vector[Colindex[i][nt]] + tmpResult_i;
			result[i] += tmp_Result;
//#pragma omp critical
//			cout << "first " << i << endl;
		}
	}
}
if(omp_get_thread_num() == 1){
	for(i=L1;i<L;i++) {
		nt = Rowcount[i]-1;
		if(nt>=0) { 
			vector_i = vector[i];
			mRow = Mat[i];
			tmpResult_i = 0;
			for(j=0;j<nt;j++) {
				if( Colindex[i][j] >= L1){
					tmpResult_i += Mat[i][j]*vector[Colindex[i][j]];
					tmp_Result=mRow[j]*vector_i;
					tmp_index=Colindex[i][j];
//#pragma omp critical
//					cout << "second " << tmp_index << endl;
					result[tmp_index] += tmp_Result;
				}
			}
			tmp_Result = mRow[nt]*vector[Colindex[i][nt]] + tmpResult_i;
//#pragma omp critical
//			cout << "second " << i << endl;
			result[i] += tmp_Result;
		}
	}
}
}

//cout << "before second section" << endl;

#pragma omp parallel shared(result,Rowcount,Mat,vector,Colindex) \
		     private(i, j, nt, vector_i, mRow, tmpResult_i, tmp_Result, tmp_index)
{
if(omp_get_thread_num() == 0){
	for(i=L1;i<L;i++) {
		nt = Rowcount[i]-1;
		if(nt>=0) { 
			vector_i = vector[i];
			mRow = Mat[i];
			tmpResult_i = 0;
			for(j=0;j<nt;j++) {
				if( Colindex[i][j] < L1){			
					tmp_Result=mRow[j]*vector_i;
					tmp_index=Colindex[i][j];
					result[tmp_index] += tmp_Result;
				}
			}
		}
	}
}
if(omp_get_thread_num() == 1){
	for(i=L1;i<L;i++) {
		nt = Rowcount[i]-1;
		if(nt>=0) { 
			vector_i = vector[i];
			mRow = Mat[i];
			tmpResult_i = 0;
			for(j=0;j<nt;j++) {
				if( Colindex[i][j] < L1){
					tmpResult_i += Mat[i][j]*vector[Colindex[i][j]];
				}
			}
			result[i] += tmpResult_i;
		}
	}
}
}


}


void ParCGMatrix::specifyNonZeroLocation(const int I, const int J)
{row[I].add(J);}

void ParCGMatrix::allocate(void)
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

#ifdef DEBUG_ParCGMatrix
ofstream Logfile;
void main(void)
{
	int i;
	ifstream is;
	ofstream os;
	is.open("GlobalStiffness.txt");
	int NumEquations;
	is >> NumEquations;
	ParCGMatrix A(NumEquations);
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
	ParCGMatrix A;
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
#endif
