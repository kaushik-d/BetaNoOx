#include "stdafx.h"

#ifdef PETSC
#include <cmath>
#include <iostream>
#include <iomanip>
#include "petscsparse.hpp"
#include "utility/excepts.hpp"
#include "utility/progressIndicator.hpp"

#include "petscksp.h"

CreateErrorHandler(PETScSparseMatrix);
//=============================================================
PETScSparseMatrix::PETScSparseMatrix(const int numEqns) :
		LargeMatrix(numEqns)
{
	WhoClass("PETScSparseMatrix");
	WhoMethod("PETScSparseMatrix(const int)");
	Assert(row = new SortedList<int>[numberOfEquations]);
	//mat = 0;
	//assemblyFlag=0;
	FilledDiagonal=0;
}
//=============================================================
PETScSparseMatrix::~PETScSparseMatrix()
{
	DeleteWhoMethod("~PETScSparseMatrix()");
//	cerr << "Deleting PETScSparseMatrix matrix" << endl;
	int L=numberOfEquations;
	if(numberOfEquations>0) {
		MatDestroy( mat );
		delete [] row;
	}
	PetscFinalize();
	if(FilledDiagonal)
		delete [] FilledDiagonal;
}
//=============================================================
void PETScSparseMatrix::binaryWrite(FILE *outfile)
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
void PETScSparseMatrix::binaryRead(FILE *infile)
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
double& PETScSparseMatrix::operator()(const int ii,
		const int jj)
{
	cerr << "cannot use operator() with PetSc solver!!!" << endl;
	WhoMethod("operator()(const unsigned Int,const unsigned Int)");
	FatalError("cannot use operator() with PetSc solver!!!");
	double val=0;
	return val;
}
//=============================================================
void PETScSparseMatrix::setaij(const int i, const int j, double val)
{
	//this can be used only for setting the diagonal !
	if(i>=numberOfEquations || j>=numberOfEquations ) {
		cerr << "numberOfEquations = " << numberOfEquations << endl;
		WhoMethod("setaij(const int i, const int j, double val)");
		cerr << "ii,jj = " << i << " , " << j << endl;
		FatalError("Bad operand...");
	}
	if(i==j){
		if(FilledDiagonal[i]==false){
			MatSetValue(mat, i, i, val, ADD_VALUES);
			FilledDiagonal[i]=true;
		}
	}else{
		cout << "PETScSparseMatrix::setaij() can be used only to set diagonal terms!" << endl;
		exit(1);
	}
}
//=============================================================
void PETScSparseMatrix::addtoaij(const int i, const int j, double val)
{
	if(i>=numberOfEquations || j>=numberOfEquations ) {
		cerr << "numberOfEquations = " << numberOfEquations << endl;
		WhoMethod("setaij(const int i, const int j, double val)");
		cerr << "ii,jj = " << i << " , " << j << endl;
		FatalError("Bad operand...");
	}
	if(i>j)
		MatSetValue(mat, j, i, val, ADD_VALUES);
	else
		MatSetValue(mat, i, j, val, ADD_VALUES);
}
//=============================================================
double PETScSparseMatrix::getaij(const int i, const int j)
{
	PetscScalar val;
	MatGetValues(mat, 1, &i, 1, &j, &val);
	return (double)val;
}
//=============================================================
void PETScSparseMatrix::solve(double *b,double *initialGuess) 
{
	WhoMethod("solve(double *)");
	ProgressIndicator pi;
	register int i=0;
//	int j,k;
	const int L=numberOfEquations;
	bool convergedSolution = false,zeroOnDiagonal = false;
	DEBUG3(printMatrix("DEBUG is ON"));
	// Check for positive definite
	/**************************************************************************************
	/* JV081406 - just checking for now... operator function needs to be fixed
	/**************************************************************************************
	for(i=0;i<L;i++)  // Quick check not a true check.
		if(operator()(i,i) == 0) zeroOnDiagonal = true; 
	*/
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

// 	cout << "Initial Guess:" << endl;
//	for(i=0;i<L;i++) cout << i << " : " << initialGuess[i] << endl;


		pi.start("Solve Status:");
  MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY); 
  MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

  double iter_error = 1e-15;
  int max_iter_num = 1000000;
  int num_of_iter;

  PetscErrorCode ierr;

	Vec rhs, x;
    VecCreateSeqWithArray(PETSC_COMM_SELF, L, b, &rhs);
    VecDuplicate(rhs, &x);

    KSP ksp;
    KSPCreate(PETSC_COMM_SELF, &ksp);

	//this block is for symmetric matrices
	KSPSetType(ksp,KSPCG); // to set CG type - symmetric matrix solver !
	PC	pc;	//pre-conditioner
	KSPGetPC(ksp,&pc);
	PCSetType(pc,PCICC); //ILU is for general matrices; for symmetric you need to use PCSetType(pc,PCICC); or -pc_type icc
	//PCSetType(pc,PCJACOBI); //Jacobi (i.e. diagonal scaling preconditioning)
	//PCJacobiSetUseAbs(pc);
	//PCJacobiSetUseRowMax(pc);
	//===============

    KSPSetTolerances(ksp, iter_error, PETSC_DEFAULT, 
      PETSC_DEFAULT, max_iter_num);
    KSPSetFromOptions(ksp);
    KSPSetOperators(ksp, mat, mat, 
      SAME_PRECONDITIONER);
    KSPSolve(ksp,rhs,x);


    PetscReal r_norm;
    KSPGetResidualNorm(ksp, &r_norm);
    KSPGetIterationNumber(ksp, &num_of_iter);

	//ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of processors = %d, rank = %d\n",size,rank);CHKERRQ(ierr);
	cout << "max_iter_num\t" << max_iter_num << endl;
	cout << "iter_error\t" << iter_error << endl;

	cout << "Matrix solver step " << num_of_iter << ", residual " << r_norm << ".\n";

	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp,&reason);
	if(reason < 0){
		cout << "PetSc did not solve the equations ! " << endl;
		exit(1);
	}

	PetscScalar *p;
    VecGetArray(x, &p);
    for(int i=0; i<L; i++) {
	b[i] = p[i];
    }
    VecRestoreArray(x, &p);

cout <<"Iterations for convergence="<< num_of_iter << " - Residual Norm = " << r_norm << endl;
		pi.done();

	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr); //instead of -ksp_view command line option 

    KSPDestroy(ksp);
    VecDestroy(rhs);
    VecDestroy(x);
}
//=============================================================
void PETScSparseMatrix::zeroYourself(void)
{
	const int L=numberOfEquations;
	//for(i=0;i<L;i++) for(j=0,limit=row[i].getNum();j<limit;j++) mat[i][j]=0;
	MatZeroEntries(mat);
}
//=============================================================
void PETScSparseMatrix::printMatrix(char *message)
{
	const int L=numberOfEquations;
	cout << message << "\nPETScSparseMatrix with " << numberOfEquations <<
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
void PETScSparseMatrix::printMatrixAddressFormat(char * message, ostream* out)
{ 
	//*out << message <<endl; 

	int numNonZeros=0;
	int i,j,limit;
	const int L=numberOfEquations;

  MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY); 
  MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
  //counting nonzeros
	/*
  MatAssemblyBegin(*(Mat*)mat, MAT_FINAL_ASSEMBLY); 
  MatAssemblyEnd(*(Mat*)mat, MAT_FINAL_ASSEMBLY);
  MatInfo info;
  MatGetInfo(*(Mat*)mat, MAT_LOCAL, &info);
  numNonZeros = (int)info.nz_used;
  */

  //remember when using symmetric storage in petsc, upper triagular is stored
  //whereas in alpha the non-zero location information is stored in lower triangle format
	for(i=0;i<L;i++) {
		limit=row[i].getNum();
		numNonZeros+=limit;
	}

	*out << numberOfEquations << " " << numNonZeros << endl;

#define FORMAT setw(30)<<setprecision(20)

	PetscScalar val;
	int ii,jj;

	for(i=0;i<L;i++) {
		for(j=0,limit=row[i].getNum();j<limit;j++){
			/*
			ii=i;
			jj=row[i][j];
			*/
			jj=i;
			ii=row[i][j];
			//flipping ii and jj to get actual non-zero location in upper triangle
			MatGetValues(mat, 1, &ii, 1, &jj, &val);
			//but writing out in the form of lower triangle to keep it consistent with other solves
			//*out << '\t' << ii << " " << jj << " " << FORMAT<< val;
			*out << '\t' << jj << " " << ii << " " << FORMAT<< val;
		}
		*out << '\n';
	}
}
#undef FORMAT 
//=============================================================
void PETScSparseMatrix::readMatrixAddressFormat(istream* instream)
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
		//this block if unsymmetric matrix
		/*
		MatSetValue(*(Mat*)mat, ii, jj, value, ADD_VALUES);
		if(ii != jj)
			MatSetValue(*(Mat*)mat, jj, ii, value, ADD_VALUES);
		//=======
		*/
	}

}
//=============================================================
void PETScSparseMatrix::vectorProduct(double *result, double *vector)
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
void PETScSparseMatrix::specifyNonZeroLocation(const int I, const int J)
{ 
	//row[I].add(J); 
	row[J].add(I); // petsc stores the upper triangle (the other solvers store the lower triangle)
					// and so the indices need to be switched.
}
//=============================================================
void PETScSparseMatrix::allocate(void)
{
	const int L=numberOfEquations;
    //PetscInitializeNoArguments();
	//PetscInitialize(0,&arguments,(char *)0,help);
	if( doesFileExist("petscOptions.txt")==false){
		ofstream options("petscOptions.txt");
		options << "-version" << endl;
		options << "-info" << endl;
		options << "-ksp_view" << endl;
		options << "-ksp_type cg" << endl;
		options << "-ksp_converged_reason" << endl;
		options.close();
	}
	PetscInitialize(0,0,"petscOptions.txt",0);

	int i;
	int *nnz;
	nnz=new int[L];
	FilledDiagonal=new bool[L];

	int maxnnz=0,maxrow=-1;
	for(i=0;i<L;i++) {
		nnz[i]=row[i].getNum();
		if(nnz[i] > maxnnz){
			maxnnz=nnz[i];
			maxrow=i;
		}
	}
	cout << "max nnz = " << maxnnz << " max row = " << maxrow << endl;
	for(i=0;i<L;i++) FilledDiagonal[i]=false;

    //MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, L, L, PETSC_DEFAULT, nnz, &mat);

	MatCreate(PETSC_COMM_WORLD,&mat);
	MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,L,L);
	MatSetType(mat,MATSBAIJ);  //symmetric
    MatSeqSBAIJSetPreallocation(mat, 1, nnz); 
	MatSetFromOptions(mat);

	MatSetOption(mat,MAT_NEW_NONZERO_ALLOCATION_ERR);
	delete [] nnz;
}
//=============================================================
#ifdef DEBUG_PETScSparseMatrix
ofstream Logfile;
void main(void)
{
	int i;
	ifstream is;
	ofstream os;
	is.open("GlobalStiffness.txt");
	int NumEquations;
	is >> NumEquations;
	PETScSparseMatrix A(NumEquations);
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
