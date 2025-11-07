#include "stdafx.h"

#include "matrix.hpp"
#include <iomanip>

#define FORMAT setw(13)<<setprecision(4)

#define STATIC_DATA static
 
//===================== CONSTRUCTORS ================
Matrix::Matrix(const Matrix& Arg)
{
	// n=#rows and m=#columns
	isSymmetric=Arg.isSymmetric;
	allocate(Arg.numRows, Arg.numCols);
	(*this)=0.0;
	copy(Arg);
}

Matrix::Matrix(Matrix& Arg, int nn, int mm)
{
	int i,j;

	numRows=nn;
	numCols=mm;
    allocate(nn,mm);
	for(i=0;i<numRows;i++) {
		skip[i] = Arg.skip[i];
		for(j=0;j<numCols;j++)
			(*this)(i,j)=Arg(i,j);
	}
}


Matrix::Matrix(int nn,int mm)
{
	isSymmetric=true;
	allocate(nn,mm);
	(*this)=0.0;
}
//===================== DESTRUCTORS ================
Matrix::~Matrix(void)
{
	delete[] skip;
    delete [] data;
}	
//===================== OVERLOADED OPERATORS ================
Matrix Matrix::operator* (const Matrix& Arg) const
{
	int i,j,k;
	if(numCols!=Arg.numRows) {
		cerr << "\nMatrix::Error in operator * - (numCols!=Arg.numRows)\n";
		cerr << "This operator can be currently used only for square matrices \n";
		exit(1);
	}
	Matrix Temp(numRows,Arg.numCols);
	for(i=0;i<numRows;i++)
		for(k=0;k<Arg.numCols;k++) {
			Temp(i,k) = 0;
			for(j=0;j<numCols;j++)
				Temp(i,k) = Temp(i,k)+(*this)(i,j)*Arg(j,k);
		}
	return Temp;
}


ostream &operator<<(ostream &o, const Matrix &Arg)
{
	int i,j;

	o << '[' << Arg.numRows << " x " << Arg.numCols << ']' << '\n';
	for(i=0;i<Arg.numRows;i++) {
		o << "| ";
		for(j=0;j<Arg.numCols;j++)
			{ o <<FORMAT<< Arg(i,j) ;}
		o << " |";
		if(Arg.skip[i]) o << " X";
		//o << '\n';
                o << endl;
	}             
	return o;
}

Matrix Matrix::operator= (const Matrix& Other_Matrix)
{
	if(numRows!=Other_Matrix.numRows || numCols!=Other_Matrix.numCols) {
		resize(Other_Matrix.numRows, Other_Matrix.numCols);
	}
	copy(Other_Matrix);
	return *this;
}

const double &Matrix::operator= (const double &val)
{
	int i;

	for(i=0;i<numRows;i++){
		skip[i]=0;
	}

    for(i=0;i<numRows*numCols;++i)
        data[i]=val;
    
	return val;
}

Matrix Matrix::operator+ (const Matrix& Arg) const
{
	int i,j;
	if(numRows!=Arg.numRows || numCols!=Arg.numCols) {
		cerr << "\nMatrixClass::Invalid Matrix Operation in operator +";
		exit(1);
	}
	Matrix Temp(numRows,numCols);
	for(i=0;i<numRows;i++)
		for(j=0;j<numCols;j++)
			Temp(i,j)=(*this)(i,j)+Arg(i,j);
	return Temp;
}

Matrix Matrix::operator- (const Matrix& Arg) const 
{
	int i,j;
	if(numRows!=Arg.numRows || numCols!=Arg.numCols) {
		cerr << "\nMatrixClass::Invalid Matrix Operation in operator -";
		exit(1);
	}
	Matrix Temp(numRows,numCols);
	for(i=0;i<numRows;i++)
		for(j=0;j<numCols;j++)
			Temp(i,j)=(*this)(i,j)-Arg(i,j);
	return Temp;
}

void Matrix::operator+= (const Matrix& Arg)
{
	int i,j;
	if(numRows!=Arg.numRows || numCols!=Arg.numCols) {
		cerr << "\nMatrixClass::Invalid Matrix Operation in operator+=";
		exit(1);
	}
	for(i=0;i<numRows;i++)
		for(j=0;j<numCols;j++)
			(*this)(i,j)+=Arg(i,j);
} 

void Matrix::operator-= (const Matrix& Arg)
{
	int i,j;
	if(numRows!=Arg.numRows || numCols!=Arg.numCols) {
		cerr << "\nMatrixClass::Invalid Matrix Operation in operator-=";
		exit(1);
	}
	for(i=0;i<numRows;i++)
		for(j=0;j<numCols;j++)
			(*this)(i,j)-=Arg(i,j);
}
//=================================================
Matrix Matrix::operator* (const double& Arg) const //multiply each term
{
	int i,j;
	Matrix Temp(numRows,numCols);

	for(i=0;i<numRows;i++)
		for(j=0;j<numCols;j++)
			Temp(i,j) = (*this)(i,j) * Arg;
	return Temp;
}
//=================================================
Matrix Matrix::multiplyWithMatrix(const Matrix& Arg, int MatrixProperty) //JV112005
// when using MATRIX_SYMMETRIC, the matrix has to be upper triangle stored
{
	int i,j,k;
	if(numCols!=Arg.numRows) {
		cerr << "\nMatrix::Error in operator * - (numCols!=Arg.numRows)\n";
		exit(1);
	}
	Matrix Temp(numRows,Arg.numCols);
	for(i=0;i<numRows;i++)
		for(k=0;k<Arg.numCols;k++) {
			Temp(i,k) = 0;
			for(j=0;j<numCols;j++){
				if(MatrixProperty==MATRIX_GENERAL){
					Temp(i,k) += (*this)(i,j)*Arg(j,k);
				}else if (MatrixProperty==MATRIX_SYMMETRIC){
					if(j<i) //i.e. if it is lower triangle , switch indices
						Temp(i,k) += (*this)(j,i)*Arg(j,k);
					else
						Temp(i,k) += (*this)(i,j)*Arg(j,k);
				}
			}
		}
	return Temp;
}
//=================================================
void Matrix::multiplyWithVector(const double* v, const int num, double *&F, int MatrixProperty) 
//multiply with vector v. output is F
// when using MATRIX_SYMMETRIC, the matrix has to be upper triangle stored
{
	int i,j;
	double sum;
	if(numCols != num) {
		cout << "wrong parameters in Matrix Matrix::operator* (const double* v, const int num, double *F)" << endl;
		exit(1);
	}
		
	for(i=0;i<numRows;i++){
		sum=0.0;
		for(j=0;j<numCols;j++){
			if(MatrixProperty==MATRIX_GENERAL){
				sum += (*this)(i,j) * v[j];
			}else if (MatrixProperty==MATRIX_SYMMETRIC){
				if(j<i) //i.e. if it is lower triangle , switch indices
					sum += (*this)(j,i) * v[j];
				else
					sum += (*this)(i,j) * v[j];
			}
		}
		F[i]=sum;
	}
}
//=================================================
void Matrix::convertSymmetricToGeneral()
{
	int i,j;
	for(i=0;i<numRows;i++){
		for(j=0;j<=i;j++){
			(*this)(i,j) = (*this)(j,i);
		}
	}
}
//=================================================
void Matrix::resize(int nn,int mm)
//Resizes and destroys contents
{
    //delete[] skip;
    //delete [] data;
    this->~Matrix();
    allocate(nn,mm);
}
//=================================================
void Matrix::allocate(int nn, int mm)
{
	numRows=nn;
	numCols=mm;
    data = new Coefficient [numRows*numCols];
	skip = new int [numRows];
}
//=================================================
void Matrix::copy(const Matrix& Other_Matrix)
{

	if(numRows!=Other_Matrix.numRows || numCols!=Other_Matrix.numCols) {
		cerr << "In = for Matrix\n";
		cerr << "numRows = " << numRows << "\nnumCols = " << numCols << "\nOther_Matrix.numRows = " << 
				Other_Matrix.numRows << "\nOther_Matrix.numCols = " << Other_Matrix.numCols << '\n';
		cerr << "\nInvalid Matrix Operation in operator= \n" << endl;
		exit(1);
	}
    //Copy skip
    memcpy(skip,Other_Matrix.skip,numRows*sizeof(int));
    //copy data
    memcpy(data,Other_Matrix.data,numRows*numCols*sizeof(double));
}
//=================================================
//========== DECOMPOSITION AND SOLVING ============
//*****************
// The following is adopted from Numerical Recipes
//*****************
STATIC_DATA const double TINY=1.0e-20;

void Matrix::ludcmp(int *indx,double *d)
{
	int i,imax,j,k,last_eqn;
	double big,dum,sum,temp;
	double *vv;

	if(numRows!=numCols) {
		cerr << "Error: Matrix::ludcmp - (numRows!=numCols)\n";
		exit(1);
	}
	for(i=0;i<numRows;i++) if (!skip[i]) last_eqn=i;
	vv = new double [numRows];
	if(vv==NULL) {
		cerr << "\n\nOut of memory in ludcmp!! Exiting.\n";
		exit(1);
	}
	*d=1.0;
	for (i=0;i<numRows;i++) {
		if(skip[i])continue;
		big=0.0;
		for (j=0;j<numRows;j++) {
			if(skip[j])continue;
			if ((temp=fabs((*this)(i,j))) > big) big=temp;
		}
		if (big == 0.0) {
			cerr << "Error: Matrix::ludcmp - Singular matrix!! Exiting.\n";
			exit(1);
		}
		vv[i]=1.0/big;
	}
	for (j=0;j<numRows;j++) {
		if(skip[j]) continue;
		for (i=0;i<j;i++) {
			if(skip[i]) continue;
			sum=(*this)(i,j);
			for (k=0;k<i;k++) {
				if(skip[k]) continue;
				sum -= (*this)(i,k)*(*this)(k,j);
			}
			(*this)(i,j)=sum;
		}
		big=0.0;
		for (i=j;i<numRows;i++) {
			if(skip[i])continue;
			sum=(*this)(i,j);
			for (k=0;k<j;k++) {
				if(skip[k]) continue;
				sum -= (*this)(i,k)*(*this)(k,j);
			}
			(*this)(i,j)=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<numRows;k++) {
				if(skip[k]) continue;
				dum=(*this)(imax,k);
				(*this)(imax,k)=(*this)(j,k);
				(*this)(j,k)=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if ((*this)(j,j) == 0.0) (*this)(j,j)=TINY;
		if (j != last_eqn /*numRows-1*/) {
			dum=1.0/((*this)(j,j));
			for (i=j+1;i<numRows;i++) {
				if(skip[i])continue;
				(*this)(i,j) *= dum;
			}
		}
	}
	delete vv;
}

Matrix Matrix::invert(void)
{
	double d;
	int i,j,*indx;
                 
	for(i=0;i<numRows;i++) skip[i]=0; // Don't skip any equations when inverting.
	indx = new int [numRows];
	ludcmp(indx,&d);
	Matrix B(numRows,numRows);
	for(i=0;i<numRows;i++)
		B[i][i]=1;
	for(i=0;i<numRows;i++)
		lubksb(indx,B[i]);

	for(i=0;i<numRows;i++)
		for(j=0;j<numCols;j++)
			(*this)(i,j) = B(j,i);

	delete indx;
	return *this;
}


void Matrix::SolverSkip(int i)
{
	skip[i]=1;
}

void Matrix::SolverDontSkip(int i)
{
	skip[i]=0;
}

void Matrix::Solve(double *vector)
{
	double d;
	int *indx;

	indx = new int [numRows];
	ludcmp(indx,&d);
	lubksb(indx,vector);
	delete indx;
}
void Matrix::lubksb(int *indx,double *b)
{
	int i,ii=-1,ip,j;
	double sum; 

// Augment b
	for(i=0;i<numRows;i++) {
		if(skip[i]) {
			for(j=0;j<numRows;j++) {
				if(skip[j])continue;
					b[j]-=(*this)(j,i)*b[i];
			}
		}
	}

	for (i=0;i<numRows;i++) {
		if(skip[i])continue;
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii>=0)
			for (j=ii;j<i;j++) {
				if(skip[j])continue;
				sum -= (*this)(i,j)*b[j];
			}
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=numRows-1;i>=0;i--) {
		if(skip[i])continue;
		sum=b[i];
		for (j=i+1;j<numRows;j++) {
			if(skip[j])continue;
			sum -= (*this)(i,j)*b[j];
		}
		b[i]=sum/(*this)(i,i);
	}
}
//=================================================
//================OUTPUT===========================
//=================================================
void Matrix::print(char *s) const
{
	cout << s << ' ' << (*this);
}
//========================================================
void Matrix::print(char *s, ostream *out) const
{
	(*out) << s << ' ' << (*this);
}

//========================================================
void Matrix::print(char *s, int m, int n) const
{
 cout << s << ' ' <<endl;
 
 
 	int i,j;

// Fix	cout << '[' << Arg.numRows << " x " << Arg.numCols << ']' << '\n';
	for(i=0;i<m;i++) {
		cout << "| ";
		for(j=0;j<n;j++)
			{ cout <<FORMAT<< (*this)(i,j) ;}
		cout << " |";
//Fix		if(Arg.skip[i]) cout << " X";
		//o << '\n';
                cout << endl;
	}             
}

//========================================================
void Matrix::print(char *s, int m, int n, ostream * out) const
{
 (*out) << s << ' ' <<endl;
 
 
 	int i,j;

// Fix	cout << '[' << Arg.numRows << " x " << Arg.numCols << ']' << '\n';
	for(i=0;i<m;i++) {
		(*out) << "| ";
		for(j=0;j<n;j++)
			{ (*out) <<FORMAT<< (*this)(i,j) ;}
		(*out) << " |";
//Fix		if(Arg.skip[i]) (*out) << " X";
		//o << '\n';
                (*out) << endl;
	}             
}
//========================================================
void Matrix::printWithRowColumnNumbering(char *s, int m, int n, int* rowNumber, int* columnNumber, ostream * out) const
{
	(*out) << s << ' ' <<endl;
	(*out) << setiosflags(ios::scientific) <<endl;
 
 	int i,j;

	(*out) << "       ";
	for(i=0;i<n;i++)
		(*out) << setw(7) << columnNumber[i] << "      ";
	(*out) << endl;

	// Fix	cout << '[' << Arg.numRows << " x " << Arg.numCols << ']' << '\n';
	for(i=0;i<m;i++) {
		(*out) << setw(5) << rowNumber[i] << "| ";
		for(j=0;j<n;j++)
			{ (*out) <<FORMAT<< (*this)(i,j) ;}
		(*out) << " |";
//Fix		if(Arg.skip[i]) (*out) << " X";
		//o << '\n';
                (*out) << endl;
	}             
}
//========================================================
void Matrix::printInMaple(char *s, int mm, int nn, ostream * out) const
{
 //(*out) << s << ' ' <<endl;
 
 (*out) << s << ":=Matrix(" <<endl;
 
 	int i,j;
		(*out) << "[ ";

	for(i=0;i<mm;i++) {
		(*out) << "[ ";
		for(j=0;j<nn;j++){
			(*out) <<FORMAT<< (*this)(i,j);
			if (j!=nn-1)
				(*out) <<", ";
		}
		(*out) << " ]";
		if (i==mm-1)
			(*out) <<"]);";
		else
			(*out) <<", ";
		(*out) << endl;
	}             
}
//=================================================
//=============MISCELLANEOUS=======================
Matrix Matrix::transpose(void)
{
	int i,j;

	if(numRows == numCols) {
		//quicker in-place transposition for square matrices
		double temp;
		for(i=1;i<numRows;i++)
			for(j=0;j<numCols-i;j++) {
				temp = (*this)(i,j);
				(*this)(i,j) = (*this)(j,i);
				(*this)(j,i) = temp;
			}
	}else{
		//non-square

		//Although it shouldn't come up, transposing a
		//non-squiare matrix will mess up skip[].  Perform
		//a check and if any element of skip is non-zero,
		//print a warning to cout.  Then resize skip appropriately
		bool skipFlag = false;
		for(i=0;i<numRows;++i){
			if(skip[i]) skipFlag = true;
		}
		if(skipFlag) cout<<"WARNING - Transposing a non-square matrix deletes skipped equation data!" << endl;
		delete [] skip;
		skip = new int [numCols];
		for(i=0;i<numCols;++i){
			skip[i] = 0;
		}
		//Switch the dimensions
		int temp = numCols;
		numCols = numRows;
		numRows = temp;

		//Transpose data
		double * oldData = data;
		data = new double [numRows*numCols];
		for(i=0;i<numRows;++i){
			for(j=0;j<numCols;++j){
				data[i*numCols+j] = oldData[j*numRows+i];
			}
		}
		delete [] oldData;
	}
	return *this;
}

Matrix Matrix::transposeSelf(void) //DG_Feb2004
{
	int i,j;

	if(numRows == numCols) {
		//quicker in-place transposition for square matrices
		double temp;
		for(i=1;i<numRows;i++)
			for(j=0;j<numCols-i;j++) {
				temp = (*this)(i,j);
				(*this)(i,j) = (*this)(j,i);
				(*this)(j,i) = temp;
			}
	}else{
		//non-square

		//Although it shouldn't come up, transposing a
		//non-squiare matrix will mess up skip[].  Perform
		//a check and if any element of skip is non-zero,
		//print a warning to cout.  Then resize skip appropriately
		bool skipFlag = false;
		for(i=0;i<numRows;++i){
			if(skip[i]) skipFlag = true;
		}
		if(skipFlag) cout<<"WARNING - Transposing a non-square matrix deletes skipped equation data!" << endl;
		delete [] skip;
		skip = new int [numCols];
		for(i=0;i<numCols;++i){
			skip[i] = 0;
		}

		//Switch the dimensions
		int temp = numCols;
		numCols = numRows;
		numRows = temp;

		//Transpose data
		double * oldData = data;
		data = new double [numRows*numCols];
		for(i=0;i<numRows;++i){
			for(j=0;j<numCols;++j){
				data[i*numCols+j] = oldData[j*numRows+i];
			}
		}
		delete [] oldData;
	}
	return *this;
}


// xtang 03202003: corrected
Matrix Matrix::transposeCopy(void) //DG_Feb2004
{
	int i,j;
	Matrix B(numCols,numRows);

	for(i=0;i<numRows;i++)
		for(j=0;j<numCols;j++) {
			B(j,i) = (*this)(i,j);
		}
//	(*this) = B;
//	return *this;
	return B;
}

Matrix* Matrix::transposeCopyNew(void) //JV102405
{
	int i,j;
	Matrix *B= new Matrix(numCols,numRows);

	for(i=0;i<numRows;i++)
		for(j=0;j<numCols;j++) {
			(*B)(j,i) = (*this)(i,j);
		}
	return B;
}
//========================================================
#ifdef DEBUG_MATRIX_CLASS
void main(void)
{ 
	double AA[4] = {5,5,5};
	Matrix A(2,2);
	Matrix B(2,2);
	Matrix C(2,2);
	Matrix E(3,3);
	Matrix F(4,4);
	Matrix *D;

	A[0][0] = 1;
	A[0][1] = 2;
	A[1][0] = 3;
	A[1][1] = 1;
	B[0][0] = 2;
	B[0][1] = 4;
	B[1][0] = 7;
	B[1][1] = -3;
	cout << "\n-------------\n";
	cout << "A = {{1,2},{3,1}}" << A;
	cout << "B = {{2,4},{7,-3}}" << B;

	C = A;
	cout << "C = A " << C;
	C = A + B;
	(A+B).print("(A+B) {{3,6},{10,-2}}");
	cout << "C = A + B = " << C;
	C = A * B;
	cout << "C = A * B = {{16,-2},{13,9}}" << C;
	C.transposeCopy();  //DG_Feb2004 // xtang 03202003:  works only for 2x2 matrix, but not correct for others
	cout << "C.transpose" << C;
	A.Solve(AA);
	cout << "A.Solve {{3,1},{.3333,1.6667}}" << A;
	cout << '{' << AA[0] << ' ' << AA[1] << "}\n";
	cout << "B = {{2,4},{7,-3}}" << B;
	B.invert();
	cout << "B.Invert = {{0.088235,0.11765},{0.20588,-.058824}}" << B;
	D=&C;
	(*D)[0][0]=50;
	cout << "D[0][0]=50" << (*D);
	E[0][0] = 1;
	E[0][1] = 2;
	E[1][0] = 3;
	E[1][1] = 1;
	E[0][2] = 3;
	E[2][0] = 600;
	E.SolverSkip(2);
	cout << "Solver will skip 2\n";
	cout << "E = " << E;
	AA[0]=20;
	AA[1]=5;
	AA[2]=5;
	E.Solve(AA);
	cout << '{' << AA[0] << ' ' << AA[1] << ' ' << AA[2] <<"}\n";
	cout << "E.Solve = " << E;
	F[0][0] = 1; F[0][1] = -2; F[0][2] = 3; F[0][3] = 4;
	F[1][0] = 3; F[1][1] =-3; F[1][2] = 4; F[1][3] =.5;
	F[2][0] =-2; F[2][1] = 2; F[2][2] = 12; F[2][3] = 1;
	F[3][0] = -3; F[3][1] =-3; F[3][2] = 6; F[3][3] =.5;
	AA[0] = 6;    AA[0] = 1; F.SolverSkip(0);
	AA[1] = 4.5;  AA[1] = 1; F.SolverSkip(1);
	AA[2] = 13;   AA[2] = 1; F.SolverSkip(2);
	AA[3] = .5;   //AA[3] = 1; F.SolverSkip(3);
	F.Solve(AA);
	cout << '{' << AA[0] << ' ' << AA[1] <<' '<< AA[2] <<' '<< AA[3] << "}\n";
	F=0;
	cout << "F = 0 "<< F << endl;
}
#endif
