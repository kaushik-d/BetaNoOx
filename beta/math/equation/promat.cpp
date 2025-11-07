#include "stdafx.h"

#include "promat.hpp"
#include "utility/excepts.hpp"
#include <iostream>
#include <cstdio>
#include <iomanip>
#include <cmath>

#define STATIC_DATA static

CreateErrorHandler(ProfileMatrix); 

ProfileMatrix::ProfileMatrix(const int numEqns) : LargeMatrix(numEqns)
{
	WhoClass("ProfileMatrix");
	WhoMethod("ProfileMatrix(const int)");
	int i;
	Assert(bandwidth = new int [numberOfEquations]);
	for(i=0;i<numberOfEquations;i++) bandwidth[i]=1;
	mat = 0;
	p = 0;
	maxBandwidth = 0;
}

ProfileMatrix::~ProfileMatrix(void)
{
	WhoMethod("~ProfileMatrix()");
	int i;
	if(mat!=0) {
		for(i=0;i<numberOfEquations;i++)
			delete [] mat[i];
		delete [] bandwidth;  //JV100404 fixed memory release
		delete [] mat;	
	}
	if(p) delete [] p;
}

ProfileMatrix::ProfileMatrix(const ProfileMatrix& a) : LargeMatrix(a) 
{
	WhoClass("ProfileMatrix");
	WhoMethod("ProfileMatrix(const ProfileMatrix& a)");
	int i,j;
	Assert(bandwidth = new int [numberOfEquations]);
	for(i=0;i<numberOfEquations;i++) bandwidth[i] = a.bandwidth[i];

	mat = 0;
	Assert(mat = new double *[numberOfEquations]);
	for(i=0;i<numberOfEquations;i++) {
		Assert(mat[i] = new double [bandwidth[i]]);
		for(j=0;j<bandwidth[i];j++) {
			mat[i][j] = a.mat[i][j];
		}
	}

	p = 0;
	// p is work vector
	Assert(p = new double [numberOfEquations]);
	for(i=0;i<numberOfEquations;i++) p[i] = a.p[i];

	maxBandwidth = a.maxBandwidth;

	factored = a.factored;
}

double& ProfileMatrix::operator() (const int i, const int j)
{
	register int tmp;
	if(i>=j) {
		tmp = j+bandwidth[i]-i-1;
		if(tmp<0 || tmp>i || i<0 || i >numberOfEquations)  {
			WhoMethod("operator()");
			cerr << "i,j,tmp = " << i << ", " <<j << ", " << tmp << endl;
			FatalError("Subscript out of range.\n");
		}
		return mat[i][tmp];
	}
	else { printf("Backwards (%d,%d)\n",i,j);
		tmp = i+bandwidth[j]-j-1;
		if(tmp<0 || tmp>j || j<0 || j>numberOfEquations) {
			WhoMethod("operator()");
			cerr << "j,i,tmp = " << j << ", " <<j << ", " << tmp << endl;
			FatalError("Subscript out of range.\n");
		}
		return mat[j][tmp];
	}
}

void ProfileMatrix::setaij(const int i, const int j, double val)
{
	operator()(i,j)=val;
}
double ProfileMatrix::getaij(const int i, const int j)
{
	return operator()(i,j);
}



void ProfileMatrix::zeroYourself(void)
{
	int i,j;
	WhoMethod("zeroYourself");
	if(mat==0) FatalError("Matrix not allocated yet. Can't zero.\n");
	for(i=0;i<numberOfEquations;i++)
		for(j=0;j<bandwidth[i];j++)
			mat[i][j]=0;
	factored = false;
}

void ProfileMatrix::solve(double *vector,double *dummy)
{
	double *t,*rowi,*rowj;
	double c;
	int i,np1,ii,j,iqq,imci;
	int firstColRowI,firstColRowJ,firstOverlapColumnT,ks,overlapLength,
		firstOverlapPointerRowJ;

	WhoMethod("Solve");

	if(factored) goto L160;

	// Check for positive definite
	CheckAndFixZeroDiagonal(); //can diagonals have zero or negative numbers if the matrix is factorized?

	cout << "--     Factoring\n";
	cout << "MaxBandwidth = " << maxBandwidth << '\n';
	t = new double [maxBandwidth];
	for(i=0;i<numberOfEquations;i++) {
		rowi = mat[i];
// Calculate First Column in Row i
		firstColRowI = i-bandwidth[i]+1;
// Store row i in t[]
		for(iqq=0;iqq<bandwidth[i];iqq++) t[iqq]=rowi[iqq]; 
// Set Diagonal Term on Row i to -1
		rowi[bandwidth[i]-1]=-1;
// Loop from left end of row to term next to row diagonal.
//	Need to from inner products of rows.  All products are zero except
//		where rows overlap. 
		for(j=firstColRowI;j<i;j++) {
			rowj = mat[j];
			firstColRowJ=j-bandwidth[j]+1;
//		Find largest starting column number.
			if(firstColRowI<firstColRowJ) ks = firstColRowJ;
			else ks = firstColRowI;
//		firstOverlapColumnT is location in vector t[] where overlap begins
//			between rows i & j.
			firstOverlapColumnT = ks-firstColRowI;
			overlapLength=j-ks+1;
//		firstOverlapPointerRowJ is pointer to location in row j where 
//			overlap starts.
			firstOverlapPointerRowJ=ks-firstColRowJ;
			c = dotProduct(overlapLength,&t[firstOverlapColumnT],
					&rowj[firstOverlapPointerRowJ]);
			firstOverlapPointerRowJ=j-firstColRowI;
			t[firstOverlapPointerRowJ]=-c;
			overlapLength=j-firstColRowI;
			rowi[overlapLength]=t[firstOverlapPointerRowJ]*p[j];
		}			
// Now Repeat Previous Block for j=i
		overlapLength = i-firstColRowI+1;
/*//JV072009 commented the check...
		// xtang 03092003 checking why there is a zero diagonal: hang node- when a node of 
		// in an element neither being constrained nor connected to any other element.
		double aaa=dotProduct(overlapLength,t,rowi);
		if (fabs(aaa) <1e-10) 
			cout << " zero diagonal, Equ. "<<i<<endl;
		p[i] = -1/aaa;
*/
		p[i] = -1/dotProduct(overlapLength,t,rowi);

		if(p[i] <= 0)
			cout << "Negative Diagonal, Eqn. " << i << ": " << 1/p[i] << endl;
	}
	delete t;
	factored=true;
L160:
	cout << "--     Substituting\n";
// Forward Substitution
	for(i=0;i<numberOfEquations;i++) {
		firstColRowI=i-bandwidth[i]+1;
//		for(j=0;j<1;j++)
//		{
			c = dotProduct(bandwidth[i],mat[i],&vector[firstColRowI]);
			vector[i] = -c;
//		} 
	}
// Diagonal
//	for(j=0;j<1;j++)
		for(iqq=0;iqq<numberOfEquations;iqq++)
			vector[iqq]*=p[iqq];

// BackwardSubstitution
	np1=numberOfEquations;
	for(ii=1;ii<numberOfEquations;ii++) {
		i=np1-ii;
		firstColRowI=i-bandwidth[i]+1;
		if(firstColRowI >= i) continue;
		imci=i-firstColRowI;
//		for(j=0;j<1;j++)
			for(iqq=0;iqq<imci;iqq++)
				vector[firstColRowI+iqq]-=mat[i][iqq]*vector[i];
	} 
}

inline double ProfileMatrix::dotProduct(int length, double* v1, double *v2)
{
	STATIC_DATA int i;
	register double val;
	val=0;
	for(i=0;i<length;i++)
		val+=v1[i]*v2[i];
	return val;			
}

void ProfileMatrix::printMatrix(char *s)
{
	int i,j;
	double sum;
	WhoMethod("printMatrix");
	cout << s<< '\n';
	for(i=0;i<numberOfEquations;i++) {
//		for(j=0;j<i-bandwidth[i]+1;j++)
//			cout << "   0  ";
		cout << "row " << i << " | ";
		sum=0;
		for(j=0;j<bandwidth[i];j++) {
			sum+=mat[i][j];
			cout << mat[i][j] << "  ";
		}
		cout << "= " << sum << '\n';
	}
}
#define FORMAT setw(13)<<setprecision(4)

void ProfileMatrix::printMatrix(char *s, ostream* out)
{
	int i,j;
	double sum;
	WhoMethod("printMatrix");
	(*out) << s<< '\n';
	int columnwidth=0,number=numberOfEquations;
	while (number > 0 ) {
		number=number/10;
		columnwidth++;
	}
	for(i=0;i<numberOfEquations;i++) {
//		for(j=0;j<i-bandwidth[i]+1;j++)
//			cout << "   0  ";

		(*out) << "row " << setw(columnwidth)<<setprecision(4) << i << " | ";
		sum=0;
		for(int jj=bandwidth[i]; jj < (i+1); jj++)	//JV072203 this line was added so that 
			(*out) << FORMAT << 0 << "  ";  //the printed out matrix is more readable											
		for(j=0;j<bandwidth[i];j++) {
			sum+=mat[i][j];
			(*out) << FORMAT << mat[i][j] << "  ";
		}
		(*out) << endl;
		//(*out) << "= " << sum << '\n';   //JV072203 - taken off for more readability 
											// uncomment if u need to print the sum !!!
	}
}
#undef FORMAT

void ProfileMatrix::printMatrixAddressFormat(char *s, ostream* out)
{
	printMatrix(s, out);
}

void ProfileMatrix::specifyNonZeroLocation(const int i, const int j)
{
	// For row i check to see if required bandwidth is large enough
	//		if not increase it to required size.
	if(i-j+1>bandwidth[i]) bandwidth[i]=i-j+1;
	//in essence, this function specifies the distance from the (i,j) to (i,i)
	//rather than specifying the 'non-zero location'
	//and memory is allocated for the whole 'distance' rather than 
	//just storing info about that non-zero point
	//this is why profile takes more memory space.
}

bool ProfileMatrix::checkNonZeroLocation(const int i, const int j)
{
	WhoClass("ProfileMatrix");
	WhoMethod("checkNonZeroLocation(const int i, const int j)");
	if(i >= j){
		if(i-j+1>bandwidth[i]) 
			return false;
		else 
			return true;
	}else{
		if(j-i+1>bandwidth[j]) 
			return false;
		else 
			return true;
	}
}

void ProfileMatrix::allocate(void)
{
	WhoMethod("allocate(void)");
	int i,totalMemory=0;
	maxBandwidth = 0;
	for(i=0;i<numberOfEquations;i++) {
		if(bandwidth[i] > maxBandwidth) maxBandwidth = bandwidth[i];
		// Calculate memory requirements
		totalMemory += bandwidth[i];
	}
	cout << "ProfileMatrix::totalMemory = " << totalMemory << " words." << endl;
	Assert(mat = new double *[numberOfEquations]);
	for(i=0;i<numberOfEquations;i++) Assert(mat[i] = new double [bandwidth[i]]);
	// p is work vector
	Assert(p = new double [numberOfEquations]);
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
