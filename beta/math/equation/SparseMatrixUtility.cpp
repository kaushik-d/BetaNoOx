#include "stdafx.h"

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "utility/utility_defines.h"

//=============================================================
void writeCSRMatrixInASCII(char* filename, int* &ia, int* &ja, double* &a, int numberOfEquations, int ncoefs)
{
	const int L=numberOfEquations;
	int i;
	ofstream os(filename);
	SET_SCIENTIFIC(os);

	os << L << " " << ncoefs << endl;
	for(i=0; i<L+1; i++) os << ia[i] << " ";
	os << endl;
	for(i=0; i<ncoefs; i++) os << ja[i] << " ";
	os << endl;
	for(i=0; i<ncoefs; i++) os << DOUBLE_FORMAT << a[i] << endl;
	os << endl;
	os.close();
}
//=============================================================
void writeCSRMatrixInBinary(char* filename, int* &ia, int* &ja, double* &a, int numberOfEquations, int ncoefs)
{
	ofstream os;
	os.open(filename,ios::out|ios::binary);
	os.write(reinterpret_cast<char *>(&numberOfEquations),sizeof(numberOfEquations));
	os.write(reinterpret_cast<char *>(&ncoefs),sizeof(ncoefs));
	os.write(reinterpret_cast<char *>(ia),(numberOfEquations+1)*sizeof(int));
	os.write(reinterpret_cast<char *>(ja),(ncoefs)*sizeof(int));
	os.write(reinterpret_cast<char *>(a),(ncoefs)*sizeof(double));
	os.close();
}
//=============================================================
void readCSRMatrixInBinary(char* filename, int* &ia, int* &ja, double* &a, int &numberOfEquations, int &ncoefs)
{
	ifstream is;
	is.open(filename,ios::in|ios::binary);
	is.read(reinterpret_cast<char *>(&numberOfEquations),sizeof(numberOfEquations));
	is.read(reinterpret_cast<char *>(&ncoefs),sizeof(ncoefs));
	cout << "numberOfEquations :" << numberOfEquations << endl;
	cout << "ncoefs :" << ncoefs << endl;
	if(ia!=0)	delete [] ia;
	if(ja!=0)	delete [] ja;
	if(a!=0)	delete [] a;
	ia	= new int[numberOfEquations+1];
	ja	= new int[ncoefs];
	a	= new double[ncoefs];
	is.read(reinterpret_cast<char *>(ia),(numberOfEquations+1)*sizeof(int));
	is.read(reinterpret_cast<char *>(ja),(ncoefs)*sizeof(int));
	is.read(reinterpret_cast<char *>(a),(ncoefs)*sizeof(double));
	is.close();
}
//=============================================================
void readCSRMatrixInASCII(char* filename, int* &ia, int* &ja, double* &a, int &n)
{
	int i;
	int numeq,numnonzeros;
	ifstream is;
	is.open(filename);
	is >> numeq >> numnonzeros;
	n=numeq;
	if(ia!=0)	delete [] ia;
	if(ja!=0)	delete [] ja;
	if(a!=0)	delete [] a;
	ia	= new int[numeq+1];
	ja	= new int[numnonzeros];
	a	= new double[numnonzeros];
	for(i=0; i<numeq+1; i++) is >> ia[i];
	for(i=0; i<numnonzeros; i++) is >> ja[i];
	for(i=0; i<numnonzeros; i++) is >> a[i];
	is.close();
}
void readLoadVector(char* filename, double * b, int n)
{
	int i;
	ifstream is;
	is.open(filename);
	for(i=0; i<n; i++) is >> b[i];
	is.close();
}
void readMatrix(char* filename, int* &ia, int* &ja, double* &a, int &n)
{
	int numeq,numnonzeros;
	ifstream is;
	is.open(filename);
	is >> numeq >> numnonzeros;
	n=numeq;

	ia	= new int[numeq+1];
	ja	= new int[numnonzeros];
	a	= new double[numnonzeros];

	ia[numeq]=numnonzeros + 1;
	ia[0]=1;

	int i,ii,jj;
	int oldii=0;
	int ctr=0;
	double value;
	for(i=0;i<numnonzeros;i++){
		is >> ii >> jj >> value;
		a[ctr]=value;
		ja[ctr]=jj+1;
		if(ii != oldii){
			oldii = ii;
			ia[ii]=ctr+1;
		}
		ctr++;
	}
	is.close();
}
