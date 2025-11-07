#include "stdafx.h"
#include<cmath>
#include<string>
using namespace std;
#include "utility/utility_defines.h"
//==========================================================
void copyVector(const double * copyFrom, double * copyTo, int length)
{
    memcpy(copyTo, copyFrom, length*sizeof(double));
}

//==========================================================
double dotProduct( const double * a, const double * b, const int &length)
{
	double sum =0.;
	for(int i=0; i<length; i++)
	{ sum += a[i] * b[i]; }
	return(sum);
}
//==========================================================
void normalizeVector( double * a, const int &length)
{
	double Val = sqrt(dotProduct(a,a,length));
	for(int i=0; i<length; i++)
	{ a[i]/=Val; }
}
//==========================================================
void sumVector( double *c, double * a, double * b, int length)
{
	for(int i=0; i<length; i++)
	{ c[i] = a[i] + b[i] ; }
}
//==========================================================
void initializeVector( double * a, double value, int length)
{
	for(int i=0; i<length; i++)
	{ a[i] = value ; }
}
//==========================================================
void ConvertNumericalZeroToZero(double *a, int length, double numericalZero)
{
	for(int i=0; i<length; i++){ 
		if( fabs(a[i]) < numericalZero)
			a[i] = 0.0; 
	}
}
//==========================================================
int ConvertNegativeNumbersToZero(double *a, int length)
{
	int index = -1;
	for(int i=0; i<length; i++){ 
		if( a[i] < 0.0){
			index=i;
			a[i] = 0.0; 
		}
	}
	return index;
}
//==========================================================
void WriteDoubleVectorToFile(string filename, double *list, int L)
{
	ofstream outfile;
	outfile.open(filename.c_str());
	outfile.setf(ios::scientific, ios::floatfield);
	outfile << L << endl;
	for(int i=0;i<L;i++) outfile << DOUBLE_FORMAT << list[i] << endl;
//	for(int i=0;i<L;i++) outfile <<  setw(23)<<setprecision(12) << list[i] << endl;
	outfile.close();

}
//==========================================================
bool ReadDoubleVectorFromFile(string filename, double *&list, int L)
{
	ifstream infile;
	infile.open(filename.c_str());
	if(!infile.good()) return false;
	for(int i=0;i<L;i++) infile >>list[i];
	infile.close();
	return true;
}
//==========================================================
void getDoubleMinMax(double *x,int length, double &min, double &max)
{
	min = x[0];
	max = x[0];
	for (int i=1;i<length;i++){
		if(x[i] > max)
			max = x[i];

		if(x[i] < min)
			min = x[i];
		
	}

}
//=========================================================================