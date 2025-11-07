#include "stdafx.h"

#include <cmath>
#define DEBUG_Level0
#include "tensor.hpp"

////////////////////////////////////////////////////////////////////////////
//
//  T e n s o r
//
CreateErrorHandler(Tensor);
const TensorValue Tensor::zeros[9] = {0};
Tensor::Tensor()
{
	WhoClass("Tensor");
	WhoMethod("Tensor()");
	TensorValue *p;
	T = (new TensorValue *[3]) - 1;
	p = new TensorValue[9]; // This insures contiguous memory.
	T[1] = p-1;
	T[2] = p+2;
	T[3] = p+5;
	zero();
}
Tensor::~Tensor()
{
	WhoMethod("~Tensor()");
	delete [] (T[1]+1);
	delete [] (T+1);
}
void Tensor::write(FILE *fp)
{
	WhoMethod("write(FILE *)");
	if(fwrite(T[1]+1,sizeof(TensorValue),9,fp)!=9) {
		WhoMethod("write(FILE *)");
		FatalError("Cannot write tensor.");
	}
}
void Tensor::read(FILE *fp)
{
	WhoMethod("read(FILE *)");
	if(fread(T[1]+1,sizeof(TensorValue),9,fp)!=9) {
		WhoMethod("read(FILE *)");
		FatalError("Cannot read tensor.");
	}
}
	
ostream& operator <<(ostream& o, Tensor& a)
{
	o << "|  " << a(1,1) << "  " << a(1,2) << "  " << a(1,3) << "|  " <<  endl;
	o << "|  " << a(2,1) << "  " << a(2,2) << "  " << a(2,3) << "|  " <<  endl;
	o << "|  " << a(3,1) << "  " << a(3,2) << "  " << a(3,3) << "|  " <<  endl;
	return o;
}
////////////////////////////////////////////////////////////////////////////
//
//  R o t a t i o n  T e n s o r
//
RotationTensor::RotationTensor()
{
	WhoClass("RotationTensor");
}
RotationTensor::RotationTensor(const double &angle, const int axis)
{
	double theta,ct,st;
	WhoClass("RotationTensor");

	// Change angle to radians;
	theta=angle*M_PI/180.0;

	// Calculate sin and cos of theta
	ct=cos(theta); st=sin(theta);

	// Zero rotation matrix
	zero();

	// Setup rotation matrix
	switch (axis) {
		case 1:
			(*this)(1,1)=1.;
			(*this)(2,2)=ct;
			(*this)(2,3)=st;
			(*this)(3,2)=-st;
			(*this)(3,3)=ct;
			break;
		case 2:
			(*this)(2,2)=1.;
			(*this)(1,1)=ct;
			(*this)(1,3)=-st;
			(*this)(3,1)=st;
			(*this)(3,3)=ct;
			break;
		case 3:
			(*this)(3,3)=1.;
			(*this)(1,1)=ct;
			(*this)(1,2)=st;
			(*this)(2,1)=-st;
			(*this)(2,2)=ct;
			break;
		default:
			cerr << "axis = " << axis << endl;
			FatalError("Bad rotation axis.");
	}
}

RotationTensor::RotationTensor(const float &T1,const float &T2, const float &T3)
{
   TensorValue s1=sin(T1),c1=cos(T1),s2=sin(T2),c2=cos(T2),s3=sin(T3),c3=cos(T3);
   T[1][1]=c2*c1;          T[1][2]=c2*s1;          T[1][3]=-s2;
   T[2][1]=s3*s2*c1-c3*s1; T[2][2]=s3*s2*s1+c3*c1; T[2][3]=s3*c2;
   T[3][1]=c3*s2*c1+s3*s1; T[3][2]=c3*s2*s1-s3*c1; T[3][3]=c3*c2;
}

void RotationTensor::determineXYZRotations(double &T1, double &T2, double &T3)
{
	// This method ....but it works... Sorta...
	double c2;
	TensorValue **lt = T;

   TensorValue c1c3=-lt[2][1]*lt[3][3]*lt[1][1]/
			(lt[2][3]*lt[1][1]*lt[1][3]-lt[2][1]*lt[3][3]);
	if(c1c3<0) {
		lt[1][1]=-lt[1][1]; 
		lt[1][2]=-lt[1][2]; 
		lt[1][3]=-lt[1][3];
	}
	if(lt[1][1]>0) { // c2 is positive
		T2 = asin(-lt[1][3]);
		c2 = cos(T2);
		if(lt[1][2]>0) T1 = asin(lt[1][2]/c2);
		else T1 = -asin(lt[1][2]/c2);
		if(lt[2][3]>0) T3 = asin(lt[2][3]/c2);
		else T3 = -asin(lt[2][3]/c2);
	} else { // c2 is negative
		T2 = M_PI-asin(-lt[1][3]);
		c2 = cos(T2);
		if(lt[1][2]>0) T1 = -asin(lt[1][2]/c2);
		else T1 = asin(lt[1][2]/c2);
		if(lt[2][3]>0) T3 = -asin(lt[2][3]/c2);
		else T3 = asin(lt[2][3]/c2);
		}
	if(c1c3<0) {
		lt[1][1]=-lt[1][1]; 
		lt[1][2]=-lt[1][2]; 
		lt[1][3]=-lt[1][3];
	}
}
///////////////////////////////////////////////////////////////////////////////
//
//  S y m m e t r i c  T e n s o r  O r d e r 2
//
const TensorValue SymmetricTensorOrder2::zeros[6] = {0};
const int SymmetricTensorOrder2::Tindex[4][4] = {{0,0,0,0},{0,1,4,6},{0,4,2,5},{0,6,5,3}};
CreateErrorHandler(SymmetricTensorOrder2);
SymmetricTensorOrder2::SymmetricTensorOrder2() 
{
	WhoClass("SymmetricTensorOrder2");
	T = new TensorValue [6] - 1;
	zero();
}
SymmetricTensorOrder2::~SymmetricTensorOrder2()
{
	delete [] (T+1);
}
SymmetricTensorOrder2 &SymmetricTensorOrder2::copy(TensorValue *t)
{
	WhoMethod("copy(TensorValue *t)");
	memcpy(T+1,t,6*sizeof(TensorValue)); 
	return (*this);
}
SymmetricTensorOrder2 &SymmetricTensorOrder2::operator=(SymmetricTensorOrder2 &
	tensor)
{
	WhoMethod("operator=(SymmetricTensorOrder2 &)");
	memcpy(T+1,tensor.T+1,6*sizeof(TensorValue));
	return (*this);
}
void SymmetricTensorOrder2::Rotate(const double &angle,const int axis)
{
	WhoMethod("Rotate(const double &,const int)");
	RotationTensor a(angle,axis);
	Rotate(a);
}
void SymmetricTensorOrder2::Rotate(RotationTensor &a)
{
	WhoMethod("Rotate(RotationTensor &)");
	VoigtRotationTensor V(a);
	V.Rotate((*this));
}
void SymmetricTensorOrder2::write(FILE *fp)
{
	WhoMethod("write(FILE *)");
	if(fwrite(T+1,sizeof(TensorValue),6,fp)!=6) {
		WhoMethod("write(FILE *)");
		FatalError("Cannot write tensor.");
	}
}
void SymmetricTensorOrder2::read(FILE *fp)
{
	WhoMethod("read(FILE *)");
	if(fread(T+1,sizeof(TensorValue),6,fp)!=6) {
		WhoMethod("read(FILE *)");
		FatalError("Cannot read tensor.");
	}
}
ostream& operator <<(ostream& o, SymmetricTensorOrder2& a)
{
	o << "|  " << a(1,1) << "  " << a(1,2) << "  " << a(1,3) << "|  " <<  endl;
	o << "|  " << a(2,1) << "  " << a(2,2) << "  " << a(2,3) << "|  " <<  endl;
	o << "|  " << a(3,1) << "  " << a(3,2) << "  " << a(3,3) << "|  " <<  endl;
	return o;
}
///////////////////////////////////////////////////////////////////////////////
//
//  V o i g t  R o t a t i o n  T e n s o r
//
VoigtRotationTensor::VoigtRotationTensor(RotationTensor &R)
{
	TensorValue **a = R.T;
	TensorValue a11 = a[1][1], a12 = a[1][2], a13 = a[1][3],
					a21 = a[2][1], a22 = a[2][2], a23 = a[2][3],
					a31 = a[3][1], a32 = a[3][2], a33 = a[3][3];
	V[0][0] = a11*a11; V[0][1] = a12*a12; V[0][2] = a13*a13; V[0][3] = 2*a11*a12; V[0][4] = 2*a12*a13; V[0][5] = 2*a11*a13;
	V[1][0] = a21*a21; V[1][1] = a22*a22; V[1][2] = a23*a23; V[1][3] = 2*a21*a22; V[1][4] = 2*a22*a23; V[1][5] = 2*a21*a23;
	V[2][0] = a31*a31; V[2][1] = a32*a32; V[2][2] = a33*a33; V[2][3] = 2*a31*a32; V[2][4] = 2*a32*a33; V[2][5] = 2*a31*a33;
	V[3][0] = a11*a21; V[3][1] = a12*a22; V[3][2] = a13*a23; V[3][3] = a11*a22+a12*a21; V[3][4] = a12*a23+a13*a22; V[3][5] = a11*a23+a13*a21;
	V[4][0] = a21*a31; V[4][1] = a22*a32; V[4][2] = a23*a33; V[4][3] = a21*a32+a22*a31; V[4][4] = a22*a33+a23*a32; V[4][5] = a21*a33+a23*a31;
	V[5][0] = a11*a31; V[5][1] = a12*a32; V[5][2] = a13*a33; V[5][3] = a11*a32+a12*a31; V[5][4] = a12*a33+a13*a32; V[5][5] = a11*a33+a13*a31;
	
}
void VoigtRotationTensor::Rotate(SymmetricTensorOrder2 &sym)
{
	int i,j;
	TensorValue r[6],sum,*lt=sym.T+1,*p;
	for(i=0;i<6;i++) {
		for(j=0,sum=0,p=V[i];j<6;j++) sum+=p[j]*lt[j];
		r[i]=sum;
	}
	memcpy(lt,r,6*sizeof(TensorValue));
}

//#define DEBUG_Tensor
#ifdef DEBUG_Tensor

main()
{
	double t[6]={1,2,3,4,5,6};
	int i;

	cout << "Testing Tensor...\n";
	cout << "Instantiating SymmetricTensorOrder2 a,b\n";
	SymmetricTensorOrder2 a,b;

	cout << "Setting up a\n";
	a(1,1)=10; a(1,2) = 20; a(1,3) = 30;
	a(2,2)=3; a(2,3) = 1;
	a(3,3)=5;
	cout << a;
	cout << "Setting up b\n";
	b(1,1)=10; b(1,2) = 20; b(1,3) = 30;
	b(2,2)=3; b(2,3) = 1;
	b(3,3)=5;
	cout << b;
	cout << "Testing vector assignment a(i)={1,2,3,4,5,6}\n";
	for(i=1;i<7;i++) cout << a(i) << endl;
	a.copy(t);
	cout << a;
	cout << "Testing rotation of a.Rotate(90,1) {90 degrees about axis 1}\n";
	a.Rotate(90,1);
	cout << a;
	cout << "Testing write/read to /tmp/junk.tensor\n";
	cout << "Write\n" << a;
	FILE *fp = fopen("/tmp/junk.tensor","wb");
	a.write(fp);
	fclose(fp);
	cout << "a = 0\n"; a = 0;
	if((fp = fopen("/tmp/junk.tensor","rb"))==NULL) cout << "Can't open file." << endl;
	a.read(fp);
	fclose(fp);
	cout << "Read\n" << a;
	

	cout << "Testing Tensor...\n";
	cout << "Instantiating Tensor a,b\n";
	Tensor A,B;

	cout << "Setting up A\n";
	A(1,1)=10; A(1,2) = 20; A(1,3) = 30;
	A(2,1)=5;  A(2,2)=3;    A(2,3) = 1;
	A(3,1)=3;  A(3,2)=17;   A(3,3)=5;
	cout << A;
	cout << "Setting up B\n";
	B(1,1)=10; B(1,2) = 20; B(1,3) = 30;
	B(2,1)=5;  B(2,2)=3;    B(2,3) = 1;
	B(3,1)=3;  B(3,2)=17;   B(3,3)=5;
	cout << B;
	cout << "Testing write/read to /tmp/junk.tensor\n";
	cout << "Write\n" << A;
	fp = fopen("/tmp/junk.tensor","wb");
	A.write(fp);
	fclose(fp);
	cout << "A = 0\n"; A = 0;
	if((fp = fopen("/tmp/junk.tensor","rb"))==NULL) cout << "Can't open file." << endl;
	A.read(fp);
	fclose(fp);
	cout << "Read\n" << A;
	cout << endl;
}
#endif

