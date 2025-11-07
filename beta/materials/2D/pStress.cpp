#include "stdafx.h"

// Merge with PlaneElastic ???

#include <iostream>
#include <iomanip>
#include <cmath>    /* find out what this does !! */

#include "pStress.hpp"
#include "math/matrix.hpp"

void getRotatedProperties(Matrix &newC, Matrix &oldC,
								  const double theta );

//====================================================
PlaneStress::PlaneStress()
{
}
//====================================================
PlaneStress::~PlaneStress(void)          
{
}
//--------------------------------------------------------------------
void PlaneStress::initialize()
{
	numStrains = 3;
	setType("planeStress");
	CommonElasticMaterial();
	// BETA_OUT<<"Created plane stress material object"<<endl;
}
//====================================================
void PlaneStress::print(ostream & outStream)
{
outStream<<"***************"<<endl;
outStream<<      name      <<endl;
outStream<<"***************"<<endl;

outStream << "E11 : "<< e11 <<endl;
outStream << "E22 : "<< e22 <<endl;
outStream << "nu12: "<< nu12<<endl;
outStream << "G12 : "<< G12 <<endl;
outStream << "a11 : "<< a11 <<endl;
outStream << "a22 : "<< a22 <<endl;
(*cMat).PrintMatrix(" Local stiffness matrix cMat  ");
}

//====================================================

bool PlaneStress::read(istream * inStream)
{
(* inStream) >>e11>>e22>>nu12>>G12>>a11>>a22;
//static int echoStatus;
//static ofstream echo;
//if(echoStatus !=999) {echo.open("echo.mat", ios::out );	}
//echoStatus =999;

BETA_OUT<<"***************"<<endl;
BETA_OUT<<      name       <<endl;
BETA_OUT<<"***************"<<endl;

BETA_OUT << "E11 : "<< e11 <<endl;
BETA_OUT << "E22 : "<< e22 <<endl;
BETA_OUT << "nu12: "<< nu12<<endl;
BETA_OUT << "G12 : "<< G12 <<endl;
BETA_OUT << "a11 : "<< a11 <<endl;
BETA_OUT << "a22 : "<< a22 <<endl;

	  REAL s11,s12,s22,s66,den;

	  s11 = 1/e11;
	  s12 =  -nu12/e11;
	  s22 = 1./e22;
	  s66 = 1./G12;
	  den =  s11 * s22 - s12 * s12;


	  (*cMat)(0,0) =  s22/den;
	  (*cMat)(1,0) = -s12/den;
	  (*cMat)(2,0) =  0.;
	  (*cMat)(0,1) = -s12/den;
	  (*cMat)(1,1) =  s11/den;
	  (*cMat)(2,1) =  0.;
	  (*cMat)(0,2) =  0.;
	  (*cMat)(1,2) =  0;
	  (*cMat)(2,2) =  1./s66;

	 // cInv(0,0) =  s11;
	 // cInv(1,0) =  s12;
	 // cInv(2,0) =  0.;
	 // cInv(0,1) =  s12;
	 // cInv(1,1) =  s22;
	 // cInv(2,1) =  0.;
	 // cInv(0,2) =  0.;
	 // cInv(1,2) =  0;
	 // cInv(2,2) =  s66;

/*
// Temporary hardwire for isotropic plane strain

      (*cMat)(0,0) =  e11*(1-nu12)/(1- nu12 - 2.*nu12*nu12);
	  (*cMat)(1,0) = nu12/(1-nu12) * (*cMat)(0,0);
	  (*cMat)(2,0) =  0.;
	  (*cMat)(0,1) = nu12/(1-nu12) * (*cMat)(0,0);
	  (*cMat)(1,1) = e11*(1-nu12)/(1- nu12 - 2.*nu12*nu12) ;
	  (*cMat)(2,1) =  0.;
	  (*cMat)(0,2) =  0.;
	  (*cMat)(1,2) =  0;
	  (*cMat)(2,2) =  e11/(2.*(1+nu12));

*/

	 // cInv.PrintMatrix(" Local compliance matrix cInv ");
	  (*cMat).PrintMatrix(" Local stiffness matrix cMat  ");

(*globalCmat) = (*cMat);  //check this !!
return(true);
}




 /* ****************************************************************** */

void getRotatedProperties(Matrix &newC, Matrix &oldC,
								  const double theta)
{
	  REAL ct, st, ct2,st2,st3,ct3,st4,ct4,thetaRad;

	  thetaRad = theta * 3.1415927 /180.;

	  /* Warning: this routine assumes that "Old" is ortho.   */

	  st  = sin(thetaRad);
	  ct  = cos(thetaRad);
	  st2 = st * st;
	  ct2 = ct * ct;
	  st3 = st * st2;
	  ct3 = ct * ct2;
	  st4 = st2 * st2;
	  ct4 = ct2 * ct2;

#define oldC(i,j) oldC[i][j]

	  newC(2,2) = ct2*(oldC(0,0) - 2*oldC(0,1) + oldC(1,1)
	  - 2*oldC(2,2))* st2 + oldC(2,2)*(ct4 + st4);
	  newC(1,2) = ct3*(oldC(0,1) - oldC(1,1) + 2*oldC(2,2))*st +
	  ct*(oldC(0,0) - oldC(0,1) - 2*oldC(2,2))*st3;
	  newC(1,1) = ct4*oldC(1,1) +
	  2*ct2*(oldC(0,1) + 2*oldC(2,2))*st2 + oldC(0,0)*st4;
	  newC(0,2) = ct3*(oldC(0,0) - oldC(0,1) - 2*oldC(2,2))*st +
	  ct*(oldC(0,1) - oldC(1,1) + 2*oldC(2,2))*st3;
	  newC(0,1) = ct2*(oldC(0,0) + oldC(1,1) - 4*oldC(2,2))*st2 +
	  oldC(0,1)*(ct4 + st4);
	  newC(0,0) = ct4*oldC(0,0) + 2*ct2*(oldC(0,1) +
	  2*oldC(2,2))*st2 + oldC(1,1)*st4;
	  newC(1,0) = newC(0,1);
	  newC(2,0) = newC(0,2);
	  newC(2,1) = newC(1,2);

#undef oldC
}

 /* ****************************************************************** */
