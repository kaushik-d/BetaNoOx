#include "stdafx.h"

#include <iostream>
#include <iomanip>
#include <cmath>    /* find out what this does !! */

#include "pStrain.hpp"
#include "math/matrix.hpp"

//====================================================
PlaneStrain::PlaneStrain()
{
}
//====================================================
PlaneStrain::~PlaneStrain(void)          
{
}
//--------------------------------------------------------------------
void PlaneStrain::initialize()
{
	numStrains = 3;
	setType("planeStrain");
	CommonElasticMaterial();
	// BETA_OUT<<"Created plane stress material object"<<endl;
}
//====================================================
void PlaneStrain::print(ostream & outStream)
{
outStream<<"***************"<<endl;
outStream<<      name      <<endl;
outStream<<"***************"<<endl;

outStream << "E11 : "<< e11 <<endl;
outStream << "E22 : "<< e22 <<endl;
outStream << "E33 : "<< e33 <<endl;
outStream << "nu12: "<< nu12<<endl;
outStream << "nu23: "<< nu23<<endl;
outStream << "nu13: "<< nu13<<endl;
outStream << "G12 : "<< G12 <<endl;
outStream << "G23 : "<< G23 <<endl;
outStream << "G13 : "<< G13 <<endl;
(*cMat).PrintMatrix(" Local stiffness matrix cMat  ");
}

//====================================================

bool PlaneStrain::read(istream * inStream)
{
(* inStream) >>e11>>e22>>e33>>nu12>>nu23>>nu13>>G12>>G23>>G13;
//static int echoStatus;
//static ofstream echo;
//if(echoStatus !=999) {echo.open("echo.mat", ios::out );	}
//echoStatus =999;

BETA_OUT << "E11 : "<< e11 <<endl;
BETA_OUT << "E22 : "<< e22 <<endl;
BETA_OUT << "E33 : "<< e33 <<endl;
BETA_OUT << "nu12: "<< nu12<<endl;
BETA_OUT << "nu23: "<< nu23<<endl;
BETA_OUT << "nu13: "<< nu13<<endl;
BETA_OUT << "G12 : "<< G12 <<endl;
BETA_OUT << "G23 : "<< G23 <<endl;
BETA_OUT << "G13 : "<< G13 <<endl;

	  REAL delta,nu21,nu32,nu31;

	  nu21 = nu12/e11*e22;
	  nu32 = nu23/e22*e33;
	  nu31 = nu13/e11*e33;

	  delta = (1.0 - nu12*nu21 - nu23*nu32 - nu13*nu31 - 2.0*nu12*nu23*nu31)/(e11*e22*e33);

	  (*cMat)(0,0) = (1.0 - nu23*nu32)/(e22*e33*delta);
	  (*cMat)(0,1) = (nu21 + nu31*nu23)/(e22*e33*delta);
	  (*cMat)(0,2) = 0.;
	  (*cMat)(1,0) = (nu12+nu13*nu32)/(e33*e11*delta);
	  (*cMat)(1,1) = (1.0 - nu31*nu13)/(e33*e11*delta);
	  (*cMat)(1,2) = 0.;
	  (*cMat)(2,0) = 0.;
	  (*cMat)(2,1) = 0.;
	  (*cMat)(2,2) = G12;

	  (*cMat).PrintMatrix(" Local stiffness matrix cMat  ");

(*globalCmat) = (*cMat);  //check this !!
return(true);
}

/* ****************************************************************** */
