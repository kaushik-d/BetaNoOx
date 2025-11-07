#include "stdafx.h"

#include <fstream>
#include <iomanip>
#include <string.h>

#include "TransientElementWorkspace.hpp"
#include "HeatTransferElement2DTransient.hpp"
#include "math/gaussQuadrature/hexQuadPoint.hpp"

//====================================================================
extern	int			verboseFlag;
//====================================================================
HeatTransferElement2Dtransient::HeatTransferElement2Dtransient(void)
{ 
	//if(traceFlag==1) {BETA_OUT<< "Creating HeatTransferElement2Dtransient element"<<endl;}
	numStrains       = 2;    // Eventually fix... these are the three temp. grad
	setIntegrationOrder(3);
	numberOfIntegrationDirections=2;
}
//====================================================================
HeatTransferElement2Dtransient::~HeatTransferElement2Dtransient(void)
{ 
	//BETA_OUT<< "Deleting HeatTransferElement2Dtransient element"<<endl; 
}
//====================================================================
void HeatTransferElement2Dtransient::addToFe(int numDof, int numStrains,
							 InterpAndSolDerivatives *id_ptr,
							 double *flux,
							 double integrationFactor,
							 double *force)
{
	InterpAndSolDerivatives &id=*id_ptr;
	double  *p_S_x1=id.p_S_x1,
		  *p_S_x2=id.p_S_x2;		
	  for(int i = 0; i<numDof; i++)
		 {
        force[i] += ( flux[0] * p_S_x1[i] + 
			          flux[1] * p_S_x2[i] ) * integrationFactor;
		 }
  }
//==================================================================
void HeatTransferElement2Dtransient::gradientsOfFieldVariables(InterpAndSolDerivatives *id_ptr)
{
	InterpAndSolDerivatives &id=*id_ptr;
	int i;
	int numberOfInterp = id.numberOfInterp;
	double *S=id.S,
		 *p_S_Lx1=id.p_S_Lx1,   *p_S_x1=id.p_S_x1, 
		 *p_S_Lx2=id.p_S_Lx2,   *p_S_x2=id.p_S_x2;

	double *xDisp=id.nodalValues_u1;

	double p_u1_x1, p_u1_x2;

	p_u1_x1 = p_u1_x2 = 0.;

	
	for(i=0; i< numberOfInterp; i++)
	{
				p_u1_x1 += p_S_x1[i] * xDisp[i];
				p_u1_x2 += p_S_x2[i] * xDisp[i];

	}
	id.p_u1_x1 = p_u1_x1;
	id.p_u1_x2 = p_u1_x2;
}
//==================================================================
void HeatTransferElement2Dtransient::addToKe(Matrix* cMat,InterpAndSolDerivatives *id_ptr,
									  double integrationFactor)       
{	
	InterpAndSolDerivatives &id=*id_ptr;
double *S=id.S,
       *p_S_x1=id.p_S_x1,    
	   *p_S_x2=id.p_S_x2;	
HeatTransferMaterial* material=(HeatTransferMaterial*)BasicElement::material;
	//saturation value already factored into the cMat

double k11,k22;
k11 = (*cMat)(0,0);
k22 = (*cMat)(1,1);

double k12;
k12 = (*cMat)(0,1);

//BETA_OUT<<" k11,k22,k33 "<<k11<<"  "<<k22<<"  "<<k33<<endl;

// Limit to orthotropic for now
// double k12=0,k13=0,k23=0;

int i,j;

double p_S_x1_i=0.0;
double p_S_x2_i=0.0;

	for(i=0; i<numDof; i++) {
		p_S_x1_i = p_S_x1[i];
		p_S_x2_i = p_S_x2[i];
		for(j=i; j<numDof; j++)	{
			(*Ke)[i][j] += (
				   (p_S_x1[j]*k11+p_S_x2[j]*k12)*p_S_x1_i+
			       (p_S_x1[j]*k12+p_S_x2[j]*k22)*p_S_x2_i )
                   * integrationFactor;
		}
	 }

}  //end of addToKe
//==========================================================
