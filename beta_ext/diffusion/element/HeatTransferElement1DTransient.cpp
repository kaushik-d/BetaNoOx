#include "stdafx.h"

#include <fstream>
#include <iomanip>
#include <string.h>

#include "TransientElementWorkspace.hpp"
#include "HeatTransferElement1DTransient.hpp"
#include "math/gaussQuadrature/hexQuadPoint.hpp"

//====================================================================
extern	int			verboseFlag;
//====================================================================
HeatTransferElement1Dtransient::
                       HeatTransferElement1Dtransient(void)
//====================================================================
{ 
	//if(traceFlag==1) {BETA_OUT<< "Creating HeatTransferElement1Dtransient element"<<endl;}
	numStrains       = 1;    // Eventually fix... these are the three temp. grad
	setIntegrationOrder(3);
	numberOfIntegrationDirections=1;
}
//====================================================================

HeatTransferElement1Dtransient::
                       ~HeatTransferElement1Dtransient(void)
{ 
	//BETA_OUT<< "Deleting HeatTransferElement1Dtransient element"<<endl; 
}

//====================================================================
void HeatTransferElement1Dtransient::
							 addToFe(int numDof, int numStrains,
							 InterpAndSolDerivatives *id_ptr,
							 double *flux,
							 double integrationFactor,
							 double *force)
{
	InterpAndSolDerivatives &id=*id_ptr;
	double  *p_S_x1=id.p_S_x1;		
	  for(int i = 0; i<numDof; i++)
		 {
        force[i] += ( flux[0] * p_S_x1[i] ) * integrationFactor;
		 }
}
//==================================================================
void HeatTransferElement1Dtransient::
									gradientsOfFieldVariables(InterpAndSolDerivatives *id_ptr)
{
	InterpAndSolDerivatives &id=*id_ptr;
	int i;
	int numberOfInterp = id.numberOfInterp;
	double *S=id.S,
		 *p_S_Lx1=id.p_S_Lx1,   *p_S_x1=id.p_S_x1;
   
	double *xDisp=id.nodalValues_u1;

	double p_u1_x1;

	p_u1_x1 = 0.0;

	
	for(i=0; i< numberOfInterp; i++)
	{
				p_u1_x1 += p_S_x1[i] * xDisp[i];
	}
	id.p_u1_x1 = p_u1_x1;
}  // end 
//==================================================================
void HeatTransferElement1Dtransient::addToKe(Matrix* cMat,InterpAndSolDerivatives *id_ptr,
									  double integrationFactor)       
{	
	InterpAndSolDerivatives &id=*id_ptr;
double *S=id.S,
       *p_S_x1=id.p_S_x1;
HeatTransferMaterial* material=(HeatTransferMaterial*)BasicElement::material;
	//saturation value already factored into the cMat

double k11;
k11 = (*cMat)(0,0);

//BETA_OUT<<" k11,k22,k33 "<<k11<<"  "<<k22<<"  "<<k33<<endl;

int i,j;
double p_S_x1_i = 0.0;
	for(i=0; i<numDof; i++) {
		p_S_x1_i = p_S_x1[i];
		for(j=i; j<numDof; j++)	{
			(*Ke)[i][j] += (p_S_x1[j]*k11*p_S_x1_i) * integrationFactor;
		}
	 }
}  //end of addToKe
//==========================================================
