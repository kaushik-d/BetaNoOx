#include "stdafx.h"

#include "ElasticityElement2D.hpp"
#include "elements/InterpAndSolDerivatives.hpp"
#include "math/gaussQuadrature/quadQuadPoint.hpp"

  
extern int verboseFlag;
//================================================================
ElasticityElement2D::ElasticityElement2D(void)
{ //if(traceFlag==1)
    //  {BETA_OUT<< "Creating ElasticityElement2D element"<<endl;}
	  numStrains       = 3;
	  setIntegrationOrder(2);
     numberOfIntegrationDirections=2;
}
//=====================================================================
ElasticityElement2D::~ElasticityElement2D(void)
{  // Finish !!
}
//====================================================================
 void ElasticityElement2D::addToKe(Matrix *cMat,InterpAndSolDerivatives *id_ptr,
				        double integrationFactor,double *stress)
{
	InterpAndSolDerivatives &id=*id_ptr;
	double *p_S_x1=id.p_S_x1,    *p_S_x2=id.p_S_x2;
	int numberOfInterp = id.numberOfInterp;
	double C_x_B[3][36];   // Dimensioned for up to 3 strains
	                       // and 36 dof (this code is only good
	                       // for 3 strains

	int    ii,jj;
	double * b0 = &b[0][0];
	double * b1 = &b[1][0];
	double * b2 = &b[2][0];

	for(ii=0; ii<numDof; ii++)
		{double b0ii = b[0][ii];
		 double b1ii = b[1][ii];
		 double b2ii = b[2][ii];
		 for(jj=0; jj<numStrains;jj++)
			  {C_x_B[jj][ii] = b0ii * (*cMat)[0][jj] +
			   b1ii * (*cMat)[1][jj] +
			   b2ii * (*cMat)[2][jj] ;
			  }
		}

		for(ii=0; ii<numDof; ii++)
			{double p0ii = C_x_B[0][ii];
			 double p1ii = C_x_B[1][ii];
			 double p2ii = C_x_B[2][ii];
			 for(jj=ii; jj<numDof; jj++)
				 {
					(*Ke)[ii][jj] += ( p0ii * b0[jj] +
									       p1ii * b1[jj] +
									       p2ii * b2[jj] ) * integrationFactor;
				 }
			 }

if(analysisType==GEOMETRIC_NONLINEAR)    //add Ksigma
 {
 double b11, b11dj,q; 
 int mb1,nb1, nodei,nodej;
 
 #define ndofpn  2
 
 for(nodei= 0;  nodei<numberOfInterp; nodei++)
    {
	  mb1= nodei * ndofpn  ; 
     for(nodej= nodei;  nodej<numberOfInterp; nodej++)
        {
         b11 =stress[0] * p_S_x1[nodei]* p_S_x1[nodej] +
			  stress[1] * p_S_x2[nodei]* p_S_x2[nodej]    ;
          q = stress[2] *(p_S_x1[nodei]*p_S_x2[nodej]+
                          p_S_x2[nodei]*p_S_x1[nodej])  ;
	       b11 = b11 + q ;
                nb1= nodej * ndofpn  ;
                b11dj = b11 * integrationFactor       ;
                (*Ke)[mb1][nb1]     +=  b11dj ;
                (*Ke)[mb1+1][nb1+1] +=  b11dj ;
			}		 
    }	 	 
	}  // end of Ksigma
#undef ndofpn
	
}  // end of addToKe
//====================================================================
 void ElasticityElement2D::calculateBmatrix(InterpAndSolDerivatives *id_ptr, double **b)
{
	InterpAndSolDerivatives &id=*id_ptr;
double *p_S_x1=id.p_S_x1, *p_S_x2= id.p_S_x2;
int i;
int c1 = 0;
int c2 = 1;
int numberOfInterp = id.numberOfInterp;
for( i=0; i< numberOfInterp; i++, c1+=2, c2+=2 )
  {
	b[0][c1] = p_S_x1[i];       b[0][c2] = 0.    ;
	b[1][c1] = 0.   ;           b[1][c2] = p_S_x2[i] ;
	b[2][c1] = p_S_x2[i];       b[2][c2] = p_S_x1[i] ;
  }


 if(analysisType == GEOMETRIC_NONLINEAR) //Add nonlinear terms
 {
  int i;
  int c1 = 0;
  int c2 = 1;

  for( i=0; i< numberOfInterp; i++, c1+=2, c2+=2 )
  {
   double dg1=p_S_x1[i];
   double dg2=p_S_x2[i];
	b[0][c1] += id.p_u1_x1*dg1;                 b[0][c2] += id.p_u2_x1*dg1 ;
	b[1][c1] += id.p_u1_x2*dg2;                 b[1][c2] += id.p_u2_x2*dg2 ;
	b[2][c1] += id.p_u1_x1*dg2 + id.p_u1_x2*dg1 ; b[2][c2] += id.p_u2_x1*dg2 + id.p_u2_x2*dg1;
  }
 }
}//end of linearBmat
//======================================================================
void ElasticityElement2D::calculateStrain(InterpAndSolDerivatives *id_ptr,double *strain)
{//engineering shear
	InterpAndSolDerivatives &id=*id_ptr;
	strain[0] = id.p_u1_x1 ;
	strain[1] = id.p_u2_x2 ;
	strain[2] = id.p_u1_x2 + id.p_u2_x1;
 if(analysisType == GEOMETRIC_NONLINEAR) //Add nonlinear terms
  {
   strain[0] = id.p_u1_x1 + .5 * ( id.p_u1_x1 * id.p_u1_x1 
		                         + id.p_u2_x1 * id.p_u2_x1 );
   strain[1] = id.p_u2_x2 + .5 * ( id.p_u1_x2 * id.p_u1_x2 
		                         + id.p_u2_x2 * id.p_u2_x2 );
   strain[2] = id.p_u1_x2 + id.p_u2_x1 + id.p_u1_x1 * id.p_u1_x2 
		                  + id.p_u2_x1 * id.p_u2_x2  ;
  }
}
//======================================================================
void ElasticityElement2D::derivativesOfFields(InterpAndSolDerivatives *id_ptr)
{
	InterpAndSolDerivatives &id=*id_ptr;
int i;
int numberOfInterp = id.numberOfInterp;
double *S=id.S, *p_S_Lx1=id.p_S_Lx1, *p_S_Lx2=id.p_S_Lx2,
		          *p_S_x1=id.p_S_x1,    *p_S_x2=id.p_S_x2;
double *xDisp = id.nodalValues_u1;
double *yDisp = id.nodalValues_u2;

double p_u1_x1, p_u1_x2, p_u2_x1, p_u2_x2;
p_u1_x1 = p_u1_x2 = p_u2_x1 = p_u2_x2 = 0.;

for(i=0; i< numberOfInterp; i++)
{
				p_u1_x1 += p_S_x1[i] * xDisp[i];
				p_u1_x2 += p_S_x2[i] * xDisp[i];
				p_u2_x1 += p_S_x1[i] * yDisp[i];
				p_u2_x2 += p_S_x2[i] * yDisp[i];
}
id.p_u1_x1 = p_u1_x1;
id.p_u1_x2 = p_u1_x2;
id.p_u2_x1 = p_u2_x1;
id.p_u2_x2 = p_u2_x2;
} // end of derivativesOfFields for 2D
//=====================================================================
void ElasticityElement2D::getVolAvgValues(double *uGlobal, double *volAvgStress,
	         double *volAvgStrain, double *SED, double *volume)
{
// Calculate volume averaged values (strain and stress) and volume
// (Jae Noh, 6/23/98)
// copied into ElasticityElement2D by JV - 060204
#define maxStresses 6
//InterpAndSolDerivatives dispDeriv; 
GaussPointList				gaussPointList;
InterpAndSolDerivatives       interpData(numNodesPerElement,numNodesPerElement); interpData.allocate();
double *elemDisp[3];
//Matrix *cMat, *sMat, *sMatOrtho, *sMatNonOrtho;
Matrix *cMat, *sMatOrtho, *sMatNonOrtho;
double strain[6];//, SEDtmp0[6], SEDtmp1[6], SEDtmp2[6], SED2[3];
double *stress;
double *rotStress;
double **quadPointStresses   = bag->quadPointStresses;
double **rotatedStresses     = bag->rotatedStresses;

 //------------------------------------------------
 sMatOrtho = new Matrix(6,6);
 sMatNonOrtho = new Matrix(6,6);
 elemDisp[0] = interpData.nodalValues_u1; 
 elemDisp[1] = interpData.nodalValues_u2; 
 elemDisp[2] = interpData.nodalValues_u3; 

 interpData.numberOfInterp = numNodesPerElement;
 extractNodalCoordinates(interpData.xCoor,interpData.yCoor,
	                      interpData.zCoor);
 if( uGlobal != NULL) {extractSolution( uGlobal, elemDisp);} 
 
 initializeArrays();
 getQuadraturePoints(gaussPointList); 
 int totalNumIPs = gaussPointList.totalNumIPs;
//...............................................................
//                      Integration loop
//...............................................................
double integrationFactor, weight;

    *volume=0;  //xtang 051302

	int is;
for(is=0; is<maxStresses; is++)
{   volAvgStress[is] = volAvgStrain[is]=0.; }
	SED[0]=SED[1]=SED[2]=0.0;

for(int ip=0;ip<totalNumIPs;ip++)
{
	interpolation(ip, gaussPointList, &interpData);
    jacobian(&interpData);
	weight = gaussPointList.weightList[ip];
	integrationFactor = interpData.detJ * weight * getIntegrationFactor();
    *volume += integrationFactor;

	ElasticMaterial* material = (ElasticMaterial*)BasicElement::material;

    //OrientationAtIP->interpolate_rotation(&interpData,ElementNodeRotations);
    OrientationAtIP->interpolate_rotation_by_angle(&interpData,ElementNodeRotations);
    //OrientationAtIP->interpolate_old_beta(&interpData,ElementNodeRotations,material->getMaterialRotation());

    material->rotateMaterialToGlobal(OrientationAtIP);

	cMat = material->getPointerToGlobalCmat();  // cMat->print("cmat");

	 //if(uGlobal!=NULL){dispGradients(interpData,dispDeriv);}   
	//calculateBmatrix(interpData,dispDeriv,b);
	if(uGlobal!=NULL){derivativesOfFields(&interpData);}               
	calculateBmatrix(&interpData,b);

	stress      = quadPointStresses[ip];
	rotStress   = rotatedStresses[ip];

	//calculateStrain(interpData,dispDeriv,strain); 
	calculateStrain(&interpData, strain); 

	// xtang 122000
//    material->calculateStress(stress,strain,getMaterialAngle(interpData),
//				materialRotationAxis); 
//    material->rotateStress(rotStress,stress,-getMaterialAngle(interpData),
//									 materialRotationAxis, true);
    
    material->calculateStress(stress,
                              strain,
                              isv[ip],
                              OrientationAtIP);

    OrientationAtIP->rotate_Voigt_stress_from_global_to_local(stress,rotStress);

   // calculate strain energy density
		//FIX LATER !!!
/*
   //cMat->print("cmat");
	sMat = cMat;     // xtang:  local coords ?? see 
	sMat->invert();
	for (i=0; i<3; i++){
      (*sMatOrtho)[i+3][i+3]=(*sMat)[i+3][i+3];
		for(j=0; j<3; j++){
	       (*sMatOrtho)[i][j]=(*sMat)[i][j];
          (*sMatNonOrtho)[i][j+3]=(*sMat)[i][j+3];
          (*sMatNonOrtho)[i+3][j]=(*sMat)[i+3][j];
		}
	}
   (*sMatNonOrtho)[3][4]=(*sMat)[3][4];   
   (*sMatNonOrtho)[3][5]=(*sMat)[3][5];   
   (*sMatNonOrtho)[4][3]=(*sMat)[4][3];   
   (*sMatNonOrtho)[4][5]=(*sMat)[4][5];   
   (*sMatNonOrtho)[5][3]=(*sMat)[5][3];   
   (*sMatNonOrtho)[5][4]=(*sMat)[5][4];   
   //sMat->print("sMat");
   //sMatOrtho->print("sMatOrtho");
   //sMatNonOrtho->print("sMatNonOrtho");

	for(i=0; i<6; i++){
		SEDtmp0[i]=SEDtmp1[i]=SEDtmp2[i]=0.0;
		for(j=0; j<6; j++){
		    SEDtmp0[i] += rotStress[j]* (*sMat)[j][i];
		    SEDtmp1[i] += rotStress[j]* (*sMatOrtho)[j][i];
		    SEDtmp2[i] += rotStress[j]* (*sMatNonOrtho)[j][i];
		}
	}
	SED2[0]=SED2[1]=SED2[2]=0.0;
	for(i=0; i<6; i++){
		SED2[0] += SEDtmp0[i] * rotStress[i];
		SED2[1] += SEDtmp1[i] * rotStress[i];
		SED2[2] += SEDtmp2[i] * rotStress[i];
	}
	for(i=0; i<3; i++) SED[i] += 0.5 * SED2[i] * integrationFactor;
*/
   for(is=0; is<maxStresses; is++){
	   volAvgStress[is] += integrationFactor * stress[is];
	   volAvgStrain[is] += integrationFactor * strain[is];
	}
}  // end of integration loop

 delete sMatOrtho;    // xt062901;
 delete sMatNonOrtho; // xt062901;

return;
}
//==================================================================
