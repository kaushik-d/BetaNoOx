#include "stdafx.h"
#include "utility/utility.h"

#include "ElasticityElement3D.hpp"

#include "elements/InterpAndSolDerivatives.hpp"
#include "math/gaussQuadrature/hexQuadPoint.hpp"
//#include "FEMLibrary/models/OutputOptions.hpp"
#include "math/tensor.hpp"

#define FF1 setw(13)<<setprecision(4)<<setiosflags(ios::scientific)

void rotateEngineeringVoigtStrain(double *rotatedStress, double *stress,
                       const double &theta,const int idir);

//====================================================================
extern	int			verboseFlag;
//====================================================================
ElasticityElement3D::ElasticityElement3D(void):
IsoElement()
{ 
		//if(traceFlag==1)
      //{BETA_OUT<< "Creating ElasticityElement3D element"<<endl;}
	numStrains       = 6;
	setIntegrationOrder(3);
	numberOfIntegrationDirections=3;

	reduceCMatFlag=0;
	b=0;

	ForceSingleNonZeroStressComponent=-1;
	bag=0;

}
//====================================================================
ElasticityElement3D::~ElasticityElement3D(void)
{ 
	//if(traceFlag==1){cout<< "Deleting ElasticityElement3D element"<<endl;}
}
//====================================================================
 void ElasticityElement3D::
              addToFe(int numDof, int numStrains, double **b,
							 double *stress,
							 double integrationFactor,
							 double *force)
//====================================================================
{
	 //Not optimized
	  for(int i = 0; i<numDof; i++)
		 {
		 double sum = 0.;
		 for(int j=0; j<numStrains; j++){sum +=  b[j][i] * stress[j];}
			 force[i] += sum * integrationFactor;
		 }
}
//=====================================================================  
void ElasticityElement3D::
               calculateBmatrix(InterpAndSolDerivatives *id_ptr,
										  double **b)
//====================================================================
{    
	InterpAndSolDerivatives &id=*id_ptr;
//Order: epsx epsy epsz epsXY epsYZ epsXZ  (Engineering shear)
double  *p_S_x1=id.p_S_x1,
		  *p_S_x2=id.p_S_x2,
		  *p_S_x3=id.p_S_x3;		      		  
int i;
int c1 = 0;
int c2 = 1;
int c3 = 2;
int numberOfInterp = id.numberOfInterp;
for( i=0; i< numberOfInterp; i++, c1+=3, c2+=3, c3+=3)
  {
	b[0][c1] = p_S_x1[i]; b[0][c2] = 0.       ; b[0][c3] = 0.        ;
	b[1][c1] = 0.		  ; b[1][c2] = p_S_x2[i]; b[1][c3] = 0.        ;
	b[2][c1] = 0.		  ; b[2][c2] = 0.       ; b[2][c3] = p_S_x3[i] ;

	b[3][c1] = p_S_x2[i]; b[3][c2] = p_S_x1[i]; b[3][c3] = 0.        ;
	b[4][c1] = 0.       ; b[4][c2] = p_S_x3[i]; b[4][c3] = p_S_x2[i] ;
	b[5][c1] = p_S_x3[i]; b[5][c2] = 0.       ; b[5][c3] = p_S_x1[i] ;
  }   


 if(analysisType == GEOMETRIC_NONLINEAR) //Add nonlinear terms
 {
  //This has not been checked !!!!!
  int i;
  int c1 = 0;
  int c2 = 1;
  int c3 = 2;

  for( i=0; i< numberOfInterp; i++, c1+=3, c2+=3, c3+=3 )
  {
   double dg1=p_S_x1[i];
   double dg2=p_S_x2[i];
   double dg3=p_S_x3[i];
	b[0][c1] += id.p_u1_x1*dg1;                
	b[1][c1] += id.p_u1_x2*dg2;                
	b[2][c1] += id.p_u1_x3*dg3 ;
	b[3][c1] += id.p_u1_x1*dg2 + id.p_u1_x2*dg1;                
	b[4][c1] += id.p_u1_x2*dg3 + id.p_u1_x3*dg2;                 
	b[5][c1] += id.p_u1_x1*dg3 + id.p_u1_x3*dg1; 

   b[0][c2] += id.p_u2_x1*dg1 ;
   b[1][c2] += id.p_u2_x2*dg2 ;
   b[2][c2] += id.p_u2_x3*dg3 ;
	b[3][c2] += id.p_u2_x1*dg2 + id.p_u2_x2*dg1;                
	b[4][c2] += id.p_u2_x2*dg3 + id.p_u2_x3*dg2;                 
	b[5][c2] += id.p_u2_x1*dg3 + id.p_u2_x3*dg1;

   b[0][c3] += id.p_u3_x1*dg1 ;
   b[1][c3] += id.p_u3_x2*dg2 ;
   b[2][c3] += id.p_u3_x3*dg3 ;
	b[3][c3] += id.p_u3_x1*dg2 + id.p_u3_x2*dg1;                
	b[4][c3] += id.p_u3_x2*dg3 + id.p_u3_x3*dg2;                 
	b[5][c3] += id.p_u3_x1*dg3 + id.p_u3_x3*dg1; 

  }//end of for i
 }//end of if

} //end of calculateBmatrix3D  

 //======================================================================
void ElasticityElement3D::calculateStrain(InterpAndSolDerivatives *id_ptr, double *strain)
{//engineering shear
	InterpAndSolDerivatives &id=*id_ptr;

	strain[0] = id.p_u1_x1 ;
	strain[1] = id.p_u2_x2 ;
	strain[2] = id.p_u3_x3 ;
	strain[3] = id.p_u1_x2 + id.p_u2_x1;
	strain[4] = id.p_u2_x3 + id.p_u3_x2;
	strain[5] = id.p_u1_x3 + id.p_u3_x1;


 if(analysisType == GEOMETRIC_NONLINEAR) //Add nonlinear terms
  {
   strain[0] += .5 * ( id.p_u1_x1 * id.p_u1_x1 
		               + id.p_u2_x1 * id.p_u2_x1
					      + id.p_u3_x1 * id.p_u3_x1);
   strain[1] += .5 * ( id.p_u1_x2 * id.p_u1_x2 
		               + id.p_u2_x2 * id.p_u2_x2 
					      + id.p_u3_x2 * id.p_u3_x2);
   strain[2] += .5 * ( id.p_u1_x3 * id.p_u1_x3 
		               + id.p_u2_x3 * id.p_u2_x3
					      + id.p_u3_x3 * id.p_u3_x3);

   strain[3] += id.p_u1_x1 * id.p_u1_x2 +
	             id.p_u2_x1 * id.p_u2_x2 + 
				    id.p_u3_x1 * id.p_u3_x2 ;  
   strain[4] += id.p_u1_x2 * id.p_u1_x3 +
	             id.p_u2_x2 * id.p_u2_x3 + 
				    id.p_u3_x2 * id.p_u3_x3 ;  
   strain[5] += id.p_u1_x1 * id.p_u1_x3 +
	             id.p_u2_x1 * id.p_u2_x3 + 
				    id.p_u3_x1 * id.p_u3_x3 ;  

  }//end of if 

}// end of method
//======================================================================
void ElasticityElement3D::
        derivativesOfFields(InterpAndSolDerivatives *id_ptr)
{
	InterpAndSolDerivatives &id=*id_ptr;
int i;
int numberOfInterp = id.numberOfInterp;
double *S=id.S,
		 *p_S_Lx1=id.p_S_Lx1,   *p_S_x1=id.p_S_x1, 
		 *p_S_Lx2=id.p_S_Lx2,   *p_S_x2=id.p_S_x2, 
		 *p_S_Lx3=id.p_S_Lx3,   *p_S_x3=id.p_S_x3;
   
double *xDisp=id.nodalValues_u1,
		 *yDisp=id.nodalValues_u2,
       *zDisp=id.nodalValues_u3;

double p_u1_x1, p_u1_x2, p_u1_x3;
double p_u2_x1, p_u2_x2, p_u2_x3;
double p_u3_x1, p_u3_x2, p_u3_x3;
p_u1_x1 = p_u1_x2 = p_u1_x3 = 0.;
p_u2_x1 = p_u2_x2 = p_u2_x3 = 0.;
p_u3_x1 = p_u3_x2 = p_u3_x3 = 0.;

for(i=0; i< numberOfInterp; i++)
{
				p_u1_x1 += p_S_x1[i] * xDisp[i];
				p_u1_x2 += p_S_x2[i] * xDisp[i];
				p_u1_x3 += p_S_x3[i] * xDisp[i];
				p_u2_x1 += p_S_x1[i] * yDisp[i];
				p_u2_x2 += p_S_x2[i] * yDisp[i];
				p_u2_x3 += p_S_x3[i] * yDisp[i];
				p_u3_x1 += p_S_x1[i] * zDisp[i];
				p_u3_x2 += p_S_x2[i] * zDisp[i];
				p_u3_x3 += p_S_x3[i] * zDisp[i];
}
id.p_u1_x1 = p_u1_x1;
id.p_u1_x2 = p_u1_x2;
id.p_u1_x3 = p_u1_x3;
id.p_u2_x1 = p_u2_x1;
id.p_u2_x2 = p_u2_x2;
id.p_u2_x3 = p_u2_x3;
id.p_u3_x1 = p_u3_x1;
id.p_u3_x2 = p_u3_x2;
id.p_u3_x3 = p_u3_x3;


}  // end of derivativesOfFields for 3D
//====================================================================
void ElasticityElement3D::printNodalStressesInLCS(ostream &stressout)
{
    stressout << elementNumber << " " << getMaterialNumber() << endl;
    extrapolateQuadValuesToNodes(bag->rotatedStresses,numStrains,bag->extrapolatedStresses);
    for(int i=0;i<numNodesPerElement;++i){
        for(int j=0;j<numStrains;++j){
            stressout << FF1 << bag->extrapolatedStresses[i][j];
        }
        stressout << "  " << node[i].nodeNum << endl;
    }
}
//====================================================================
void ElasticityElement3D::printNodalStressesInGCS(ostream &stressout)
{
    stressout << elementNumber << " " << getMaterialNumber() << endl;
    extrapolateQuadValuesToNodes(bag->quadPointStresses,numStrains,bag->extrapolatedStresses);
    for(int i=0;i<numNodesPerElement;++i){
        for(int j=0;j<numStrains;++j){
            stressout << FF1 << bag->extrapolatedStresses[i][j];
        }
        stressout << "  " << node[i].nodeNum << endl;
    }
}
//====================================================================
void ElasticityElement3D::printNodalStrainsInLCS(ostream &strainout)
{
    strainout << elementNumber << " " << getMaterialNumber() << endl;
    extrapolateQuadValuesToNodes(bag->rotatedStrains,numStrains,bag->extrapolatedStrains);
    for(int i=0;i<numNodesPerElement;++i){
        for(int j=0;j<numStrains;++j){
            strainout << FF1 << bag->extrapolatedStrains[i][j];
        }
        strainout << "  " << node[i].nodeNum << endl;
    }
}
//====================================================================
void ElasticityElement3D::printNodalStrainsInGCS(ostream &strainout)
{
    strainout << elementNumber << " " << getMaterialNumber() << endl;
    extrapolateQuadValuesToNodes(bag->quadPointStrains,numStrains,bag->extrapolatedStrains);
    for(int i=0;i<numNodesPerElement;++i){
        for(int j=0;j<numStrains;++j){
            strainout << FF1 << bag->extrapolatedStrains[i][j];
        }
        strainout << "  " << node[i].nodeNum << endl;
    }
}
//==================================================================
void ElasticityElement3D::printQuadStrains(int totalNumIPs, ostream &SO)  // in global coordinate system
//JV052304 - this function does not 'calculate' the stresses/strains - it only prints them
// the only thing it calculates is the vol at the quad points.
{
double **quadPointStresses   = bag->quadPointStresses;
double **rotatedStresses     = bag->rotatedStresses;
double **rotatedStrains      = bag->rotatedStrains;
double **quadPointStrains    = bag->quadPointStrains;

	GaussPointList	   gaussPointList;
    InterpAndSolDerivatives interpData(numNodesPerElement,numNodesPerElement); interpData.allocate();

    interpData.numberOfInterp = numNodesPerElement;
    extractNodalCoordinates(interpData.xCoor,interpData.yCoor,interpData.zCoor);
 
    getQuadraturePoints(gaussPointList); 

	double weight,vol;
	int group,ip,i;
    group=getMaterialNumber();		
    
	//SO << elementNumber << "  " << totalNumIPs << endl;
	for(ip=0;ip<totalNumIPs;ip++)
	{
		interpolation(ip, gaussPointList, &interpData);
		jacobian(&interpData);
		weight = gaussPointList.weightList[ip];
		vol = interpData.detJ * weight * getIntegrationFactor();
		for (i=0; i<numStrains; i++){
			SO << DOUBLE_FORMAT << quadPointStrains[ip][i] << " ";   // global coords.
		}
	  for (i=0; i<numStrains; i++)
	  { SO << DOUBLE_FORMAT <<rotatedStrains[ip][i] << " "; }  // local coords.
		SO << DOUBLE_FORMAT << vol << " ";
		SO << INT_FORMAT << group << " " << INT_FORMAT << elementNumber << endl;
	}
}
//==================================================================
//this is a function that prints an old type of format  for the quad stress file
/*
void ElasticityElement3D::printQuadStresses(int totalNumIPs, ostream &SO)  // in local coordinate system
{
//JV052304 - this function does not 'calculate' the stresses/strains - it only prints them
// the only thing it calculates is the vol at the quad points.

//	if(elementNumber==0)
//	{
//	SO<<"Format:"           <<endl;
//	SO<<"ElementNumber MaterialGroup NumQuadPoints ElementVolume"<<endl;
//	SO<<"For Each QuadPoint: QuadStresses QuadPointVolume MaterialGroup ElementNumber" <<endl;
//	}

double **quadPointStresses   = bag->quadPointStresses;
double **rotatedStresses     = bag->rotatedStresses;
double **rotatedStrains      = bag->rotatedStrains;
double **quadPointStrains    = bag->quadPointStrains;

	GaussPointList	   gaussPointList;
    InterpAndSolDerivatives interpData;

    interpData.numberOfInterp = numNodesPerElement;
    extractNodalCoordinates(interpData.xCoor,interpData.yCoor,interpData.zCoor);
 
    getQuadraturePoints(gaussPointList); 

	double weight,vol;
	int group,ip,i;
    group=getMaterialNumber();		
//    if(materialRotationAxis ==0) {stress=quadPointStresses;} 
  //  else stress=rotatedStresses; 

	for(ip=0;ip<totalNumIPs;ip++)
	{
	  interpolation(ip, gaussPointList, interpData);
      jacobian(interpData);
	  weight = gaussPointList.weightList[ip];
	  vol = interpData.detJ * weight * getIntegrationFactor();
      SO << setiosflags(ios::scientific);
	  SO << setprecision(10);
	  for (i=0; i<numStrains; i++)
	  { SO << F2 << quadPointStresses[ip][i]*vol  << " "; }  // global coords.
//	  for (i=0; i<numStrains; i++)
//	  { SO << F2 <<rotatedStrains[ip][i] << " "; }  // local coords.
	    SO << F2 << vol << " ";
	    SO << F1 << group << " " << elementNumber << endl;
	}
}
*/
//==================================================================
void ElasticityElement3D::printQuadStresses(int totalNumIPs, ostream &SO)
{
//JV052304 - this function does not 'calculate' the stresses/strains - it only prints them
// the only thing it calculates is the vol at the quad points.
	/*
	if(elementNumber==0)
	{
	SO<<"Format:"           <<endl;
//	SO<<"ElementNumber MaterialGroup NumQuadPoints ElementVolume"<<endl;
	SO<<"For Each QuadPoint: QuadStresses QuadPointVolume MaterialGroup ElementNumber" <<endl;
	}
	*/

double **quadPointStresses   = bag->quadPointStresses;
double **rotatedStresses     = bag->rotatedStresses;
double **rotatedStrains      = bag->rotatedStrains;
double **quadPointStrains    = bag->quadPointStrains;

    double x,y,z;
    double **stress;


	GaussPointList	   gaussPointList;
    InterpAndSolDerivatives interpData(numNodesPerElement,numNodesPerElement); interpData.allocate();

    interpData.numberOfInterp = numNodesPerElement;
    extractNodalCoordinates(interpData.xCoor,interpData.yCoor,interpData.zCoor);
 
    getQuadraturePoints(gaussPointList); 

	double weight,vol;
	int group,ip,i;
    group=getMaterialNumber();
    stress=rotatedStresses;
	for(ip=0;ip<totalNumIPs;ip++){
		interpolation(ip, gaussPointList, &interpData);
		jacobian(&interpData);

		x=y=z=0.0;
		for (i=0; i<numNodesPerElement; i++){
			x=x+interpData.S[i]*interpData.xCoor[i];
			y=y+interpData.S[i]*interpData.yCoor[i];
			z=z+interpData.S[i]*interpData.zCoor[i];
		}
		weight = gaussPointList.weightList[ip];
		vol = interpData.detJ * weight * getIntegrationFactor();
		for (i=0; i<numStrains; i++){
			SO << DOUBLE_FORMAT << quadPointStresses[ip][i] << " "; // global coords.
		}  
		for (i=0; i<numStrains; i++){
			SO << DOUBLE_FORMAT << rotatedStresses[ip][i] << " "; // local coords.
		}  
		SO << DOUBLE_FORMAT << vol << " ";
		SO << INT_FORMAT << group << " ";
		SO << DOUBLE_FORMAT << x << " " << y << " " << z << " " << endl;
	}
}
//====================================================================
void ElasticityElement3D::printNodalForcesByElement(ostream &forcesout)
{
    forcesout << elementNumber << " " << getMaterialNumber() << endl;
    int Index = 0;
    for(int i=0;i<numNodesPerElement;++i){
        for(int j=0;j<node[i].numDof;++j){
            forcesout << Fe[Index]  << '\t';
            ++Index;
        }
        forcesout << node[i].nodeNum << endl;
    }
}
//====================================================================
void ElasticityElement3D::getVolAvgValues(double *uGlobal, double *volAvgStress,
	         double *volAvgStrain, double *SED, double *volume)
{
// Calculate volume averaged values (strain and stress) and volume
// (Jae Noh, 6/23/98)
#define maxStresses 6
//InterpAndSolDerivatives dispDeriv; 
GaussPointList				gaussPointList;
InterpAndSolDerivatives       interpData(numNodesPerElement,numNodesPerElement); interpData.allocate();
int     i;
double *elemDisp[3];
Matrix *cMat, /**sMat,*/ *sMatOrtho, *sMatNonOrtho;
double strain[6];
//double SEDtmp0[6], SEDtmp1[6], SEDtmp2[6], SED2[3];
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
    //double MaterialAngle = getMaterialAngle(&interpData);

    material->calculateStress(stress,
                              strain,
                              isv[ip],
                              OrientationAtIP);
    OrientationAtIP->rotate_Voigt_stress_from_global_to_local(stress,rotStress);

   // calculate strain energy density
   //cMat->print("cmat");
		/* //commenting this entire block - replacing with code from xtelastic element for calculating SED
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
	*/ //end block

	for(i=0; i<6; i++) SED[0] += 0.5 * stress[i]*strain[i] * integrationFactor;

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
void ElasticityElement3D::addToKe(Matrix* cMat,InterpAndSolDerivatives *id_ptr,
									  double integrationFactor,double *stress)       
{
	InterpAndSolDerivatives &id=*id_ptr;

	double *S=id.S,*p_S_x1=id.p_S_x1,    *p_S_x2=id.p_S_x2,     *p_S_x3=id.p_S_x3;	
	int    ii,jj;

   double C_x_B[6][96];

	double * b0 = &b[0][0];
	double * b1 = &b[1][0];
	double * b2 = &b[2][0];
	double * b3 = &b[3][0];
	double * b4 = &b[4][0];
	double * b5 = &b[5][0];

	for(ii=0; ii<numDof; ii++)
		{double b0ii = b[0][ii];
		 double b1ii = b[1][ii];
		 double b2ii = b[2][ii];
		 double b3ii = b[3][ii];
		 double b4ii = b[4][ii];
		 double b5ii = b[5][ii];
		 for(jj=0; jj<numStrains;jj++)
			  {C_x_B[jj][ii] = b0ii * (*cMat)[0][jj] +
									 b1ii * (*cMat)[1][jj] +
									 b2ii * (*cMat)[2][jj] +
									 b3ii * (*cMat)[3][jj] +
									 b4ii * (*cMat)[4][jj] +
									 b5ii * (*cMat)[5][jj] ;
			  }
		}

		for(ii=0; ii<numDof; ii++)
			{double p0ii = C_x_B[0][ii];
			 double p1ii = C_x_B[1][ii];
			 double p2ii = C_x_B[2][ii];
			 double p3ii = C_x_B[3][ii];
			 double p4ii = C_x_B[4][ii];
			 double p5ii = C_x_B[5][ii];
			 for(jj=ii; jj<numDof; jj++)
				 {
					(*Ke)[ii][jj] += ( p0ii * b0[jj] +
									       p1ii * b1[jj] +
											 p2ii * b2[jj] +
											 p3ii * b3[jj] +
											 p4ii * b4[jj] +
											 p5ii * b5[jj]  ) * integrationFactor;
				 }
			 }


if(analysisType==GEOMETRIC_NONLINEAR)    //add Ksigma............check this code!!
 {
 double b11, b11dj; 
 int mb1,nb1, nodei,nodej;
 
 #define ndofpn  3
 int numberOfInterp = numDof/3;
 for(nodei= 0;  nodei<numberOfInterp; nodei++)
    {
	  mb1= nodei * ndofpn  ; 
     for(nodej= nodei;  nodej<numberOfInterp; nodej++)
        {
         b11 =stress[0] * p_S_x1[nodei]* p_S_x1[nodej] +
			  stress[1] * p_S_x2[nodei]* p_S_x2[nodej] + 
			  stress[2] * p_S_x3[nodei]* p_S_x3[nodej] + 
			  stress[3] *(p_S_x1[nodei]* p_S_x2[nodej] + p_S_x2[nodei]*p_S_x1[nodej]) +
			  stress[4] *(p_S_x2[nodei]* p_S_x3[nodej] + p_S_x3[nodei]*p_S_x2[nodej]) +
			  stress[5] *(p_S_x1[nodei]* p_S_x3[nodej] + p_S_x3[nodei]*p_S_x1[nodej]) ;

                nb1= nodej * ndofpn  ;
                b11dj = b11 * integrationFactor       ;
                (*Ke)[mb1  ][nb1  ] +=  b11dj ;
                (*Ke)[mb1+1][nb1+1] +=  b11dj ;
				(*Ke)[mb1+2][nb1+2] +=  b11dj ;
			}		 
    }	 	 
	}  // end of Ksigma
#undef ndofpn
}  //end of addToKe
//======================================================================
void ElasticityElement3D::initializeArrays()
{
	for(int i=0;i<numDof;i++) {
		Fe[i] = 0.;
		initialFe[i] = 0.;
		for(int j=i;j<numDof;j++)
			(*Ke)[i][j] = 0.;
	}
}
//====================================================================
bool ElasticityElement3D::geometryCheck()
{
	GaussPointList			gaussPointList;
	InterpAndSolDerivatives       id(numNodesPerElement,numNodesPerElement); id.allocate();
	id.numberOfInterp = numNodesPerElement;
	extractNodalCoordinates(id.xCoor,id.yCoor,
	                      id.zCoor);
	getQuadraturePoints(gaussPointList); 
	int totalNumIPs = gaussPointList.totalNumIPs;
	bool GeoOK=true;
	for(int ip=0;ip<totalNumIPs;ip++)
	{
		interpolation(ip, gaussPointList, &id);
		jacobian(&id);
		if (id.detJ < 0) GeoOK=false;
	}  // end of integration loop
	return GeoOK;
}
//=====================================================================
void ElasticityElement3D::setPointers(ElementWorkspace * s)
{
	b                    = s->b;
	bag                  = s;
	BasicElement::setPointers(s);
}
//=====================================================================
void ElasticityElement3D::initializeSummary()
{
	bag->totalVolume = 0.;
	for(int i=0; i<6; i++) {
		bag->strainVolume_total[i] = 0.; 
		bag->stressVolume_total[i] = 0.;
		bag->localStressVolume_total[i] =0.;  
		bag->energyVolume_total[i] =0; 
	}
}
//=====================================================================
void ElasticityElement3D::addToSummary()
{
}
//====================================================================
void ElasticityElement3D::outputSummary()
{
//double &volume = bag->volume;
double *strainVolume_total = bag->strainVolume_total;
double *stressVolume_total = bag->stressVolume_total;
double *energyVolume_total = bag->energyVolume_total;
double *localStressVolume_total = bag->localStressVolume_total;



double totalVolume = bag->totalVolume;
#define F1 setw(5)
//#define F2 setw(13)<<setprecision(4)
int j;

BETA_OUT<<"The current calculation assumes there are no cracks or holes in the body\n";
BETA_OUT<<"... if there are, the volume averages are incorrect!\n";

BETA_OUT<<"Total volume of model="<<totalVolume <<endl;
//cout<<"Total volume of model="<<bag->strainVolume_total <<endl; exit(1);

BETA_OUT<<"Volume averaged strains"<<endl;
for(j=0; j<6; j++) 
{BETA_OUT<<DOUBLE_FORMAT << strainVolume_total[j]/totalVolume<<"  ";}
BETA_OUT<<endl;
  
BETA_OUT<<"Volume averaged stresses"<<endl;
for(j=0; j<6; j++) 
{BETA_OUT<<DOUBLE_FORMAT << stressVolume_total[j]/totalVolume<<"  ";}
BETA_OUT<<endl;
 
		  //BETA_OUT<<"Volume averaged localStress"<<endl;
          //for(j=0; j<6; j++) 
          //{BETA_OUT<<F2<<localStressVolume_total[j]/totalVolume<<"  ";}
          //BETA_OUT<<endl;
 

		  BETA_OUT<<"Components of strain energy density"<<endl;
          for(j=0; j<6; j++) 
          {BETA_OUT<<DOUBLE_FORMAT << energyVolume_total[j]/totalVolume<<"  ";}
          BETA_OUT<<endl;

		 double strainEnergyDensity=0;
		 for(j=0; j<6; j++) strainEnergyDensity += energyVolume_total[j]/totalVolume;
		  BETA_OUT<<"Strain energy density = "<<DOUBLE_FORMAT << strainEnergyDensity<<endl;

	  BETA_OUT<<"Effective Exx= "<<DOUBLE_FORMAT << stressVolume_total[0]/strainVolume_total[0]<<endl;

/* where does this go ???? FIX !!!!
if(effectiveModulusCase > -1){
		//The parameter "ios::app" causes output to be appended
//		BETA_OUTPROP= new ofstream(name, ios::app);
		REPORT<<"weaveDescription="<< weaveDescription<<endl;
	    REPORT<<"effectiveModulusCase  "<<effectiveModulusCase <<endl;
//If we also output information about the weave type, #elements, etc., the
// data could be "harvested" for plotting, etc
	}

if(effectiveModulusCase==1){
	REPORT<<"          E11 = "<< stressVolume_total[0]/strainVolume_total[0]<<endl;
	REPORT<<"          NU12= "<<-strainVolume_total[1]/strainVolume_total[0]<<endl;
	REPORT<<"          NU13= "<<-strainVolume_total[2]/strainVolume_total[0]<<endl;
	double eps11 = strainVolume_total[0]/totalVolume;
	REPORT<<"          E11 = "<<2.*strainEnergyDensity/(eps11*eps11)<<endl;
} // 1

if(effectiveModulusCase==2){
	REPORT<<"          E22 = "<< stressVolume_total[1]/strainVolume_total[1]<<endl;
	REPORT<<"          NU21= "<<-strainVolume_total[0]/strainVolume_total[1]<<endl;
	REPORT<<"          NU23= "<<-strainVolume_total[2]/strainVolume_total[1]<<endl;
} // 2

if(effectiveModulusCase==3){
	REPORT<<"          E33 = "<< stressVolume_total[2]/strainVolume_total[2]<<endl;
	REPORT<<"          NU31= "<<-strainVolume_total[0]/strainVolume_total[2]<<endl;
	REPORT<<"          NU32= "<<-strainVolume_total[1]/strainVolume_total[2]<<endl;
} // 3

if(effectiveModulusCase==4){
	REPORT<<"          G12 = "<< stressVolume_total[3]/strainVolume_total[3]<<endl;
} // 4

if(effectiveModulusCase==5){
	REPORT<<"          G23 = "<< stressVolume_total[4]/strainVolume_total[4]<<endl;
} // 5

if(effectiveModulusCase==6){
	REPORT<<"          G13 = "<< stressVolume_total[5]/strainVolume_total[5]<<endl;
} // 6
*/

}
//==========================================================
void ElasticityElement3D::readSpecialCommand(istream &inStream, ElementGroup *element,
								    int numElements, char * command)
{

BETA_OUT<<"ElasticityElement3D::readSpecialCommand"<<endl;

//================================================================
if (COMPARE(command, "setEffectiveModulusCase")==0 ) {
	inStream>> effectiveModulusCase;
	inStream>> weaveDescription;
	BETA_OUT<<"EFFECTIVE_MODULUS_CASE = "<< effectiveModulusCase<<endl;
	BETA_OUT<<"Weave Description = "<< weaveDescription<<endl;
	return;
}; //setEffectiveModulusCase



// If no match, try super class
BETA_OUT<<"No match... try super class"<<endl;
IsoElement::readSpecialCommand(inStream,element,numElements, command);
}//readSpecialCommand
//==========================================================
void ElasticityElement3D::setTemperature(double temperature)
{
    int numIPs=getTotalNumberOfIPs();
    //If isv isn't allocated, allocate it
    if(!checkForAllocatedISVs()){
        material->allocateISVs(isv,numIPs);
    }
    for(int i=0;i<numIPs;++i){
        //Typecast - need IsA function
        ((ElasticISV*) isv[i])->setTemperature(temperature);
    }
}
//==========================================================
void ElasticityElement3D::setMoisture(double moisture)
{
    int numIPs=getTotalNumberOfIPs();
    //If isv isn't allocated, allocate it
    if(!checkForAllocatedISVs()){
        material->allocateISVs(isv,numIPs);
    }
    for(int i=0;i<numIPs;++i){
        //Typecast - need IsA function
        ((ElasticISV*) isv[i])->setMoisture(moisture);
    }
} 
//==========================================================
/*
void ElasticityElement3D::openOutputFiles(FileManager *fm, int flags)
{
	if(flags == 0){
		BETA_OUT << "No Optional Output Parameters found!" << endl;
		return;
	}
	if(flags & stress){
		stressout=fm->OpenOutputStream("stress");
		*stressout << "stress" << endl;
	}
	if(flags & strain){
		strainout=fm->OpenOutputStream("strain");
		*strainout << "stress" << endl;

	}
	if(flags & quadstress){
		quadstressout=fm->OpenOutputStream("quadstress");
	}
	if(flags & quadstrain){
		quadstrainout=fm->OpenOutputStream("quadstrain");
	}
	if(flags & volumeAverage){
		volumeAverageout=fm->OpenOutputStream("volumeAverage");
	}
}
//=========================================================================
void ElasticityElement3D::closeOutputFiles(FileManager *fm, int flags)
{
	if(flags & stress){
		fm->CloseOutputStream((ofstream*)stressout);
	}
	if(flags & strain){
		fm->CloseOutputStream((ofstream*)strainout);
	}
	if(flags & quadstress){
		fm->CloseOutputStream((ofstream*)quadstressout);
	}
	if(flags & quadstrain){
		fm->CloseOutputStream((ofstream*)quadstrainout);
	}
	if(flags & volumeAverage){
		fm->CloseOutputStream((ofstream*)volumeAverageout);
	}
	if(flags == 0){
		BETA_OUT << "No Optional Output Parameters found!" << endl;
	}
}
*/
//==================================================================
void ElasticityElement3D::InitializeOutputSummary()
{
	if(!ElasticityElement3D::initialized){
		ElasticityElement3D::elematfile=new ofstream("~ElasticityElement3D.elemat");
		//*ElasticityElement3D::elematfile << "readmaterialgroup()" << endl;

		ElasticityElement3D::manglesfile=new ofstream("~ElasticityElement3D.mangles");
		//*ElasticityElement3D::manglesfile << "setElementMaterialAngles" << endl;

		ElasticityElement3D::initialized=true;
	}
}
//==================================================================
void ElasticityElement3D::OutputSummary()
{
	if(!activeElementFlag) return;
	//elemat file
	*ElasticityElement3D::elematfile << newElementNumber << " " << newElementNumber << " 1 " << materialGroup << endl;

	////mangles file
	//if(materialRotationAxis!=0){
	//	*ElasticityElement3D::manglesfile << newElementNumber << " ";
	//	*ElasticityElement3D::manglesfile << materialRotationAxis << " ";
	//	for(int iii=0; iii<numNodesPerElement; iii++){
	//			*ElasticityElement3D::manglesfile << materialRotationAnglesAtNodes[iii] << " ";
	//	}
	//	*ElasticityElement3D::manglesfile << endl;
	//}
}
//==================================================================
void ElasticityElement3D::FinalizeSummary()
{
	if(ElasticityElement3D::initialized){
		//*ElasticityElement3D::elematfile << "-1 0 0" << endl;
		((ofstream*)ElasticityElement3D::elematfile)->close();

		((ofstream*)ElasticityElement3D::manglesfile)->close();
		ElasticityElement3D::initialized=false;
	}
}
//==================================================================
void ElasticityElement3D::ConsolidateSummary(string outputtype, ofstream &ofile)
{
	if(ElasticityElement3D::consolidated) return;

	ifstream is;
	ChangeToUpper(outputtype)

	if(outputtype=="ELEMAT"){
		is.open("~ElasticityElement3D.elemat");
	}
	else if(outputtype=="MANGLES"){
		is.open("~ElasticityElement3D.mangles");
	}else{
		return;
	}

	char line[255];
	while(!is.eof()){
		is.getline(line, 255);
		ofile << line << endl;
	}

	is.close();
	ElasticityElement3D::consolidated=true;
}
//==================================================================
void ElasticityElement3D::setConsolidate(bool flag)
{
	ElasticityElement3D::consolidated=flag;
}
//===========================================================================
void ElasticityElement3D::processActionCode(char * action, IntegrationDataType &idt)
{
	int calculateKFlag=0, 
		calculateForceFlag=0, 
		calculateInitialForceFlag=0, 
		calculateStressFlag=0;

	if(strstr(action,"K") !=NULL) calculateKFlag             = 1;
	if(strstr(action,"F") !=NULL) calculateForceFlag         = 1;
	if(strstr(action,"I") !=NULL) calculateInitialForceFlag  = 1;
	if(strstr(action,"S") !=NULL) calculateStressFlag        = 1;

	if(analysisType==GEOMETRIC_NONLINEAR || calculateForceFlag==1)
		{calculateStressFlag=1;}

	if(calculateForceFlag==1 && calculateInitialForceFlag==1){
		EXIT_BETA("Cannot calculate current element forces and initial forces during the same call");
    }

	calculateElementVolume(idt);
	if (calculateStressFlag==1)			calculateStress(idt); //this needs to be the first option (before calculateFe and calculateKe)
	if (calculateForceFlag==1)			calculateFe(idt);
	if (calculateInitialForceFlag==1)	calculateInitialForce(idt);
	if (calculateKFlag==1)				calculateKe(idt);
}
//===========================================================================
void ElasticityElement3D::loopOverIntegrationPoints(char * action, IntegrationDataType &idt)
{
	initializeArrays();
	bag->elementVolume=0.;

	double *strainVolume_element = bag->strainVolume_element;
	double *stressVolume_element = bag->stressVolume_element;
	double *energyVolume_element = bag->energyVolume_element;
	double *localStressVolume_element = bag->localStressVolume_element;

	for(int i=0; i<numStrains; i++){
		strainVolume_element[i]		=0.;
		stressVolume_element[i]		=0.;
		energyVolume_element[i]		=0.;
		localStressVolume_element[i]=0.;
	}

	ElasticMaterial* material = (ElasticMaterial*)BasicElement::material;

	int totalNumIPs = idt.gaussPointList.totalNumIPs;
	double weight;
	for(int ip=0;ip<totalNumIPs;ip++) {
		// set up data or the current point: 
		// these data must be obtained for all routines and 
		// are basically for the integration factor
		interpolation(ip, idt.gaussPointList, idt.id);
		jacobian(idt.id);
		weight = idt.gaussPointList.weightList[ip];
		idt.integrationFactor = idt.id->detJ * weight * getIntegrationFactor();
		idt.ip=ip;
        if(idt.uGlobal!=NULL){derivativesOfFields(idt.id);}
        calculateBmatrix(idt.id,b);
        //OrientationAtIP->interpolate_rotation(idt.id,ElementNodeRotations);
        OrientationAtIP->interpolate_rotation_by_angle(idt.id,ElementNodeRotations);
        //OrientationAtIP->interpolate_old_beta(idt.id,ElementNodeRotations,material->getMaterialRotation());
		material->rotateMaterialToGlobal(OrientationAtIP);
		processActionCode(action, idt);
	}

	bag->totalVolume += bag->elementVolume;
	for(int i=0; i<numStrains; i++){
		bag->strainVolume_total[i]     += strainVolume_element[i];
		bag->stressVolume_total[i]     += stressVolume_element[i];
		bag->energyVolume_total[i]     += energyVolume_element[i];
		bag->localStressVolume_total[i]+= localStressVolume_element[i];
	}
	if(strstr(action,"K") !=NULL) checkForNegativeDiagonal(); 
}
//===========================================================================
void ElasticityElement3D::update(char * action, double *uGlobal)
{
	IntegrationDataType   idt;
	idt.id=new InterpAndSolDerivatives(numNodesPerElement,numNodesPerElement); idt.id->allocate();

	idt.elemDisp[0] = idt.id->nodalValues_u1; 
	idt.elemDisp[1] = idt.id->nodalValues_u2;
	idt.elemDisp[2] = idt.id->nodalValues_u3;     
	idt.id->numberOfInterp = numNodesPerElement;
	idt.uGlobal=uGlobal;

	extractNodalCoordinates(idt.id->xCoor,idt.id->yCoor,idt.id->zCoor);
	if( uGlobal != NULL) {extractSolution( uGlobal, idt.elemDisp);} 
	getQuadraturePoints(idt.gaussPointList); 

    //Will be allocated if mangles has been called.
    if(ElementNodeRotations == 0){AllocateRotations();}

	loopOverIntegrationPoints(action, idt);
}
//=========================================================================
void ElasticityElement3D::calculateFe(IntegrationDataType &idt)
{
	double **quadPointStresses   = bag->quadPointStresses;
	double *stress	= quadPointStresses[idt.ip];

	// JV051906 - forcing only one stress component to be non zero
	if(ForceSingleNonZeroStressComponent != -1){
		for(int i=0; i<numStrains; i++)
			if(i!=ForceSingleNonZeroStressComponent) stress[i]=0;
	}

	addToFe(numDof,numStrains,b,stress,idt.integrationFactor,Fe);
}
//===========================================================================
void ElasticityElement3D::calculateInitialForce(IntegrationDataType &idt)
{
	double *stress	= bag->quadPointStresses[idt.ip];

	ElasticMaterial* material=(ElasticMaterial*)BasicElement::material;
	material->calculateInitialStress(stress, isv[idt.ip], OrientationAtIP);
	addToFe(numDof,numStrains,b,stress,idt.integrationFactor,initialFe);
}
//===========================================================================
void ElasticityElement3D::calculateKe(IntegrationDataType &idt)
{
	// local properties that require transformations	   
	ElasticMaterial* material=(ElasticMaterial*)BasicElement::material;
	Matrix *cMat = material->getPointerToGlobalCmat();
		
	double *stress	= bag->quadPointStresses[idt.ip];
	addToKe(cMat,idt.id,idt.integrationFactor,stress);
}
//===========================================================================
void ElasticityElement3D::calculateStress(IntegrationDataType &idt)
{
	// local properties that require transformations	   
	ElasticMaterial* material=(ElasticMaterial*)BasicElement::material;

	double	*strain		= bag->quadPointStrains[idt.ip];
	double	*stress		= bag->quadPointStresses[idt.ip];
	double	*rotStress	= bag->rotatedStresses[idt.ip];
	double	*rotStrain	= bag->rotatedStrains[idt.ip];
    
	calculateStrain(idt.id, strain); //in global coordinate system
	material->calculateStress(stress, 
                              strain, 
                              isv[idt.ip],
                              OrientationAtIP);
    OrientationAtIP->rotate_Voigt_stress_from_global_to_local(stress,rotStress);
    OrientationAtIP->rotate_Voigt_engineering_strain_from_global_to_local(strain,rotStrain);

	double *strainVolume_element		= bag->strainVolume_element;
	double *stressVolume_element		= bag->stressVolume_element;
	double *energyVolume_element		= bag->energyVolume_element;
	double *localStressVolume_element	= bag->localStressVolume_element;

	for(int i=0; i<numStrains; i++){
		strainVolume_element[i]		+= strain[i] *idt.integrationFactor;
		stressVolume_element[i]		+= stress[i] *idt.integrationFactor;
		energyVolume_element[i]		+= .5 * stress[i] * strain[i] *idt.integrationFactor;
		localStressVolume_element[i]+= rotStress[i] *idt.integrationFactor;
	}
}
//===========================================================================
//                             Static Data                        
//===========================================================================
int			ElasticityElement3D::effectiveModulusCase;
char		ElasticityElement3D::weaveDescription[80];
ofstream*	ElasticityElement3D::elematfile=0;
ofstream*	ElasticityElement3D::manglesfile=0;
bool		ElasticityElement3D::initialized=false;
bool		ElasticityElement3D::consolidated=false;
