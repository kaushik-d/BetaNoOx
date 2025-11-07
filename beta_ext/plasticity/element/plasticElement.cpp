#include "stdafx.h"

#include <fstream>
#include <string.h>
#include <iostream>
#include <iomanip>
#include "math/gaussQuadrature/hexQuadPoint.hpp"
#include "elements/InterpAndSolDerivatives.hpp"

#include "plasticElement.hpp"

void deriv3D(int numberOfInterp,
              double xi,  double eta, double zeta,				  
				  double * Sxi, double * Seta, double * Szeta);

void  shape3D(int numberOfInterp,
              double xi,  double eta, double zeta,				  
				  double * S);   //  this function is defined in lib: "elementLib"

// +xtang 981104
/*void calculate_extrSTRAIN(double **quadPointStresses,
	int numNodesPerElement, int numStrains , int integrationOrder,
        int totalNumIPs ); */  
// -xtang 981104   

extern int verboseFlag;

//====================================================================
PlasticElement::PlasticElement()
{
	isv=0;
}
//====================================================================
PlasticElement::~PlasticElement()
{
	((PlasticMaterial *) material)->deleteISVSpace(isv);
}
//====================================================================
void PlasticElement::
                      update(char * action, double *uGlobal)
//====================================================================
{
 
	InterpAndSolDerivatives       id(numNodesPerElement,numNodesPerElement); id.allocate();
//	DisplacementDerivatives dispDeriv; 
//	InterpolationData	interpData;
	GaussPointList		gaussPointList;
	int		calculateStressFlag,
			calculateInitialForceFlag,
			calculateKFlag, 
			calculateForceFlag,
			calculateISVFlag,
			outQuadInfoFlag,
	    updateGaussISVFlag,
	checkStrengthRatioFlag;

//	int     printStressFlag;
//	int     printStrainFlag; // xtang 981104
	double	*elemDisp[3];
	Matrix	*cMat;
	double	strain[10],matStrain[10];
	double	*stress;
	double	*rotStress;

	double **quadPointStresses   = bag->quadPointStresses;
	double **rotatedStresses     = bag->rotatedStresses;
	//double **extrapolatedStresses= bag->extrapolatedStresses;

	//double **rotatedStrains      = bag->rotatedStrains;
	double **quadPointStrains    = bag->quadPointStrains;

    //double x,y,z, vol;  // 990717

	//
	elemDisp[0] = id.nodalValues_u1;//dispDeriv.xDisp; 
	elemDisp[1] = id.nodalValues_u2;//dispDeriv.yDisp;
	elemDisp[2] = id.nodalValues_u3;//dispDeriv.zDisp;     

	id.numberOfInterp = numNodesPerElement;

//	printStressFlag = getPrintStressFlag();
//	printStrainFlag = getPrintStrainFlag();  // xtang 981104
 
	extractNodalCoordinates(id.xCoor,id.yCoor,id.zCoor);
 
	if( uGlobal != NULL) {extractSolution( uGlobal, elemDisp);} 
 
	calculateForceFlag = calculateInitialForceFlag=0;
	calculateKFlag     = calculateStressFlag=0;
	calculateInitialForceFlag  = 0;
	calculateStressFlag        = 0;
 	calculateISVFlag        = 0;
 	outQuadInfoFlag        = 0;
	updateGaussISVFlag        = 0;
	checkStrengthRatioFlag    =0;  // 991016 xtang
 
	if(strstr(action,"K") !=NULL) calculateKFlag             = 1;
	if(strstr(action,"F") !=NULL) calculateForceFlag         = 1;
	if(strstr(action,"I") !=NULL) calculateInitialForceFlag  = 1;
	if(strstr(action,"S") !=NULL) calculateStressFlag        = 1;
	if(strstr(action,"V") !=NULL) calculateISVFlag           = 1;
	if(strstr(action,"U") !=NULL) updateGaussISVFlag         =  1;
	if(strstr(action,"Q") !=NULL) outQuadInfoFlag            = 1;
	if(strstr(action,"R") !=NULL) checkStrengthRatioFlag     = 1;  // check maximum strength ratio and store in XTCommon

	if(strstr(action,"A") !=NULL) { // xtang 052801: set up integration factor and material angle for each Gauss point
		setIntegrationFactorInISV();
		return;
	}


	if(analysisType==GEOMETRIC_NONLINEAR || calculateForceFlag==1)
		{calculateStressFlag=1;}

	if(calculateForceFlag==1 && calculateInitialForceFlag==1)
		{ BETA_OUT<< "Cannot calculate current element forces and"<<endl;
		BETA_OUT<< "initial forces during the same call"<<endl;
		exit(1);
    }

	initializeArrays();
	getQuadraturePoints(gaussPointList); 
	int totalNumIPs = gaussPointList.totalNumIPs;

	double **qStrains;
	qStrains = new double * [totalNumIPs];
	
	/////////// cout << " Test in" << elementNumber << endl;
	// in consistent with Strain[10] defined at the begining
	int i;

	for(i=0; i<totalNumIPs; i++) qStrains[i] = new double [10];  
	// -xtang  -ending

	double integrationFactor;
	double volume =0.;
	double weight;
	int group=getMaterialNumber();		

//double strengthRatio;

	for(int ip=0;ip<totalNumIPs;ip++)  {

		interpolation(ip, gaussPointList, &id);
        jacobian(&id);
		weight = gaussPointList.weightList[ip];
		integrationFactor = id.detJ * weight * getIntegrationFactor();
        volume += integrationFactor;

             
        //if(uGlobal!=NULL){dispGradients(id,dispDeriv);}               
		if(uGlobal!=NULL){derivativesOfFields(&id);}
        //calculateBmatrix(interpData,dispDeriv,b);
		calculateBmatrix(&id,b);
	        
		stress      = quadPointStresses[ip];
		rotStress   = rotatedStresses[ip];
		
		// dStrain
		calculateStrain(&id,strain); // xtang 98-10-28, in global coord system
		
		for (int xt=0; xt<numStrains; xt++) qStrains[ip][xt]=strain[xt];  // xtang 98-10-28

        // +xtang 990420
  
		PlasticMaterial* material = (PlasticMaterial*)BasicElement::material;

        //OrientationAtIP->interpolate_rotation(&id,ElementNodeRotations);
        OrientationAtIP->interpolate_rotation_by_angle(&id,ElementNodeRotations);
        //OrientationAtIP->interpolate_old_beta(&id,ElementNodeRotations,material->getMaterialRotation());

        OrientationAtIP->rotate_Voigt_engineering_strain_from_global_to_local(strain,matStrain);
   //     material->rotateEngineeringStrain(matStrain,strain, -getMaterialAngle(&id),
			//materialRotationAxis, true);
		
		// 062001
        if(calculateISVFlag==1){ // "V"
			// get dStress from dStrain 
			((PlasticMaterial*) material)->updateISV(&isv[ip], matStrain); 
		}

        if(updateGaussISVFlag==1){  //"U"
			// update dStress and dState
			isv[ip].updateISV();
		}

        if(calculateForceFlag==1){ // "F"
			// get stress from ISV and tranform it to global coordinate system
            OrientationAtIP->rotate_Voigt_stress_from_local_to_global(isv[ip].dStress,stress);
			//material->rotateStress(rotStress,isv[ip].dStress, getMaterialAngle(&id),
			//    materialRotationAxis, false);
		    addToFe(numDof,numStrains,b,stress,integrationFactor,Fe);
		}

		if(calculateInitialForceFlag==1) { // "I"
			// don't know how to handle this 
			// xtang: do nothing
			// material->calculateInitialStress(stress, stateVariable);
			// addToFe(numDof,numStrains,b,stress,integrationFactor,initialFe);
        }

		if(calculateKFlag==1) { // "K"
			// +xtang 990421:  update Material properties here
			((PlasticMaterial *) material)->updateLocalizedCmat(OrientationAtIP, isv[ip]);  // xt: 102300 elastic material
			cMat = material->getPointerToGlobalCmat();  // xt 102300: elastic material
			addToKe(cMat,&id,integrationFactor,stress);
		}

	}  // end of integration loop


    if(calculateInitialForceFlag==1  )  {
    }

    if(calculateKFlag==1) { 
//		checkForNegativeDiagonal(); 
    }

	
	//if(printStressFlag==1){ //JV071107
	    //extrapolateStresses(totalNumIPs);//, quadPointStresses);
	//}

	// +xtang 981104
//	if(printStrainFlag==1){
		// calculate_extrSTRAIN(qStrains,numNodesPerElement,numStrains , 
		// integrationOrder, totalNumIPs );
//	}
    // -xtang 981104

	// xtang 98-10-11
//	if(printStressFlag==1){ 
		//printQuadStresses(totalNumIPs);
//	}
    // -xtang 98-10-11

	// xtang 98-10-28
//	if(printStrainFlag==1){ 
		//printQuadStrains(totalNumIPs, qStrains);
//	}
    // -xtang 98-10-28

	// +xtang 98-10-28 
	// release memory
	for(i=0; i<totalNumIPs; i++) delete [] qStrains[i]; 
    delete [] qStrains; 
	// -xtang 98-10-28

return;
}

//====================================================================
void PlasticElement::printElementalStresses(ostream * outStream) //LCS
//this function does the following:
//1. take the quad stresses from the ISV (stress in LCS) and stores it in the element workspace
//2. extrapolates the quad stresses to the nodes and prints it to the file
{
	double **quadPointStresses   = bag->quadPointStresses;
    int totalNumIPs = getTotalNumberOfIPs();

	for(int ip=0;ip<totalNumIPs;ip++)  {
		for (int i=0; i<numStrains;i++) {
			quadPointStresses[ip][i]=isv[ip].stress[i]; //LCS
		}
		quadPointStresses[ip][numStrains]=isv[ip].effectiveStress;
		quadPointStresses[ip][numStrains+1]=isv[ip].effectivePlasticStrain;
	}
    int NumPlasticStresses = numStrains+2;
	
    double ** ExtrapolatedStresses = new double*[numNodesPerElement];
    double * ExtrapStressesPtr = new double[numNodesPerElement*NumPlasticStresses];
    for(int i=0;i<numNodesPerElement;++i){
        ExtrapolatedStresses[i] = &ExtrapStressesPtr[i*NumPlasticStresses];
    }
    *outStream << elementNumber << " " << getMaterialNumber() << endl;
    extrapolateQuadValuesToNodes(quadPointStresses,NumPlasticStresses,ExtrapolatedStresses);
    for(int i=0;i<numNodesPerElement;++i){
        for(int j=0;j<NumPlasticStresses;++j){
            *outStream << ExtrapolatedStresses[i][j] << " ";
        }
        *outStream << node[i].nodeNum << endl;
    }
    delete [] ExtrapStressesPtr;
    delete [] ExtrapolatedStresses;

}
//==================================================================
void PlasticElement::printQuadStresses(int totalNumIPs, ostream &SO)
{
//JV052304 - this function does not 'calculate' the stresses/strains - 
	//it only rotates and prints them
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
    //if(materialRotationAxis ==0) {stress=quadPointStresses;} 
    //else 
    stress=rotatedStresses; 
	
	for(ip=0;ip<totalNumIPs;ip++)
	{
		interpolation(ip, gaussPointList, &interpData);
		jacobian(&interpData);

		for (i=0; i<numStrains;i++) {
			rotatedStresses[ip][i]=isv[ip].stress[i]; //LCS
		}
        //OrientationAtIP->interpolate_rotation(&interpData,ElementNodeRotations);
        OrientationAtIP->interpolate_rotation_by_angle(&interpData,ElementNodeRotations);
        //OrientationAtIP->interpolate_old_beta(&interpData,ElementNodeRotations,material->getMaterialRotation());
        OrientationAtIP->rotate_Voigt_stress_from_local_to_global(rotatedStresses[ip],quadPointStresses[ip]);
		PlasticMaterial* material = (PlasticMaterial*)BasicElement::material;
		//material->rotateStress(quadPointStresses[ip],rotatedStresses[ip], getMaterialAngle(&interpData),
		//	materialRotationAxis, false);


	  x=y=z=0.0;
	  for (i=0; i<numNodesPerElement; i++)
	  {
		  x=x+interpData.S[i]*interpData.xCoor[i];
		  y=y+interpData.S[i]*interpData.yCoor[i];
		  z=z+interpData.S[i]*interpData.zCoor[i];
	  }
	  weight = gaussPointList.weightList[ip];
	  vol = interpData.detJ * weight * getIntegrationFactor();
      SO << setiosflags(ios::scientific);
	  SO << setprecision(16); //full precision
	  for (i=0; i<numStrains; i++)
	  { SO << quadPointStresses[ip][i] << " "; }  // global coords.
	  for (i=0; i<numStrains; i++)
	  { SO << rotatedStresses[ip][i] << " "; }  // local coords.
	    SO << vol << " ";
	    SO << group << " ";
		SO << x << " " << y << " " << z << " " << endl;
	}
}
//====================================================================
void PlasticElement::getVolAvgValues(double *uGlobal, 
									 double *volAvgStress,  // output
									 double *volAvgStrain,  // output
									 double *SED,   // output as plastic strain: 052801
									 double *volume)
//====================================================================

{
// Calculate volume averaged values (strain and stress) and volume

#define maxStresses 6
	GaussPointList				gaussPointList;
//	InterpolationData			interpData;
	InterpAndSolDerivatives       id(numNodesPerElement,numNodesPerElement); id.allocate();
	//int		i, j;
	double strain[maxStresses];
	double stress[maxStresses];
	//double plasticStrain[maxStresses];

	id.numberOfInterp = numNodesPerElement;
	extractNodalCoordinates(id.xCoor,id.yCoor,id.zCoor);
 
	initializeArrays();
	getQuadraturePoints(gaussPointList); 
	int totalNumIPs = gaussPointList.totalNumIPs;

	double integrationFactor, weight;

    *volume=0;  //xtang 051302
	int is;
	
	for(is=0; is<maxStresses; is++)
	{   volAvgStress[is] = volAvgStrain[is]=0.; }

	for(int ip=0;ip<totalNumIPs;ip++)
	{
		interpolation(ip, gaussPointList, &id);
		jacobian(&id);
		weight = gaussPointList.weightList[ip];
		integrationFactor = id.detJ * weight * getIntegrationFactor();
		*volume += integrationFactor;

        //OrientationAtIP->interpolate_rotation(&id,ElementNodeRotations);
        OrientationAtIP->interpolate_rotation_by_angle(&id,ElementNodeRotations);
        //OrientationAtIP->interpolate_old_beta(&id,ElementNodeRotations,material->getMaterialRotation());

        OrientationAtIP->rotate_Voigt_stress_from_local_to_global(isv[ip].stress,stress);
        OrientationAtIP->rotate_Voigt_engineering_strain_from_local_to_global(isv[ip].strain,strain);
    
		//PlasticMaterial* material = (PlasticMaterial*)BasicElement::material;
	 //   material->rotateStress(stress,isv[ip].stress,getMaterialAngle(&id),
		//						 materialRotationAxis, false);
	 //   material->rotateEngineeringStrain(strain,isv[ip].strain,getMaterialAngle(&id),
		//						 materialRotationAxis, false);

		for(is=0; is<maxStresses; is++){
			volAvgStress[is] += integrationFactor * stress[is];
			volAvgStrain[is] += integrationFactor * strain[is];
		}
	}  // end of integration loop

	return;
}
//==================================================================



// xtang 052801
void PlasticElement::setIntegrationFactorInISV()
{
	GaussPointList				gaussPointList;
//	InterpolationData			interpData;
	InterpAndSolDerivatives       id(numNodesPerElement,numNodesPerElement); id.allocate();
	//int		i, j;
	
	id.numberOfInterp = numNodesPerElement;
	extractNodalCoordinates(id.xCoor,id.yCoor,id.zCoor);
 
	getQuadraturePoints(gaussPointList); 
	int totalNumIPs = gaussPointList.totalNumIPs;
	double weight;

	for(int ip=0;ip<totalNumIPs;ip++)
	{
		interpolation(ip, gaussPointList, &id);
		jacobian(&id);
		weight = gaussPointList.weightList[ip];
		isv[ip].integrationFactor = id.detJ * weight * getIntegrationFactor();
        //isv[ip].materialOrientation.interpolate_rotation(&id,ElementNodeRotations);
        isv[ip].materialOrientation.interpolate_rotation_by_angle(&id,ElementNodeRotations);
        //isv[ip].materialOrientation.interpolate_old_beta(&id,ElementNodeRotations,material->getMaterialRotation());
		//isv[ip].materialAngle = getMaterialAngle(&id);
	} 
}


//====================================================================
void PlasticElement::getVolAvgISV(double *volAvgStress, 
								  double *volAvgStrain, 
								  double *volAvgPlasticStrain,
								  double &effectiveStress,
								  double &effectivePlasticStrain,
								  double &volume)
//====================================================================

{
// Calculate volume averaged values (strain and stress) and volume

#define maxStresses 6
	double strain[maxStresses];
	double stress[maxStresses];
	double plasticStrain[maxStresses];

	// data for gaussian loop
	GaussPointList				gaussPointList;
	getQuadraturePoints(gaussPointList); 
	int totalNumIPs = gaussPointList.totalNumIPs;

	// initialization
	int is;
	for(is=0; is<maxStresses; is++)
	{   volAvgStress[is] = volAvgStrain[is]=volAvgPlasticStrain[is]=0.; }

	volume=effectiveStress=effectivePlasticStrain=0.0;

	// gauss loop
	for(int ip=0;ip<totalNumIPs;ip++)
	{
		volume += isv[ip].integrationFactor;
		effectiveStress += isv[ip].integrationFactor*isv[ip].effectiveStress;
		effectivePlasticStrain += isv[ip].integrationFactor*isv[ip].effectivePlasticStrain;
	
		PlasticMaterial* material = (PlasticMaterial*)BasicElement::material;

        isv[ip].materialOrientation.rotate_Voigt_stress_from_local_to_global(isv[ip].stress,stress);
        isv[ip].materialOrientation.rotate_Voigt_engineering_strain_from_local_to_global(isv[ip].strain,strain);
        isv[ip].materialOrientation.rotate_Voigt_engineering_strain_from_local_to_global(isv[ip].plasticStrain,plasticStrain);

	 //   material->rotateStress(stress,isv[ip].stress,isv[ip].materialAngle, materialRotationAxis, false);
	 //   material->rotateEngineeringStrain(strain,isv[ip].strain,isv[ip].materialAngle,materialRotationAxis, false);
		//material->rotateEngineeringStrain(plasticStrain,isv[ip].plasticStrain,isv[ip].materialAngle, materialRotationAxis, false);
		for(is=0; is<maxStresses; is++){
			volAvgStress[is] += isv[ip].integrationFactor * stress[is];
			volAvgStrain[is] += isv[ip].integrationFactor * strain[is];
			volAvgPlasticStrain[is] += isv[ip].integrationFactor * plasticStrain[is];  // xtang 052801: this is plastic strain !
		}
	}  // end of integration loop

	return;
}
//==================================================================





////////////////////////////////////////////////////////////////////
///
///   xtang's version of update
///   01-10-01
///
/// 
////////////////////////////////////////////////////////////////////

// we need this struct for passing data among different 

//====================================================================
void PlasticElement::updateXTversion(char * action, double *uGlobal)
//====================================================================
{ 
}
//==========================================================
void PlasticElement::openOutputFiles(FileManager *fm, int flags)
{
	/*
	if(flags == 0){
		BETA_OUT << "No Optional Output Parameters found!" << endl;
		return;
	}
	*/
	eVolDist=fm->OpenBinaryOutputStream("volDist");
}
//=========================================================================
void PlasticElement::closeOutputFiles(FileManager *fm, int flags)
{
	fm->CloseOutputStream((ofstream*)eVolDist);
}
//=========================================================================

//===============================================
//==========static data==========================
//===============================================
std::ostream * PlasticElement::eVolDist = NULL ;  // null-pointer 
