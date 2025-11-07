#include "stdafx.h"

#include "math/gaussQuadrature/hexQuadPoint.hpp"
#include "elements/InterpAndSolDerivatives.hpp"
    
#include "dmElasticElement.hpp"
#include "../material/dmElastic.hpp"
#include "../material/BasicDegradationModel.hpp"
//====================================================================
DMElasticElement::DMElasticElement()
{
	hasNewDamage=false;
}
//===========================================================================
DMElasticElement::~DMElasticElement()
{
}
//===========================================================================
void DMElasticElement::processActionCode(char * action, IntegrationDataType &idt)
{
	int calculateKFlag=0, 
		calculateForceFlag=0, 
		calculateInitialForceFlag=0, 
		calculateStressFlag=0, 
		calculateISVFlag=0, 
		outQuadInfoFlag=0, 
		updateFailureIndexFlag=0;

	if(strstr(action,"K") !=NULL) calculateKFlag             = 1;
    if(strstr(action,"F") !=NULL) {calculateForceFlag        = 1;calculateStressFlag=1;}
	if(strstr(action,"I") !=NULL) calculateInitialForceFlag  = 1;
	if(strstr(action,"S") !=NULL) calculateStressFlag        = 1;
	if(strstr(action,"V") !=NULL) {calculateISVFlag          = 1;calculateStressFlag=1;}
	if(strstr(action,"Q") !=NULL) {outQuadInfoFlag           = 1;calculateStressFlag=1;}
	if(strstr(action,"R") !=NULL) {updateFailureIndexFlag	 = 1;calculateStressFlag=1;}

	if(analysisType==GEOMETRIC_NONLINEAR)
		{calculateStressFlag=1;}

	if(calculateForceFlag==1 && calculateInitialForceFlag==1){
		EXIT_BETA("Cannot calculate current element forces and initial forces during the same call");
    }

	calculateElementVolume(idt);  // point-wise routine
	if (calculateStressFlag==1)			calculateStress(idt); //this needs to be the first option (before calculateFe and calculateKe)
	if (calculateForceFlag==1)			calculateFe(idt);
	if (calculateInitialForceFlag==1)	calculateInitialForce(idt);
	if (calculateKFlag==1)				calculateKe(idt);
	if (calculateISVFlag==1)			calculateISV(idt);
	if (outQuadInfoFlag==1)				outQuadInfo(idt);
	if (updateFailureIndexFlag==1)		updateFailureIndex(idt);
}
//===========================================================================
void DMElasticElement::loopOverIntegrationPoints(char * action, IntegrationDataType &idt)
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

	hasNewDamage=false;

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
		// calculate strain in global coords (not material coords)
        if(idt.uGlobal!=NULL){derivativesOfFields(idt.id);}
        calculateBmatrix(idt.id,b);

        //OrientationAtIP->interpolate_rotation(idt.id,ElementNodeRotations);
        OrientationAtIP->interpolate_rotation_by_angle(idt.id,ElementNodeRotations);
        //OrientationAtIP->interpolate_old_beta(idt.id,ElementNodeRotations,material->getMaterialRotation());
		dmElasticMaterial* material = (dmElasticMaterial*)BasicElement::material;
		double	*strain		= bag->quadPointStrains[idt.ip];

        material->updateLocalizedCmat(OrientationAtIP,isv[idt.ip]);
        //Note, Jamming/Compression Modification no longer works because of the way it interfaces with 
        //the ISV is obsolete. If it's needed in the future, its functionality should be incorporated 
        //into a new DegradationModel
        //The methods will be retained in DMElasticMaterial (but commented out)
        //if (DMElasticElement::compressionModification  &&  !(DMElasticElement::jamming) )//w/ compression modification
			//material->updateLocalizedCmat2(MaterialAngle, materialRotationAxis, isv[idt.ip], strain);
		
		processActionCode(action, idt);
	}
//	if(strstr(action,"K") !=NULL) checkForNegativeDiagonal(); 

}
//===========================================================================
/*
// Jae 981023
void DMElasticElement::getVolAvgValues(double *uGlobal, double *volAvgStress,
	         double *volAvgStrain, double *SED, double *volume)
//====================================================================
{
// Calculate volume averaged values (strain and stress) and volume
// (Jae Noh, 6/23/98)
#define maxStresses 6
//DisplacementDerivatives dispDeriv; 
//InterpolationData       interpData;
	InterpAndSolDerivatives       id(numNodesPerElement,numNodesPerElement); id.allocate();
GaussPointList				gaussPointList;
int     i, j;
double *elemDisp[3];
Matrix *cMat, *sMat, *sMatOrtho, *sMatNonOrtho;
double strain[6], SEDtmp0[6], SEDtmp1[6], SEDtmp2[6], SED2[3];
double *stress;
double *rotStress;

	double **quadPointStresses   = bag->quadPointStresses;
	double **rotatedStresses     = bag->rotatedStresses;

 //------------------------------------------------
 sMatOrtho = new Matrix(6,6);
 sMatNonOrtho = new Matrix(6,6);
 elemDisp[0] = id.nodalValues_u1;
 elemDisp[1] = id.nodalValues_u2;
 elemDisp[2] = id.nodalValues_u3;

 id.numberOfInterp = numNodesPerElement;
 extractNodalCoordinates(id.xCoor,id.yCoor,id.zCoor);
 if( uGlobal != NULL) {extractSolution( uGlobal, elemDisp);} 
 
 initializeArrays();
 getQuadraturePoints(gaussPointList); 
 int totalNumIPs = gaussPointList.totalNumIPs;
double integrationFactor, weight;

    *volume=0;  //xtang 051302


int is;

	for(is=0; is<maxStresses; is++){   volAvgStress[is] = volAvgStrain[is]=0.; }
	SED[0]=SED[1]=SED[2]=0.0;

//...............................................................
//                      Integration loop
//...............................................................
for(int ip=0;ip<totalNumIPs;ip++)
{
	interpolation(ip, gaussPointList, &id);
    jacobian(&id);
	weight = gaussPointList.weightList[ip];
	integrationFactor = id.detJ * weight * getIntegrationFactor();
    *volume += integrationFactor;

	dmElasticMaterial* material = (dmElasticMaterial*)BasicElement::material;
    material->updateLocalizedCmat(getMaterialAngle(&id),
  	                              materialRotationAxis, &isv[ip]);
    cMat = material->getPointerlocalizedCmat();  // cMat->print("cmat");
                
    if(uGlobal!=NULL){derivativesOfFields(&id);}
	calculateBmatrix(&id,b);
	stress      = quadPointStresses[ip];
	rotStress   = rotatedStresses[ip];

	calculateStrain(&id,strain); 
	// xtang 122000
    material->calculateStress(stress,strain,stateVariable,getMaterialAngle(&id),
				materialRotationAxis); 
    material->rotateStress(rotStress,stress,-getMaterialAngle(&id),
									 materialRotationAxis, true);
   // calculate strain energy density
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

   for(is=0; is<maxStresses; is++){
	   volAvgStress[is] += integrationFactor * stress[is];
	   volAvgStrain[is] += integrationFactor * strain[is];
	}
}  // end of integration loop

bag->totalVolume += *volume;  //JV061407 

 delete sMatOrtho;    // xt062901;
 delete sMatNonOrtho; // xt062901;
}
*/
//=========================================================================
void DMElasticElement::getVolAvgValues(double *uGlobal, double *volAvgStress,
	         double *volAvgStrain, double *SED, double *volume)
//====================================================================
{
// Calculate volume averaged values (strain and stress) and volume

//JV111408 : removed code related to calculating SED, put it back when actually needed
	update("F", uGlobal);

	for(int i=0; i<numStrains; i++){
		volAvgStrain[i] = bag->strainVolume_element[i];
		volAvgStress[i] = bag->stressVolume_element[i];
	}
	*volume=bag->elementVolume;

	bag->totalVolume += *volume;  //JV061407 
}
//=========================================================================
void DMElasticElement::printDamageState(ostream *ost)
{
    int num_of_ip=getTotalNumberOfIPs();
    int mat=getMaterialNumber();

    dmElasticMaterial* dmMaterial = (dmElasticMaterial*)material;

    (*ost) << elementNumber << "\t" <<mat<<endl;
    dmMaterial->printDamageState(isv, num_of_ip, ost);
}
//=========================================================================
void DMElasticElement::printDamageStateBinary(ostream *ost)
{
        int numIPs=getTotalNumberOfIPs();
		
		int mat=getMaterialNumber();
		ost->write((char*) &elementNumber, sizeof(int));
		ost->write((char*) &mat, sizeof(int));

		
		//New version: +DG_oct_13_2006
		dmElasticMaterial* dmMaterial = (dmElasticMaterial*)material;
		for (int j=0; j<numIPs; j++) 
		{
            DamageISV* isvcp=(DamageISV*) isv[j];  //Creating a copy of isv
			unsigned FailureStatus=0;
			char info;
			info=dmMaterial->DegradationModel->getDamageFlag(isvcp->TriggeredFailureModes[0], isvcp->damFactorForE11); ost->write((char*) &info, sizeof(char));
			info=dmMaterial->DegradationModel->getDamageFlag(isvcp->TriggeredFailureModes[1], isvcp->damFactorForE22); ost->write((char*) &info, sizeof(char));
            info=dmMaterial->DegradationModel->getDamageFlag(isvcp->TriggeredFailureModes[2], isvcp->damFactorForE33); ost->write((char*) &info, sizeof(char));
			info=dmMaterial->DegradationModel->getDamageFlag(isvcp->TriggeredFailureModes[3], isvcp->damFactorForG12); ost->write((char*) &info, sizeof(char));
			info=dmMaterial->DegradationModel->getDamageFlag(isvcp->TriggeredFailureModes[4], isvcp->damFactorForG23); ost->write((char*) &info, sizeof(char));
			info=dmMaterial->DegradationModel->getDamageFlag(isvcp->TriggeredFailureModes[5], isvcp->damFactorForG13); ost->write((char*) &info, sizeof(char));
			FailureStatus = isvcp->TriggeredFailureModes[0] + 
				isvcp->TriggeredFailureModes[1] + 
				isvcp->TriggeredFailureModes[2] + 
				isvcp->TriggeredFailureModes[3] + 
				isvcp->TriggeredFailureModes[4] + 
				isvcp->TriggeredFailureModes[5] ;
			if (FailureStatus != 0) {info='1';} else {info='0';}  
			ost->write((char*) &info, sizeof(char));
        }
		//New version: -DG_oct_13_2006
}
//=========================================================================
void DMElasticElement::LoadDamageState(istream *is)
{
	int num_of_ip=getTotalNumberOfIPs();
	dmElasticMaterial* dmMaterial = (dmElasticMaterial*)material;

	int ele, matID;
	*is>>ele>>matID;

	dmMaterial->LoadDamageState(isv, num_of_ip, is);
}
//=========================================================================
void DMElasticElement::calculateInitialForce(IntegrationDataType &idt)
{
	double *stress	= bag->quadPointStresses[idt.ip];

	dmElasticMaterial* material=(dmElasticMaterial*)BasicElement::material;
	material->calculateInitialStress(stress, isv[idt.ip], OrientationAtIP); 
	addToFe(numDof,numStrains,b,stress,idt.integrationFactor,initialFe);
}
//===========================================================================
void DMElasticElement::calculateKe(IntegrationDataType &idt)
{
	// local properties that require transformations	   
	dmElasticMaterial* material=(dmElasticMaterial*)BasicElement::material;
	Matrix *cMat= material->getPointerToGlobalCmat(); 
	
	double *stress	= bag->quadPointStresses[idt.ip];
	addToKe(cMat, idt.id, idt.integrationFactor,stress);
}
//===========================================================================
void DMElasticElement::calculateISV(IntegrationDataType &idt)
{
	double	*rotStress   = bag->rotatedStresses[idt.ip];
	dmElasticMaterial* mat=(dmElasticMaterial*) material;

	hasNewDamage =( mat->UpdateISV(isv[idt.ip], rotStress) | hasNewDamage );//the function is the first argument in the 'or' parameter list because 
	//I have heard that in some cases if the first parameter evaluates to true, the program will not evaluate the other parameters.
	//it appears that atleast in VC8 (debug), that is not the case. dunno about release version or supercomputer version. 
}
//===========================================================================
void DMElasticElement::outQuadInfo(IntegrationDataType &idt)
{
	double x,y,z;
	x=idt.id->x;
	y=idt.id->y;
	z=idt.id->z;
	double weight = idt.gaussPointList.weightList[idt.ip];
	double vol = idt.integrationFactor;

	double	*rotStress;
	double **rotatedStresses     = bag->rotatedStresses;
	rotStress   = rotatedStresses[idt.ip];

	int group=getMaterialNumber();		

	//(*eQuadInfo).write((char*)(isv[ip].ratio),sizeof(isv[ip].ratio));  // this is a binary output

	(*eQuadInfo).write((char*)(&rotStress[0]),sizeof(double));  // this is a binary output
	(*eQuadInfo).write((char*)(&rotStress[1]),sizeof(double));  // this is a binary output
	(*eQuadInfo).write((char*)(&rotStress[2]),sizeof(double));  // this is a binary output
	(*eQuadInfo).write((char*)(&rotStress[3]),sizeof(double));  // this is a binary output
	(*eQuadInfo).write((char*)(&rotStress[4]),sizeof(double));  // this is a binary output
	(*eQuadInfo).write((char*)(&rotStress[5]),sizeof(double));  // this is a binary output
	(*eQuadInfo).write((char*)(&rotStress[0]),sizeof(double));  // this is a binary output //REINPUT


	(*eQuadInfo).write((char*) &group,sizeof(group));  // this is a binary output
	(*eQuadInfo).write((char*) &vol,sizeof(vol));  // this is a binary output
	(*eQuadInfo).write((char*) &x,sizeof(x));  // this is a binary output
	(*eQuadInfo).write((char*) &y,sizeof(y));  // this is a binary output
	(*eQuadInfo).write((char*) &z,sizeof(z));  // this is a binary output
}
//===========================================================================
void DMElasticElement::updateFailureIndex(IntegrationDataType &idt)
{
    double	*rotStress   = bag->rotatedStresses[idt.ip];
	dmElasticMaterial* mat=(dmElasticMaterial*) material;
    ((DamageISV*)isv[idt.ip])->SwapRatioAndLastRatio(); //put current Failure Index into Last Failure Index (swaps pointers)
	mat->updateFailureIndex(isv[idt.ip], rotStress);
}
//===========================================================================
void DMElasticElement::calculateStress(IntegrationDataType &idt)
{
	// local properties that require transformations	   
	dmElasticMaterial* material = (dmElasticMaterial*)BasicElement::material;

	double	*strain		= bag->quadPointStrains[idt.ip];
	double	*stress		= bag->quadPointStresses[idt.ip];
	double	*rotStress	= bag->rotatedStresses[idt.ip];
	double	*rotStrain	= bag->rotatedStrains[idt.ip];
    
	calculateStrain(idt.id, strain); //in global coordinate system
   	material->calculateStress(stress, 
                              strain, 
                              isv[idt.ip],
                              OrientationAtIP); // include Hygrothermal stress: in global coord
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

void DMElasticElement::swapRatioAndLastRatioInISV()
{
    for(int j=0;j<getTotalNumberOfIPs();++j){
        DamageISV* Disv = (DamageISV*) isv[j];
        Disv->SwapRatioAndLastRatio();
	}
}

//===========================================================================

void DMElasticElement::getDeltaToNextFailure(const double &trialDelta, double& currentMinDeltaLoadFactor, const double &failureRatio) const
{
    for(int j=0;j<getTotalNumberOfIPs();++j){
        DamageISV* Disv = (DamageISV*) isv[j];
        Disv->getDeltaToNextFailure(trialDelta,currentMinDeltaLoadFactor,failureRatio);
	}//End Loop over Quad Points
}

//===========================================================================

//===============================================
//==========static data==========================
//===============================================
bool	DMElasticElement::compressionModification=false;
bool	DMElasticElement::jamming=false;