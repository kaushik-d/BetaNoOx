#include "stdafx.h"

#include <fstream>
#include <iomanip>
#include <string.h>

#include "TransientElementWorkspace.hpp"
#include "HeatTransferElement3DTransient.hpp"
#include "math/gaussQuadrature/hexQuadPoint.hpp"

//====================================================================
extern	int			verboseFlag;
//====================================================================

//====================================================================
HeatTransferElement3Dtransient::
                       HeatTransferElement3Dtransient(void)
//====================================================================
{ 
	//if(traceFlag==1) {BETA_OUT<< "Creating HeatTransferElement3Dtransient element"<<endl;}
	numStrains       = 3;    // Eventually fix... these are the three temp. grad
	setIntegrationOrder(3);
	numberOfIntegrationDirections=3;

	F1=0;
	F2=0;
	M_deltaq=0;
	Me=0;
	MeAlreadyCalculated=false;
}
//====================================================================

HeatTransferElement3Dtransient::
                       ~HeatTransferElement3Dtransient(void)
{ 
	//BETA_OUT<< "Deleting HeatTransferElement3Dtransient element"<<endl; 
}

//====================================================================
 
 void HeatTransferElement3Dtransient::
							 addToFe(int numDof, int numStrains,
							 InterpAndSolDerivatives *id_ptr,
							 double *flux,
							 double integrationFactor,
							 double *force)
{
	InterpAndSolDerivatives &id=*id_ptr;
	double  *p_S_x1=id.p_S_x1,
		  *p_S_x2=id.p_S_x2,
		  *p_S_x3=id.p_S_x3;		
	  for(int i = 0; i<numDof; i++)
		 {
        force[i] += ( flux[0] * p_S_x1[i] + 
			          flux[1] * p_S_x2[i] + 
					  flux[2] * p_S_x3[i] ) * integrationFactor;
		 }
  }
//==================================================================
void HeatTransferElement3Dtransient::addToMe(Matrix* cMat,InterpAndSolDerivatives *id_ptr,
									  double integrationFactor)       
{	
	InterpAndSolDerivatives &id=*id_ptr;
double *S=id. S;
HeatTransferMaterial* material=(HeatTransferMaterial*)BasicElement::material;
double saturationValue = material->getSaturationValue();

int i,j;
double Si=0.0;
// mass matrix: xtang 121200
	for(i=0; i<numDof; i++) {
		Si = S[i];
		for(j=i; j<numDof; j++)	{
			(*Me)[i][j] += Si*S[j]*integrationFactor*saturationValue;
		}
	 }
}  //end of addToMe
//======================================================================
void HeatTransferElement3Dtransient::
									gradientsOfFieldVariables(InterpAndSolDerivatives *id_ptr)
{
	InterpAndSolDerivatives &id=*id_ptr;
	int i;
	int numberOfInterp = id.numberOfInterp;
	double *S=id.S,
		 *p_S_Lx1=id.p_S_Lx1,   *p_S_x1=id.p_S_x1, 
		 *p_S_Lx2=id.p_S_Lx2,   *p_S_x2=id.p_S_x2, 
		 *p_S_Lx3=id.p_S_Lx3,   *p_S_x3=id.p_S_x3;
   
	double *xDisp=id.nodalValues_u1;

	double p_u1_x1, p_u1_x2, p_u1_x3;

	p_u1_x1 = p_u1_x2 = p_u1_x3 = 0.;

	
	for(i=0; i< numberOfInterp; i++)
	{
				p_u1_x1 += p_S_x1[i] * xDisp[i];
				p_u1_x2 += p_S_x2[i] * xDisp[i];
				p_u1_x3 += p_S_x3[i] * xDisp[i];				

	}
	id.p_u1_x1 = p_u1_x1;
	id.p_u1_x2 = p_u1_x2;
	id.p_u1_x3 = p_u1_x3;
	
}  // end 


// xtang 011001: only one value -- concentration
double HeatTransferElement3Dtransient::
          pointValueOfFieldVariables(InterpAndSolDerivatives *id_ptr, double *nodeValue)
{
	InterpAndSolDerivatives &id=*id_ptr;
	int i;
	int numberOfInterp = id.numberOfInterp;
	double *S=id.S; 
	double u1=0.0;
	
	for(i=0; i< numberOfInterp; i++){u1 +=S[i] * nodeValue[i];}

	return (u1);
} 


//====================================================================
void HeatTransferElement3Dtransient::printNodalFluxesInLCS(ostream &outstream)
{
    outstream << elementNumber << " " << getMaterialNumber() << endl;
    extrapolateQuadValuesToNodes(bag->rotatedFluxes,numStrains,bag->extrapolatedFluxes);
    for(int i=0;i<numNodesPerElement;++i){
        for(int j=0;j<numStrains;++j){
            outstream << bag->extrapolatedFluxes[i][j] << " ";
        }
        outstream << node[i].nodeNum << endl;
    }
}

//====================================================================
void HeatTransferElement3Dtransient::printNodalFluxesInGCS(ostream &outstream)
{
    outstream << elementNumber << " " << getMaterialNumber() << endl;
    extrapolateQuadValuesToNodes(bag->quadPointFluxes,numStrains,bag->extrapolatedFluxes);
    for(int i=0;i<numNodesPerElement;++i){
        for(int j=0;j<numStrains;++j){
            outstream << bag->extrapolatedFluxes[i][j] << " ";
        }
        outstream << node[i].nodeNum << endl;
    }
}
//====================================================================================
void HeatTransferElement3Dtransient::addToKe(Matrix* cMat,InterpAndSolDerivatives *id_ptr,
									  double integrationFactor)       
{	
	//saturation value already factored into the cMat
	InterpAndSolDerivatives &id=*id_ptr;
	double *S=id.S,
       *p_S_x1=id.p_S_x1,    
	   *p_S_x2=id.p_S_x2,     
	   *p_S_x3=id.p_S_x3;	

//	int    ii,jj;
double k11,k22,k33;
k11 = (*cMat)(0,0);
k22 = (*cMat)(1,1);
k33 = (*cMat)(2,2);

double k12,k13,k23;
k12 = (*cMat)(0,1);
k13 = (*cMat)(0,2);
k23 = (*cMat)(1,2);

//BETA_OUT<<" k11,k22,k33 "<<k11<<"  "<<k22<<"  "<<k33<<endl;

// Limit to orthotropic for now
// double k12=0,k13=0,k23=0;

int i,j;
double p_S_x1i=0.0;
double p_S_x2i=0.0;
double p_S_x3i=0.0;

	for(i=0; i<numDof; i++) {
		p_S_x1i = p_S_x1[i];
		p_S_x2i = p_S_x2[i];
		p_S_x3i = p_S_x3[i];
		for(j=i; j<numDof; j++)	{
			(*Ke)[i][j] += (
				   (p_S_x1[j]*k11+p_S_x2[j]*k12+p_S_x3[j]*k13)*p_S_x1i+
			       (p_S_x1[j]*k12+p_S_x2[j]*k22+p_S_x3[j]*k23)*p_S_x2i+
			  	   (p_S_x1[j]*k13+p_S_x2[j]*k23+p_S_x3[j]*k33)*p_S_x3i )
                   * integrationFactor;
		}
	 }

}  //end of addToKe
//====================================================================================
void HeatTransferElement3Dtransient::readSpecialCommand(istream &inStream, ElementGroup *element,
								    int numElements, char * command)
{int i,count;

HeatTransferElement3Dtransient *p_elem;


BETA_OUT<<"HeatTransferElement3Dtransient::readSpecialCommand"<<endl;

//================================================================
if (COMPARE(command, "SETEFFECTIVEPROPERTYCASE")==0 ) {
	inStream>> effectivePropertyCase;
	inStream>> configDescription;
	BETA_OUT<<"EFFECTIVE_MODULUS_CASE = "<< effectivePropertyCase<<endl;
	BETA_OUT<<"Configuration Description = "<< configDescription<<endl;
	return;
}; //setEffectiveModulusCase



// If no match, try super class
BETA_OUT<<"No match... try super class"<<endl;
IsoElement::readSpecialCommand(inStream,element,numElements, command);
}//readSpecialCommand
//==========================================================
/*
//#include "vocabulary.hpp"
void HeatTransferElement3Dtransient::readSpecialCommand(int command, istream * inStream)
{

int i,count;

BETA_OUT<<"Reading special commands for element type HeatTransferElement3Dtransient"<<endl;

switch(command)
{
case SET_MATERIAL_ANGLES:

  //One rotation direction, but rotation can be different at each node
  (*inStream)>> materialRotationAxis;
  // ZOUT<<"Rotation axis = "<< materialRotationAxis<<endl;
  
   materialRotationAnglesAtNodes = new double[numNodesPerElement];
   count=0;
   for(i=0; i<numNodesPerElement; i++)
      { (*inStream)>> materialRotationAnglesAtNodes[i];
      //  ZOUT<<materialRotationAnglesAtNodes[i]<<" ";
       count++;  // if(count==10){ZOUT<<endl; count=0;}
      } // ZOUT<<endl;
   break;


case SET_EFFECTIVE_PROPERTY_CASE:
  (*inStream)>> effectivePropertyCase;
  (*inStream)>> configDescription;
  BETA_OUT<<"EFFECTIVE_MODULUS_CASE = "<< effectivePropertyCase<<endl;
  BETA_OUT<<"Configuration Description = "<< configDescription<<endl;
  break;
 
default:
   BETA_OUT<<"No match..."<<endl; exit(1);
 } // end of switch


}
*/
//======================================================================
void HeatTransferElement3Dtransient::initializeArrays()
{
	for(int i=0;i<numDof;i++) {
		Fe[i] = 0.;
		initialFe[i] = 0.;
		for(int j=0;j<numDof;j++) {
			(*Ke)[i][j] = 0.;
		    (*Me)[i][j] = 0.;
		}
	}
}
//=====================================================================
void HeatTransferElement3Dtransient::setPointers(ElementWorkspace * s)
{	
	BasicElement::setPointers(s);
	TransientElementWorkspace* workspace=(TransientElementWorkspace*)s;
	F1				= workspace->F1;
	F2				= workspace->F2;
	M_deltaq		= workspace->M_deltaq;

	quadPointFluxes  = s->quadPointStresses;
	
	if(!requiresParallelESAConsideration){
		Ke_orig          = s->Ke_original; // xtang 07082003

		if(workspace->StoreAllMe)
			Me			= &((*s->MeList)[elementNumber]);
		else
			Me			= s->Me;  //xtang 121200

	}
	
	F_orig           = s->F_original; // xtang 07082003
	bag              = s;
};
//=====================================================================
void HeatTransferElement3Dtransient::initializeSummary()
{
	int i;

	// xtang 011001 : why initializing here ?
	bag->totalVolume = 0.;  
    bag->temperatureVolume = 0.;  //xtang022001
	for(i=0; i<numStrains; i++) {
		bag->fluxVolume_total[i] = 0.; 
		bag->gradientVolume_total[i] = 0.;
	}
}
//=====================================================================
void HeatTransferElement3Dtransient::addToSummary()
{
	
}

//**************************************************
// xtang:  temporary for transient

void HeatTransferElement3Dtransient::formK_s(Array<double> &AnalysisParameters, double timeStep)
{
	double a1=AnalysisParameters[0];
	double a2=AnalysisParameters[1];

	int i,j;

	if(verboseFlag == Max){
		BETA_OUT <<"Element #" << elementNumber << endl;
		Ke->print("Ktangent_e", numDof,numDof,&BETA_OUT);
		Me->print("M_e", numDof,numDof,&BETA_OUT);
	}

	//for (i=0;i<numDof; i++) for (j=0;j<numDof;j++) (*Ke_orig)[i][j]=0.0; // xtang 07082003 //commented out BCO 010810 already zeroed in loopOverIntegPoints()

	//maybe optimize later... for cases where a1==0
	for(i=0; i<numDof; i++) {
		for(j=0; j<numDof; j++)	{
			(*Ke_orig)[i][j] = (a1*(*Ke)[i][j]+(*Me)[i][j])/timeStep; // 07082003
//			(*Ke)[i][j] = (a1*(*Ke_orig)[i][j]+(*Me)[i][j]);
		}
	 }
	if(verboseFlag == Max){
//		Ke->print("Mbar_e", numDof,numDof,&BETA_OUT);
		int *doflist=getDofList();
		Ke_orig->printWithRowColumnNumbering("Mbar_e", getNumDof(), getNumDof(), doflist, doflist, &BETA_OUT);
	}
}
//====================================================================
void HeatTransferElement3Dtransient::formF_s(Array<double> &AnalysisParameters, double *u, double *ut, double timeStep)
{
	double a1=AnalysisParameters[0];
	double a2=AnalysisParameters[1];

	int i,j;
	double dump;

	for (i=0;i<numDof; i++) Fe[i]=0.0; // xtang 07082003

	for(i=0; i<numDof; i++) {
		dump=0.0;
		for(j=0; j<numDof; j++)	{
			dump+=( (*Me)[i][j]- a2*(*Ke_orig)[i][j] )*u[dofList[j]];  //07082003
		}
		Fe[i] += dump/timeStep; 
//		Fe[i] += dump;
	 }
}
//====================================================================
void HeatTransferElement3Dtransient::formF_s_new(Array<double> &AnalysisParameters, double *u_s, double *delta_u, double *u_s1, double timeStep)
{
	double a1=AnalysisParameters[0];
	double a2=AnalysisParameters[1];

	int i,j;
	double dump;

	for (i=0;i<numDof; i++) Fe[i]=0.0; // xtang 07082003

	for(i=0; i<numDof; i++) {
		dump=0.0;
		for(j=0; j<numDof; j++)	{
			double Ku = 0.0;
			if(j>=i){
			Ku = (*Ke)[i][j] * u_s[dofList[j]];
			}
			else{
			Ku = (*Ke)[j][i] * u_s[dofList[j]];
			}
//			dump+=( (*Me)[i][j]- a2*(*Ke_orig)[i][j] )*u[dofList[j]];  //07082003
			dump+= - (a1+a2) * Ku;  //new K and old u
		}
		Fe[i] += dump/timeStep; 
//		Fe[i] += dump;
	 }

	if(verboseFlag == Max){
		BETA_OUT <<"Element #" << elementNumber << endl;
		BETA_OUT << "FBar_e" << endl;
		for (i=0;i<numDof; i++) BETA_OUT << "| " << Fe[i] << endl;
	}

}
//====================================================================
void HeatTransferElement3Dtransient::formEnhancedResultantVector(Array<double> &AnalysisParameters, double *u_s, double *delta_u, double *u_s1, double timeStep)
{
	double a1=AnalysisParameters[0];
	double a2=AnalysisParameters[1];

	int i;
	double dump;

	for (i=0;i<numDof; i++) F_orig[i] = -F_orig[i]; // this is very important
	// this is because when the force vector is calculated, it already assumes the expression is on the right hand side of the equation.
	// it calculates Nx[i] X flux... whereas in the documentation, we are supposed to calculate Nx[i] x D x gradient (which is negative of the first expression)

//	BETA_OUT << "\nF(a1) " <<endl; 	for (i=0;i<numDof; i++) BETA_OUT <<"| " << F_orig[i] << endl;
	for (i=0;i<numDof; i++) Fe[i]=0.0;

	for(i=0; i<numDof; i++) {
		dump=  M_deltaq[i] + a1*F_orig[i];
		Fe[i] += dump/timeStep; 
//		Fe[i] += dump;
	 }
}
//====================================================================
void HeatTransferElement3Dtransient::formResultantVector(Array<double> &AnalysisParameters, double *u_s, double *delta_u, double *u_s1, double timeStep)
{
	double a1=AnalysisParameters[0];
	double a2=AnalysisParameters[1];

	int i,j;
	double dump;

	for (i=0;i<numDof; i++) F_orig[i] = -F_orig[i]; // this is very important
	// this is because when the force vector is calculated, it already assumes the expression is on the right hand side of the equation.
	// it calculates Nx[i] X flux... whereas in the documentation, we are supposed to calculate Nx[i] x D x gradient (which is negative of the first expression)

//	BETA_OUT << "\nF(a1) " <<endl; 	for (i=0;i<numDof; i++) BETA_OUT <<"| " << F_orig[i] << endl;
	for (i=0;i<numDof; i++) Fe[i]=0.0;

	for(i=0; i<numDof; i++) {
		dump=0.0;
		for(j=0; j<numDof; j++)	{
			dump+=  (*Me)[i][j]*delta_u[dofList[j]];
		}
		dump+=  a1*F_orig[i];
		Fe[i] += dump/timeStep; 
//		Fe[i] += dump;
	 }
}
//====================================================================
void HeatTransferElement3Dtransient::addA2TermToResultantVector(Array<double> &AnalysisParameters, double *u_s, double *delta_u, double *u_s1, double timeStep)
{
	double a1=AnalysisParameters[0];
	double a2=AnalysisParameters[1];

	int i;
	double dump;

	for (i=0;i<numDof; i++) F_orig[i] = -F_orig[i]; // this is very important
	// this is because when the force vector is calculated, it already assumes the expression is on the right hand side of the equation.
	// it calculates Nx[i] X flux... whereas in the documentation, we are supposed to calculate Nx[i] x D x gradient (which is negative of the first expression)


//	BETA_OUT << "\nF(a2) " <<endl; 	for (i=0;i<numDof; i++) BETA_OUT <<"| " << F_orig[i] << endl;

	for(i=0; i<numDof; i++) {
		dump=  a2*F_orig[i];
		Fe[i] += dump/timeStep; 
//		Fe[i] += dump;
	 }
}
//====================================================================
void HeatTransferElement3Dtransient::formFlowVector(Array<double> &AnalysisParameters, double *u_s, double *delta_u, double *u_s1, double timeStep)
{
	//this function will not divide the vector by the timestep (as is done in all the other functions)
	//might have to avoid dividing by timestep in all the other functions (it was something XT did)
	double a1=AnalysisParameters[0];
	double a2=AnalysisParameters[1];

	int i,j;
	double dump;

	for (i=0;i<numDof; i++) F_orig[i] = -F_orig[i]; // this is very important
	// this is because when the force vector is calculated, it already assumes the expression is on the right hand side of the equation.
	// it calculates Nx[i] X flux... whereas in the documentation, we are supposed to calculate Nx[i] x D x gradient (which is negative of the first expression)

//	BETA_OUT << "\nF(a1) " <<endl; 	for (i=0;i<numDof; i++) BETA_OUT <<"| " << F_orig[i] << endl;
	for (i=0;i<numDof; i++) Fe[i]=0.0;

	for(i=0; i<numDof; i++) {
		dump=0.0;
		for(j=0; j<numDof; j++)	{
			dump+=  (*Me)[i][j]*delta_u[dofList[j]] /a1;
		}
		dump+=  F_orig[i];
		Fe[i] += dump;
	}
}
//====================================================================
void HeatTransferElement3Dtransient::addA2TermToFlowVector(Array<double> &AnalysisParameters, double *u_s, double *delta_u, double *u_s1, double timeStep)
{
	//this function will not divide the vector by the timestep (as is done in all the other functions)
	//might have to avoid dividing by timestep in all the other functions (it was something XT did)
	double a1=AnalysisParameters[0];
	double a2=AnalysisParameters[1];

	int i;
	double dump;

	for (i=0;i<numDof; i++) F_orig[i] = -F_orig[i]; // this is very important
	// this is because when the force vector is calculated, it already assumes the expression is on the right hand side of the equation.
	// it calculates Nx[i] X flux... whereas in the documentation, we are supposed to calculate Nx[i] x D x gradient (which is negative of the first expression)

//	BETA_OUT << "\nF(a2) " <<endl; 	for (i=0;i<numDof; i++) BETA_OUT <<"| " << F_orig[i] << endl;

	for(i=0; i<numDof; i++) {
		dump=  a2*F_orig[i];
		Fe[i] += dump;
	 }
}
//====================================================================
void HeatTransferElement3Dtransient::outputSummary()
{
#define FMT1 setw(5)
#define FMT2 setw(13)<<setprecision(4)
double totalVolume = bag->totalVolume;

/*
BETA_OUT<<"Total volume of model="<<totalVolume <<endl;

BETA_OUT<<"The current calculation assumes there are no cracks or holes in the body\n";
BETA_OUT<<"... if there are, the volume averages are incorrect!\n";
 
		  BETA_OUT<<"Volume averaged fluxes"<<endl;
          for(j=0; j<3; j++) 
          {BETA_OUT<<FMT2<<bag->fluxVolume_total[j]/totalVolume<<"  ";}
          BETA_OUT<<endl;    
		  
		  BETA_OUT<<"Volume averaged gradients p_T_x1, p_T_x2, p_T_x3"<<endl;
          for(j=0; j<3; j++) 
          {BETA_OUT<<FMT2<<bag->gradientVolume_total[j]/totalVolume<<"  ";}
          BETA_OUT<<endl;  
		  
		  BETA_OUT<<"<flux> / <gradient"<<endl;
		  for(j=0; j<3; j++) 
          {BETA_OUT<<FMT2<<bag->fluxVolume_total[j]/bag->gradientVolume_total[j]<<"  ";}
          BETA_OUT<<endl;   
		  

BETA_OUT<<"effectivePropertyCase="<<effectivePropertyCase<<endl;
*/
if(effectivePropertyCase > -1){
	//FIX!!
		//REPORT<<"configDescription="<< configDescription<<endl;
	    //REPORT<<"effectivePropertyCase  "<<effectivePropertyCase <<endl;
	}

/*
if(effectivePropertyCase==1){
	REPORT<<"          kxx = "<< -bag->fluxVolume_total[0]/bag->gradientVolume_total[0]<<endl;
} // 
if(effectivePropertyCase==2){
	REPORT<<"          kyy = "<< -bag->fluxVolume_total[1]/bag->gradientVolume_total[1]<<endl;
} // 
if(effectivePropertyCase==3){
	REPORT<<"          kzz = "<< -bag->fluxVolume_total[2]/bag->gradientVolume_total[2]<<endl;
} // 
*/
}//end of method
//===========================================================================
void HeatTransferElement3Dtransient::calculateFe(IntegrationDataType &idt)
{
   addToFe(numDof,numStrains,
	   idt.id,quadPointFluxes[idt.ip],
	   idt.integrationFactor,
	   F_orig);
}
//===========================================================================
void HeatTransferElement3Dtransient::calculateMe(IntegrationDataType &idt)
{
	addToMe(0, idt.id, idt.integrationFactor);
}
//===========================================================================
void HeatTransferElement3Dtransient::calculateKe(IntegrationDataType &idt)
{
	// local properties that require transformations	   
	Matrix *cMat;
	HeatTransferMaterial* material=(HeatTransferMaterial*)BasicElement::material;
    //OrientationAtIP->interpolate_rotation(idt.id,ElementNodeRotations);
    OrientationAtIP->interpolate_rotation_by_angle(idt.id,ElementNodeRotations);
    //OrientationAtIP->interpolate_old_beta(idt.id,ElementNodeRotations,material->getMaterialRotation());
    material->rotateMaterialToGlobal(OrientationAtIP);

	double normc=pointValueOfFieldVariables(idt.id, idt.id->nodalValues_u1);
	material->UpdateGlobalCmat(normc);

	cMat = material->getPointerToGlobalCmat(); 
	addToKe(cMat, idt.id, idt.integrationFactor);
}
//===========================================================================
void HeatTransferElement3Dtransient::calculateFlux(IntegrationDataType &idt)
{
	//this function also increments the concentration in the element
	double *flux;  //local
	double *gradientVolume_element  = bag->gradientVolume_element;
	double *fluxVolume_element      = bag->fluxVolume_element;
	double &concentrationVolume_element= bag->temperatureVolume_element;
	double normc=pointValueOfFieldVariables(idt.id, idt.id->nodalValues_u1);

    if(idt.uGlobal!=NULL) { gradientsOfFieldVariables(idt.id); } ;
	flux  = quadPointFluxes[idt.ip]; 	

	HeatTransferMaterial* material=(HeatTransferMaterial*)BasicElement::material;
    //OrientationAtIP->interpolate_rotation(idt.id,ElementNodeRotations);
    OrientationAtIP->interpolate_rotation_by_angle(idt.id,ElementNodeRotations);
    //OrientationAtIP->interpolate_old_beta(idt.id,ElementNodeRotations,material->getMaterialRotation());
    material->rotateMaterialToGlobal(OrientationAtIP);
	material->UpdateGlobalCmat(normc);

	material->calculateFlux(flux, idt.id);	
	double saturationValue = material->getSaturationValue();

	for(int i=0; i<numStrains; i++) { fluxVolume_element[i]+= flux[i] *idt.integrationFactor;}

	gradientVolume_element[0] += idt.id->p_u1_x1 *idt.integrationFactor *saturationValue;
	gradientVolume_element[1] += idt.id->p_u1_x2 *idt.integrationFactor *saturationValue;
	gradientVolume_element[2] += idt.id->p_u1_x3 *idt.integrationFactor *saturationValue;

    concentrationVolume_element += normc*idt.integrationFactor*saturationValue;
}
//===========================================================================
//.....................................................
//  This is the place that user need to put new actions
//.....................................................

void HeatTransferElement3Dtransient::processActionCode(char * action, IntegrationDataType &idt)
{
	int calculateKFlag=0, calculateForceFlag=0, calculateFluxFlag=0, calculateMFlag=0;
	if(strstr(action,"K") !=NULL) calculateKFlag = 1;
	if(strstr(action,"M") !=NULL) calculateMFlag = 1;
	if(strstr(action,"F") !=NULL) {calculateForceFlag = 1;
								   calculateFluxFlag  = 1;}
	if(strstr(action,"S") !=NULL) calculateFluxFlag = 1;

	// 1: volume of quadrature point
	calculateElementVolume(idt);  // point-wise routine
	// 2: flux for the current point
	if (calculateFluxFlag==1) calculateFlux(idt);
	// 3: Fe for the current point and then assembled to the element force vector
	if (calculateForceFlag==1) calculateFe(idt);
	// 4: Ke for the current point and then assembled to the element force vector
	if (calculateKFlag==1) calculateKe(idt);
	// 5: Me for the current point 
	if (calculateMFlag==1) calculateMe(idt);
}

//.................................................
//   This routine will not be changed 
//.................................................

void HeatTransferElement3Dtransient::loopOverIntegrationPoints(char * action, IntegrationDataType &idt){
	int totalNumIPs = idt.gaussPointList.totalNumIPs;
	// location
	double weight;
	
	// xtang 07082003
	int i,j;
	if(strstr(action,"K") !=NULL) {
		for (i=0;i<numDof; i++) for (j=0;j<numDof;j++) (*Ke_orig)[i][j] = (*Ke)[i][j]=0.0;
	}
	if(strstr(action,"M") !=NULL) {
		for (i=0;i<numDof; i++) for (j=0;j<numDof;j++) (*Me)[i][j]=0.0;
	}
	if(strstr(action,"F") !=NULL) {
		for (i=0;i<numDof; i++) F_orig[i] =  0.0;
	}
	if(strstr(action,"I") !=NULL) {
		initialFe=0;
		return;
	}



	for(int ip=0;ip<totalNumIPs;ip++) {
		// set up data or the current point: 
		// these data must be obtained for all routines and 
		// are basically for the integration factor
		interpolation(ip, idt.gaussPointList, idt.id);
		jacobian(idt.id);
		weight = idt.gaussPointList.weightList[ip];
		idt.integrationFactor = idt.id->detJ * weight * getIntegrationFactor();
		idt.ip=ip;
		processActionCode(action, idt);
	}

	// xtang 07082003
	//if(strstr(action,"K") !=NULL) (*Ke_orig) = (*Ke); //commented out BCO 010810
	//if(strstr(action,"F") !=NULL) for (i=0;i<numDof; i++) F_orig[i] =  Fe[i]; //exactly why was this commented ??
	if(strstr(action,"F") !=NULL) for (i=0;i<numDof; i++) Fe[i] = F_orig[i]; //find out how introducing this will affect the rest of the code
}

//====================================================================
void HeatTransferElement3Dtransient::update(char * action, double *uGlobal)
//====================================================================
{
	//...............................................................
	//    these data are for the element and to be updated
	//...............................................................
	double &elementVolume           = bag->elementVolume;
	double *fluxVolume_element      = bag->fluxVolume_element;
	double *gradientVolume_element  = bag->gradientVolume_element;
	//xtang 022091
	double &concentrationVolume_element= bag->temperatureVolume_element;

	// clear room 
	elementVolume =0.;
	concentrationVolume_element=0;
	int i;

	for( i=0; i<numStrains; i++) {
		fluxVolume_element[i] =0.;  
		gradientVolume_element[i] = 0.;
	}

	//...............................................................
	//    these data are integration calculation only
	//...............................................................
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

    if(ElementNodeRotations == 0){AllocateRotations();}

	//...............................................................
	//    Loop over integration points
	//...............................................................
	loopOverIntegrationPoints(action, idt);
	
	//...............................................................
	//    post process for the update
	//...............................................................
	if(strstr(action,"S") !=NULL )
	{	
		bag->totalVolume += elementVolume;
		bag->temperatureVolume += concentrationVolume_element;   //
		for(i=0; i<3; i++){
			 bag->fluxVolume_total[i] += fluxVolume_element[i];
			 bag->gradientVolume_total[i] += gradientVolume_element[i];
		}
	}

    if(strstr(action,"K") !=NULL) {
        //checkForNegativeDiagonal(); 
    }
}
//====================================================================
void HeatTransferElement3Dtransient::getVolAvgValues(double *uGlobal, double *volAvgFlux,
	         double *volAvgGradient, double *SED, double *volume)
{
// Calculate volume averaged values (gradient and flux) and volume
// and total concentration (stored in SED)
	double *gradientVolume_element  = bag->gradientVolume_element;
	double *fluxVolume_element      = bag->fluxVolume_element;
	update("FS", uGlobal);

	*volume = bag->elementVolume;
	SED[0]=bag->temperatureVolume_element;

	for(int is=0; is<numStrains; is++){
	   volAvgFlux[is]	 = fluxVolume_element[is];
	   volAvgGradient[is]= gradientVolume_element[is];
	}
}
//====================================================================
void HeatTransferElement3Dtransient::PrintActualConcentrations(double *uGlobal, ostream &SO)
{
	int i;
	IntegrationDataType   idt;
	idt.id=new InterpAndSolDerivatives(numNodesPerElement,numNodesPerElement); idt.id->allocate();

	idt.elemDisp[0] = idt.id->nodalValues_u1; 
	idt.elemDisp[1] = idt.id->nodalValues_u2;
	idt.elemDisp[2] = idt.id->nodalValues_u3;     
	idt.id->numberOfInterp = numNodesPerElement;
	idt.uGlobal=uGlobal;
	if( uGlobal != NULL) {extractSolution( uGlobal, idt.elemDisp);} 

	double saturationValue = ((HeatTransferMaterial*)material)->getSaturationValue();

	for(i=0; i<numNodesPerElement; i++){
		SO << idt.elemDisp[0][i] * saturationValue << endl;
	}
}
//======================================================================
void HeatTransferElement3Dtransient::PrintQuadFluxes(ostream &SO,int totalNumIPs)
{
	#define FORMAT setw(18)<<setprecision(12)
	int numFluxes = numStrains;
	
	SO<<elementNumber<<'\t'<<materialGroup<<'\t'<<totalNumIPs<<'\t'<<numStrains<<endl;
	for(int i=0;i<totalNumIPs;i++){
		for(int j=0;j<numFluxes;j++){
			SO<<FORMAT<<quadPointFluxes[i][j]<<'\t';
		}
		SO<<endl;
	}
}
//======================================================================
void HeatTransferElement3Dtransient::allocateElementalDataMemory(int WorkspaceMaxNumberDof)
{
Ke_orig = new Matrix(WorkspaceMaxNumberDof,WorkspaceMaxNumberDof);
Me = new Matrix(WorkspaceMaxNumberDof,WorkspaceMaxNumberDof);

BasicElement::allocateElementalDataMemory(WorkspaceMaxNumberDof);
}
/*
void HeatTransferElement1Dtransient::extractSolution( double * uGlobal, double ** elementDisp)
{
	//this is to convert negative values in the solution to zero
	int i;

   if(uGlobal !=NULL) 
     {
      int firstDof;
      int numberOfDof;
      for(i=0 ; i<numNodesPerElement;i++) {
	      firstDof = node[i].getFirstDof();
         numberOfDof = node[i].getNumDof();              
         switch (numberOfDof){
            case 1:
		{
			if(uGlobal[firstDof] < 0) 
				elementDisp[0][i] = 0.0;
			else
				elementDisp[0][i] = uGlobal[firstDof];   
			break;
		}
            case 2:
		{elementDisp[0][i] = uGlobal[firstDof];
		 elementDisp[1][i] = uGlobal[++firstDof]; break;}
            case 3:
		{elementDisp[0][i] = uGlobal[firstDof];
		 elementDisp[1][i] = uGlobal[++firstDof];
       elementDisp[2][i] = uGlobal[++firstDof];
			 break;}
            default:
                {BETA_OUT<<"Case not implemented yet"<<endl;
                 exit(1);
                }
            }//end of switch
	  }
     }
}
*/

//Static data

int HeatTransferElement3Dtransient::effectivePropertyCase;
char HeatTransferElement3Dtransient::configDescription[80];

//==========================================================
