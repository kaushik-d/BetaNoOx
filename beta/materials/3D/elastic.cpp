#include "stdafx.h"
#include "utility/utility.h"
#include "defines/beta_core_defines.h"

#include "elastic.hpp"
#include "math/matrix.hpp"
#include "math/tensor.hpp"
#include "materials/ElasticISV.hpp"

//void rotateEngineeringVoigtStrain(double *rotatedStress, double *stress,
//                       const double &theta,const int idir);

//void rotateVoigtMaterial(Matrix &Cr, Matrix &cMat,const double &theta, 
//	const int idir);
//void rotateVoigtStress(double *rotatedStress, double *stress,
//                       double theta,int idir);


#define LINE "--------------------------------------------------"
#define FORMAT setw(15)<<setprecision(4)

//==============================================================
/*
void ElasticMaterial::rotateStress(double *newStress,
				                      double *oldStress,
                                  double angle,
											 int rotationAxis, bool G2M)
{ 
	double *dd=new double [numStrains];
	int i;
	for (i=0; i<numStrains; i++) {dd[i]=oldStress[i];}

	if (G2M) 
	{
 		if(xtNumAngles>0) rotateVoigtStress(dd,oldStress,-xtAngles[0], xtAxis[0]);
	} else { // xtang 0132003: Need to check the order of transformation. Reverse order of cMat?
 		if(xtNumAngles>0) rotateVoigtStress(dd,oldStress,xtAngles[0], xtAxis[0]);
	}

//	rotateVoigtStress(newStress,dd,angle,rotationAxis);	
//JV052104 - dirty fix - this does not correctly rotate 2d stresses/strains
	if(numStrains==3){
//		BETA_OUT << " *** Not actually rotating stress/strain !! -see code for more details\n";
		for (i=0; i<numStrains; i++) {newStress[i]=dd[i];}
	}
	else{
		//rotateEngineeringVoigtStrain should not be used. don't know how it got there in the first place
		//but does not look like you need to account for 'engineering' factors here. 
		//rotateEngineeringVoigtStrain(newStress,dd,angle,rotationAxis);	

		rotateVoigtStress(newStress,dd,angle,rotationAxis);	
	}

	delete [] dd;
}
*/
//==============================================================
//JV031709 : new version fixing xtangles
//void ElasticMaterial::rotateStress(double *newStress,
//									double *oldStress,
//									double angle,
//									int rotationAxis, bool G2M)
//{ 
//	if(G2M){
//		rotateStressFromGlobalToMaterial(newStress, oldStress, angle, rotationAxis);
//	}else{
//		rotateStressFromMaterialToGlobal(newStress, oldStress, angle, rotationAxis);
//	}
//}
//==============================================================
//void ElasticMaterial::rotateStressFromGlobalToMaterial(double *newStress,
//									double *oldStress,
//									double angle,
//									int rotationAxis)
//{ //first mangles, then xtAngles
//	double *dd=new double [numStrains];
//	int i;
//
////JV052104 - dirty fix - this does not correctly rotate 2d stresses/strains
//	if(numStrains==3){
////		BETA_OUT << " *** Not actually rotating stress/strain !! -see code for more details\n";
//		for (i=0; i<numStrains; i++) {dd[i]=oldStress[i];}
//	}
//	else{
//		//rotateEngineeringVoigtStrain should not be used. don't know how it got there in the first place
//		//but does not look like you need to account for 'engineering' factors here. 
//		//rotateEngineeringVoigtStrain(newStress,dd,angle,rotationAxis);	
//		rotateVoigtStress(dd,oldStress,angle,rotationAxis);	
//	}
//
//	if(xtNumAngles>0) rotateVoigtStress(newStress, dd, -xtAngles[0], xtAxis[0]);
//	else {
//		for (i=0; i<numStrains; i++) {newStress[i]=dd[i];}
//	}
//	delete [] dd;
//}
//==============================================================
//void ElasticMaterial::rotateStressFromMaterialToGlobal(double *newStress,
//												double *oldStress,
//												double angle,
//												int rotationAxis)
//{ 
//	double *dd=new double [numStrains];
//	int i;
//	for (i=0; i<numStrains; i++) {dd[i]=oldStress[i];}
//
//	if(xtNumAngles>0) rotateVoigtStress(dd,oldStress,xtAngles[0], xtAxis[0]);
//
////JV052104 - dirty fix - this does not correctly rotate 2d stresses/strains
//	if(numStrains==3){
////		BETA_OUT << " *** Not actually rotating stress/strain !! -see code for more details\n";
//		for (i=0; i<numStrains; i++) {newStress[i]=dd[i];}
//	}
//	else{
//		//rotateEngineeringVoigtStrain should not be used. don't know how it got there in the first place
//		//but does not look like you need to account for 'engineering' factors here. 
//		//rotateEngineeringVoigtStrain(newStress,dd,angle,rotationAxis);	
//		rotateVoigtStress(newStress,dd,angle,rotationAxis);	
//	}
//
//	delete [] dd;
//}
//==============================================================
/* //does not work when using xtangles AND mangles
void ElasticMaterial::rotateStrain(double *newStress,
									double *oldStress,
									double angle,
									int rotationAxis)
{ 
	// xtang 050102: 
	int i;
	double *dd=new double [numStrains];
	for (i=0; i<numStrains; i++) {dd[i]=oldStress[i];}
	if (xtNumAngles > 0) 
	{
		rotateEngineeringVoigtStrain(dd,oldStress,-xtAngles[0], xtAxis[0]);
	}
	//rotateEngineeringVoigtStrain(newStress,dd,angle,rotationAxis);	
//JV052104 - dirty fix - this does not correctly rotate 2d stresses/strains
	if(numStrains==3){
//		BETA_OUT << " *** Not actually rotating stress/strain !! -see code for more details\n";
		for (i=0; i<numStrains; i++) {newStress[i]=dd[i];}
	}
	else{
		rotateEngineeringVoigtStrain(newStress,dd,angle,rotationAxis);	
	}
	delete [] dd;
}
*/
//---------------------------------------------------------------------------
//void ElasticMaterial::rotateStrain(double *newStress,
//									double *oldStress,
//									double angle,
//									int rotationAxis)
//{ 
////JV031709 : to rotate from Global to material coordinate system
////first mangles, then xtAngles
//	double *dd=new double [numStrains];
//	int i;
//
////JV052104 - dirty fix - this does not correctly rotate 2d stresses/strains
//	if(numStrains==3){
////		BETA_OUT << " *** Not actually rotating stress/strain !! -see code for more details\n";
//		for (i=0; i<numStrains; i++) {dd[i]=oldStress[i];}
//	}
//	else{
//		rotateEngineeringVoigtStrain(dd,oldStress,angle,rotationAxis);	
//	}
//
//	if(xtNumAngles>0) rotateEngineeringVoigtStrain(newStress, dd, -xtAngles[0], xtAxis[0]);
//	else {
//		for (i=0; i<numStrains; i++) {newStress[i]=dd[i];}
//	}
//	delete [] dd;
//}
//---------------------------------------------------------------------------
//need to verify if the following function is doing the correct transformation (check order)
//void ElasticMaterial::rotateEngineeringStrain(double *newStrain,
//				                      double *oldStrain,
//                                  double angle,
//											 int rotationAxis, bool G2M=true)
//{
//	double *dd=new double [numStrains];
//	EXIT_BETA("This function: ElasticMaterial::rotateEngineeringStrain() needs to be fixed before it can be used!!!")
//
//	//JV052104 - this rotateVoigtStress function is messed up cos it does not 
//	// differentiate between 2d and 3d. need to fix !!!
//
//	for (int i=0; i<numStrains; i++) {dd[i]=oldStrain[i];}
//	if (G2M) {
//		if(xtNumAngles>0) rotateEngineeringVoigtStrain(dd,oldStrain,-xtAngles[0], xtAxis[0]);
//		rotateEngineeringVoigtStrain(newStrain,dd,angle,rotationAxis);	
//	} else {  // xtang 0132003: Need to check the order of transformation. Reverse order of cMat?
//		if(xtNumAngles>0) rotateEngineeringVoigtStrain(dd,oldStrain,xtAngles[0], xtAxis[0]);
//		rotateEngineeringVoigtStrain(newStrain,dd,angle,rotationAxis);	
//	}
//	delete [] dd;
//}
//---------------------------------------------------------------------------
ElasticMaterial::ElasticMaterial()
{
	cMat=0;
	S=0;
	globalCmat=0;
	rotatedCmat=0;
	localizedCmat=0;
}
//--------------------------------------------------------------------
ElasticMaterial::ElasticMaterial(const ElasticMaterial& arg):
Material(arg)
{
    numberOfRotationAngles=arg.numberOfRotationAngles;   //max of 4 angles
	for(int i=0;i<4;i++){rotationAxis[i]=arg.rotationAxis[i];}
	for(int i=0;i<4;i++){rotationAngle[i]=arg.rotationAngle[i];}
    angleZ=arg.angleZ;
	currentRotationAngle=arg.currentRotationAngle;
	numStrains=arg.numStrains;
	FlagSingleAngle=arg.FlagSingleAngle;

	cMat			= new Matrix(*arg.cMat);
	S				= new Matrix(*arg.S);
	rotatedCmat		= new Matrix(*arg.rotatedCmat);
	globalCmat		= new Matrix(*arg.globalCmat);
	localizedCmat	= new Matrix(*arg.localizedCmat);

	materialState=arg.materialState;
	for(int i=0;i<6;i++){alpha[i]=arg.alpha[i];}
	for(int i=0;i<6;i++){beta[i]=arg.beta[i];}

	//xtNumAngles=arg.xtNumAngles;
	//for(int i=0;i<4;i++){xtAxis[i]=arg.xtAxis[i];}
	//for(int i=0;i<4;i++){xtAngles[i]=arg.xtAngles[i];}
}
//--------------------------------------------------------------------
ElasticMaterial::~ElasticMaterial(void)
{
	if(cMat) delete cMat;
	if(S) delete S;
	if(globalCmat) delete globalCmat;
	if(rotatedCmat) delete rotatedCmat;
	if(localizedCmat) delete localizedCmat;
}
//--------------------------------------------------------------------
void ElasticMaterial::initialize()
{
	numStrains    = 6;  // default xtang 980329
	materialState = clear;
	CommonElasticMaterial();

	// xtang 050102
	//xtNumAngles=0;
	//xtAxis[0]=xtAxis[1]=xtAxis[2]=xtAxis[3]=0;   //valid number should 1/2/3 only
	//xtAngles[0]=xtAngles[1]=xtAngles[2]=xtAngles[3]=0.0;
	//numberOfRotationAngles=0;  // must be initialized...xtang 073002

	setType("elastic");
//		 BETA_OUT<<"I'm here in file elastic.cpp"<<endl;
}
//--------------------------------------------------------------------
void ElasticMaterial::allocateISVs(std::vector<BasicISV*> &isv, const int &numISVs)const
{
    //For this case, allocate ISVs out of a continuous block of memory...
    //isv[0] points to the allocated block
    if(!isv.empty()){
        deallocateISVs(isv);
    }
    isv.resize(numISVs);
    isv[0] = new ElasticISV[numISVs];
    ElasticISV* start = (ElasticISV*)isv[0];
    for(int i=1;i<numISVs;++i){
        isv[i]=&(start[i]);
    }
}
//===========================================================
void ElasticMaterial::deallocateISVs(std::vector<BasicISV*> &isv) const 
{
    if(isv[0]){
        delete [] (ElasticISV*)isv[0];
        isv.assign(isv.size(),NULL);
    }
}
//===========================================================
void ElasticMaterial::CommonElasticMaterial(void)
{
	cMat		= new Matrix(numStrains,numStrains);
	S			= new Matrix(numStrains,numStrains);
	globalCmat	= new Matrix(numStrains,numStrains);
	rotatedCmat	= new Matrix(numStrains,numStrains);
	localizedCmat	= new Matrix(numStrains,numStrains);

	// xtang 050102
	//xtNumAngles=0;
	//xtAxis[0]=xtAxis[1]=xtAxis[2]=xtAxis[3]=0;   //valid number should 1/2/3 only
	//xtAngles[0]=xtAngles[1]=xtAngles[2]=xtAngles[3]=0.0;

    //numberOfRotationAngles=0;  // must be initialized...xtang 073002 

	for (int i=0; i<6; i++) {alpha[i]=beta[i]=0.0;}

}
//==============================================================
// Order of elastic properties...
//    0->3         E11          E22          E33           NU12
//	   4->7        NU23         NU13          G12            G23
//    8            G13 

//       ANGLE       ALPHA1         ALPHA2
///   12->15    ALPHA3        BETA1        BETA2          BETA3

//Note: NUij = -eps_j/eps_i   where load is applied in x_i direction

bool ElasticMaterial::read(istream * inStream)
{
	int i;
	double E[16];
	char *token;
	const char comment[]= "//" ;	
	int  numberOfTokens;
	char *tokenList[20];	// XTANG 980329
    double nu21, nu32, nu31;
   	ios_base:: fmtflags oldIOS;
	int oldPrecision;

    alpha[0] = alpha[1] = alpha[2] = alpha[3] = alpha[4] = alpha[5] = 0.;
    beta[0] = beta[1] = beta[2] = beta[3] = beta[4] = beta[5] = 0.;

	while(2==2) { 
		if(getLineAndTokenize(inStream,"exitElasticMaterial",tokenList, numberOfTokens)  ==1){
			BETA_OUT<<"End of input for this material"<<LINE<<"\n\n"<<endl;
			return true; 
		}
		token = tokenList[0];
		banner(token, BETA_OUT);	

		if(COMPARE(token,"readModuli")==0 )
		{
			for(i=0;i<9;i++){ (*inStream)>>E[i];}
			nu21 = E[3] * E[1]/E[0] ;
			nu32 = E[4] * E[2]/E[1] ;
			nu31 = E[5] * E[2]/E[0] ;

			BETA_OUT<<"Input data"<<endl;		 
			oldIOS = BETA_OUT.setf(ios::scientific);
			oldPrecision = BETA_OUT.precision(4);
			BETA_OUT << "\tE11 = " << E[0] << "\tnu12 = " << E[3] <<"\tG12 = " << E[6] <<endl;
			BETA_OUT << "\tE22 = " << E[1] << "\tnu23 = " << E[4] <<"\tG23 = " << E[7] <<endl;
			BETA_OUT << "\tE33 = " << E[2] << "\tnu13 = " << E[5] <<"\tG13 = " << E[8] <<endl; 

			BETA_OUT << "-------------------" << endl;


			BETA_OUT<<"For your information, the following is computed"<<endl;
			BETA_OUT<<"NU21, NU32, NU31   = "<<nu21<<"  "<<nu32<<"  "<<nu31<<endl;
			BETA_OUT.setf(oldIOS);
			BETA_OUT.precision(oldPrecision);
		 
		
			(*S)(0,0) = 1/E[0];						// E[0] = E11
			(*S)(0,1) = (*S)(1,0) = -E[3]/E[0];     // E[3] = nu12
			(*S)(0,2) = (*S)(2,0) = -E[5]/E[0];     // E[4] = nu23
			(*S)(1,1) = 1/E[1];						// E[1] = E22
			(*S)(2,1) = (*S)(1,2) = -E[4]/E[1];	    // E[5] = nu13
			(*S)(2,2) = 1/E[2];						// E[2] = E33
			(*S)(3,3) = 1/E[6];						// E[6] = G12
			(*S)(4,4) = 1/E[7];						// E[7] = G23
			(*S)(5,5) = 1/E[8];						// E[8] = G13
			materialState |= haveS;
			(*cMat) = (*S);
			cMat->invert();
			materialState |= haveC;
			//print(BETA_OUT);
			cMat->print("Constitutive matrix: local coordinate system",&BETA_OUT);
			(*globalCmat) = (*cMat);  //check this !!

			//Should add check for transverse isotropy
		} //end of read moduli

        if(COMPARE(token,"readAngles")==0 || COMPARE(token,"readTransformationAngles")==0){
            // This code creates a rotation object that will be applied to all instances of this material after mangles is applied.
            if(COMPARE(token,"readTransformationAngles")==0){
                cout << "WARNING - readTransformationAngles is depreciated and will be removed." << endl;
                cout << "Use readAngles instead." << endl;
            }

            int NumAngles, axis;
            double angle;
            (*inStream) >> NumAngles;
            BETA_OUT<<"Specified " << NumAngles << " sequential material rotations"<<endl;
            for(int i=0;i<NumAngles;++i){
                (*inStream) >> axis >> angle;
                BETA_OUT << "Rotating about axis " << axis << " by angle " << angle << endl;
                // General Note about Rotation Class and Angles:
                // Angles are input into Beta in terms of rotating a material about
                // a coordinate system.  The rotations class is defined in terms of 
                // rotating a coordinate system about a material.  One is the opposite
                // of the other, so when constructing rotation objects, input the opposite
                // of the angle that is used to input the rotation into Beta.
                MaterialRotation.add_after_self(Rotation(axis,-angle));
            }
        }

		if(COMPARE(token,"readThermalExpansionCoefficients")==0){
			(*inStream)>>alpha[0] >>alpha[1] >> alpha[2];
			BETA_OUT<<FORMAT<<alpha[0]<<FORMAT<<alpha[1]<<FORMAT<<alpha[2]<<endl;
		}

	} //end of while 
}

//==========================================================================
void ElasticMaterial::calculateStress(double *stress,
                                      double *strain,
                                      BasicISV* isv,
                                      const Rotation* Orientation)
{
    //input strain in GCS
    //output stress in GCS
	int i,j;
	double *C_i;

    double HygroThermStrainGCS[] = {0.0,0.0,0.0,0.0,0.0,0.0}, 
           HygroThermStrainLCS[] = {0.0,0.0,0.0,0.0,0.0,0.0};
	
    ElasticISV* Eisv = (ElasticISV*) isv;
    // Get the Hygrothermal strain in the local coordinate systems
	HygroThermStrainLCS[0]=alpha[0] * Eisv->getTemperature() + beta[0] * Eisv->getMoisture();
    HygroThermStrainLCS[1]=alpha[1] * Eisv->getTemperature() + beta[1] * Eisv->getMoisture();
	HygroThermStrainLCS[2]=alpha[2] * Eisv->getTemperature() + beta[2] * Eisv->getMoisture();
    HygroThermStrainLCS[3]=alpha[3] * Eisv->getTemperature() + beta[3] * Eisv->getMoisture();
    HygroThermStrainLCS[4]=alpha[4] * Eisv->getTemperature() + beta[4] * Eisv->getMoisture();
    HygroThermStrainLCS[5]=alpha[5] * Eisv->getTemperature() + beta[5] * Eisv->getMoisture();

    // Rotate the Hygrothermal strain to the global coordinate system
    Orientation->rotate_Voigt_engineering_strain_from_local_to_global(HygroThermStrainLCS,HygroThermStrainGCS);

	//########################################################################
	//Now we have Hygrothermal Strain in GCS
	//and strain due to our Mechanical loading in GCS

	//########################################################################
	//This is calculating the Total(Mechanical+Thermal) stress in GCS
	for(i=0;i<numStrains;i++) 
	{
		stress[i]=0;
		C_i = (*globalCmat)[i];
		for(j=0;j<numStrains;j++) 
		{
			stress[i] += C_i[j] * ( strain[j] - HygroThermStrainGCS[j]);
		}
	}//Stress in GCS is done being calculated 
}

//==========================================================================
void ElasticMaterial::calculateInitialStress(double *stress,
                                             BasicISV* isv,
                                             const Rotation* Orientation)
{
    // Look at CalculateStress (the one with rotations) to see what this is doing.
	int i,j;
	double *C_i;

    double HygroThermStrainGCS[] = {0.0,0.0,0.0,0.0,0.0,0.0}, 
           HygroThermStrainLCS[] = {0.0,0.0,0.0,0.0,0.0,0.0};
	
    ElasticISV* Eisv = (ElasticISV*) isv;
    // Get the Hygrothermal strain in the local coordinate systems
	HygroThermStrainLCS[0]=alpha[0] * Eisv->getTemperature() + beta[0] * Eisv->getMoisture();
    HygroThermStrainLCS[1]=alpha[1] * Eisv->getTemperature() + beta[1] * Eisv->getMoisture();
	HygroThermStrainLCS[2]=alpha[2] * Eisv->getTemperature() + beta[2] * Eisv->getMoisture();
    HygroThermStrainLCS[3]=alpha[3] * Eisv->getTemperature() + beta[3] * Eisv->getMoisture();
    HygroThermStrainLCS[4]=alpha[4] * Eisv->getTemperature() + beta[4] * Eisv->getMoisture();
    HygroThermStrainLCS[5]=alpha[5] * Eisv->getTemperature() + beta[5] * Eisv->getMoisture();

    // Rotate the Hygrothermal strain to the global coordinate system
    Orientation->rotate_Voigt_engineering_strain_from_local_to_global(HygroThermStrainLCS,HygroThermStrainGCS);

	//########################################################################
	//Now we have Hygrothermal Strain in GCS
	//and strain due to our Mechanical loading in GCS

	for(i=0;i<numStrains;i++) {
		stress[i]=0;
		C_i = (*globalCmat)[i];
		for(j=0;j<numStrains;j++) {
            //Note - this SHOULD be positive.  It will eventually be ADDED to the pseudoloadvector in equation::addmatrix()
	          stress[i] += C_i[j] * HygroThermStrainGCS[j];
		}
	}
}	
//==========================================================================
void ElasticMaterial::print(ostream &ostrm)
{
	double E11,E22,E33,nu12,nu13,nu23,G12,G23,G13;
	//long oldIOS;
	ios_base:: fmtflags oldIOS;
	int oldPrecision;

// See earlier printing...consolidate !!!!

	if(!(materialState & haveS)) {
		if((materialState & haveC)) {
			(*S) = (*cMat);
			(*S).invert();
			materialState |= haveC;
		}
		else {
			BETA_OUT<<"No material properties to print..."<<endl;
			return;
		}
	}
		
	ostrm << "-------------------\nElastic Material: " <<
			getGroupName() << "\nGroup Number: " << getGroupNum() << '\n';
	E11 = 1/(*S)[0][0]; 
	E22 = 1/(*S)[1][1]; nu12 = -(*S)[0][1]*E11;
	E33 = 1/(*S)[2][2]; nu13 = -(*S)[0][2]*E11; nu23 = -(*S)[2][1]*E22;

	G12 = 1/(*S)[3][3];
	G23 = 1/(*S)[4][4];
	G13 = 1/(*S)[5][5];
	oldIOS = ostrm.setf(ios::scientific);
	oldPrecision = ostrm.precision(4);
	ostrm << "\tE11 = " << E11 << "\tnu12 = " << nu12 << "\tnu13 = " << nu13 << '\n';
	ostrm << "\tE22 = " << E22 << "\tnu23 = " << nu23 << '\n';
	ostrm << "\tE33 = " << E33 << '\n';
	ostrm << "\tG12 = " << G12 << "\tG23 = " << G23 << "\tG13 = " << G13 << '\n';
	ostrm << "-------------------" << endl;
	ostrm.setf(oldIOS);
	ostrm.precision(oldPrecision);
}

//==========================================================================
void ElasticMaterial::rotateMaterialToGlobal(const Rotation* OrientationAtPoint)
{
    OrientationAtPoint->rotate_Cmat_from_local_to_global(*cMat,*globalCmat);
}

//==========================================================================
//void rotateVoigtStress(double *rotatedStress, double *stress,
//                       double theta,int idir)
//{
// // newStress = T * oldStress
//
//// What is the sign convention ?   What is order of stresses?????
//
//if(theta == 0.  || idir == 0)
//{
// for(int j=0; j<6; j++)rotatedStress[j] = stress[j];
// return;
//}
//
//	RotationTensor a(theta,idir);  //  This is likely slow!!!
//	Matrix t(6,6),ttc(6,6);       
//	int i,k;
// 
//	t=0;
//	ttc=0;
//	t(0,0) = a(1,1) * a(1,1);
//	t(1,0) = a(1,2) * a(1,2);
//	t(2,0) = a(1,3) * a(1,3);
//	t(3,0) = a(1,1) * a(1,2);
//	t(4,0) = a(1,2) * a(1,3);
//	t(5,0) = a(1,1) * a(1,3);
// 
//	t(0,1) = a(2,1) * a(2,1);
//	t(1,1) = a(2,2) * a(2,2);
//	t(2,1) = a(2,3) * a(2,3);
//	t(3,1) = a(2,1) * a(2,2);
//	t(4,1) = a(2,2) * a(2,3);
//	t(5,1) = a(2,1) * a(2,3);
//
//	t(0,2) = a(3,1) * a(3,1);
//	t(1,2) = a(3,2) * a(3,2);
//	t(2,2) = a(3,3) * a(3,3);
//	t(3,2) = a(3,1) * a(3,2);
//	t(4,2) = a(3,2) * a(3,3);
//	t(5,2) = a(3,1) * a(3,3);
//
//	t(0,3) = 2. * a(1,1) * a(2,1);
//	t(1,3) = 2. * a(1,2) * a(2,2);
//	t(2,3) = 2. * a(1,3) * a(2,3);
//	t(3,3) = a(1,1) * a(2,2) + a(1,2) * a(2,1);
//	t(4,3) = a(1,2) * a(2,3) + a(1,3) * a(2,2);
//	t(5,3) = a(1,1) * a(2,3) + a(1,3) * a(2,1);
//
//	t(0,4) = 2. * a(2,1) * a(3,1);
//	t(1,4) = 2. * a(2,2) * a(3,2);
//	t(2,4) = 2. * a(2,3) * a(3,3);
//	t(3,4) = a(2,1) * a(3,2) + a(2,2) * a(3,1);
//	t(4,4) = a(2,2) * a(3,3) + a(2,3) * a(3,2);
//	t(5,4) = a(2,1) * a(3,3) + a(2,3) * a(3,1);
//
//	t(0,5) = 2. * a(1,1) * a(3,1);
//	t(1,5) = 2. * a(1,2) * a(3,2);
//	t(2,5) = 2. * a(1,3) * a(3,3);
//	t(3,5) = a(1,1) * a(3,2) + a(1,2) * a(3,1);
//	t(4,5) = a(1,2) * a(3,3) + a(1,3) * a(3,2);
//	t(5,5) = a(1,1) * a(3,3) + a(1,3) * a(3,1);
//  
//	// Perform transformation  ( new = T * old )
//
//		register double sum;
//		for(i=0;i<6;i++) {
//				sum = 0;
//				for(k=0;k<6;k++) {sum += t[i][k] * stress[k];}
//		 rotatedStress[i] = sum;}
//
//}
//==========================================================================
//void rotateEngineeringVoigtStrain(double *rotatedStress, double *stress,
//                       const double &theta,const int idir)
//{
////BETA_OUT<<"input for rotateEngineeringVoigtStrain"<<endl;
////   	 BETA_OUT<<stress[0]<<" "<<stress[1]<<"  "<<stress[2]<<endl;   	                      
////    	 BETA_OUT<<stress[3]<<" "<<stress[4]<<"  "<<stress[5]<<endl;
//    	 
// stress[3] = stress[3]/2.;
// stress[4] = stress[4]/2.;
// stress[5] = stress[5]/2.;
// 
// rotateVoigtStress(rotatedStress, stress,theta,idir);  ///xtang
// 
// stress[3] = 2. * stress[3];
// stress[4] = 2. * stress[4];
// stress[5] = 2. * stress[5];                  
// rotatedStress[3] = 2. * rotatedStress[3];
// rotatedStress[4] = 2. * rotatedStress[4];
// rotatedStress[5] = 2. * rotatedStress[5]; 
// 
//// BETA_OUT<<"output for rotateEngineeringVoigtStrain"<<endl;
////   	 BETA_OUT<<rotatedStress[0]<<" "<<rotatedStress[1]<<"  "<<rotatedStress[2]<<endl;   	                      
////    	 BETA_OUT<<rotatedStress[3]<<" "<<rotatedStress[4]<<"  "<<rotatedStress[5]<<endl; 
//}


//==========================================================================
//void ElasticMaterial::rotateMaterialToGlobalSingleAngleXYZ(double *singleAngleXYZ)
//{
////Note that this is set up for one rotation only
//// Change this to use a saved transformation matrix ...
//	//same for stresses
//// cout<<"FlagSingleAngle="<<getFlagSingleAngle()<<endl;
//Matrix *Tmp1, *Tmp2;
//Tmp1 = new Matrix(6,6);
//Tmp2 = new Matrix(6,6);
////cMat->print("cMat, local"); // Jae, test 7/22/01
//
// if(fabs(singleAngleXYZ[0])<1.0e-6){(*Tmp1) = (*cMat);}
// else {rotateVoigtMaterial( (*Tmp1),(*cMat),-1.0*singleAngleXYZ[0],1);}
//
// if(fabs(singleAngleXYZ[1])<1.0e-6){(*Tmp2) = (*Tmp1);}
// else {rotateVoigtMaterial( (*Tmp2),(*Tmp1),-1.0*singleAngleXYZ[1],2);}
//
//// cout<<"theta= "<<-1.0*singleAngleXYZ[1]<<"  idir="<<2<<endl;
////Tmp2->print("cMat rotated1");
//
// if(fabs(singleAngleXYZ[2])<1.0e-6){(*globalCmat) = (*Tmp2); }
// else {rotateVoigtMaterial( (*globalCmat),(*Tmp2),-1.0*singleAngleXYZ[2],3);}
//
// //cout<<"theta= "<<-1.0*singleAngleXYZ[2]<<"  idir="<<3<<endl;
// //globalCmat->print("cMat rotated2");
////exit(1);
//
// delete Tmp1;
// delete Tmp2;
// return;
//}
//
//
////==========================================================================
//
//#define LINE "--------------------------------------------------"
//#define FORMAT setw(15)<<setprecision(4)
////==============================================================
//void ElasticMaterial::rotateStressSingleAngleXYZ(double *newStress,
//				         double *oldStress, double *singleAngleXYZ)
//{ 
//	double *Tmp1, *Tmp2;
//	Tmp1 = new double [numStrains];
//	Tmp2 = new double [numStrains];
//	if(fabs(singleAngleXYZ[2])<1.0e-6)
//	{	for(int i=0; i<numStrains; i++)Tmp1[i] = oldStress[i];  }
//    else {rotateVoigtStress( Tmp1,oldStress,-singleAngleXYZ[2],3);}
//
//	if(fabs(singleAngleXYZ[1])<1.0e-6)
//	{	for(int i=0; i<numStrains; i++)Tmp2[i] = Tmp1[i]; 	}
//    else {rotateVoigtStress( Tmp2,Tmp1,-singleAngleXYZ[1],2);}
//
//	if(fabs(singleAngleXYZ[0])<1.0e-6)
//	{	for(int i=0; i<numStrains; i++)newStress[i] = Tmp2[i];  }
//    else {rotateVoigtStress( newStress,Tmp2,-singleAngleXYZ[0],1);}
//	delete Tmp1;
//	delete Tmp2;
//	return;
//		//for(int i=0; i<numStrains; i++)newStress[i] = oldStress[i];
//}
////==============================================================
//void ElasticMaterial::rotateStrainSingleAngleXYZ(double *newStress,
//				         double *oldStress, double *singleAngleXYZ)
//{ 
//
//	double *Tmp1, *Tmp2;
//	Tmp1 = new double [numStrains];
//	Tmp2 = new double [numStrains];
//
//	if(fabs(singleAngleXYZ[2])<1.0e-6)
//	{	for(int i=0; i<numStrains; i++)Tmp1[i] = oldStress[i]; 	}
//    else {rotateEngineeringVoigtStrain( Tmp1,oldStress,-singleAngleXYZ[2],3);}
//
//	if(fabs(singleAngleXYZ[1])<1.0e-6)
//	{	for(int i=0; i<numStrains; i++)Tmp2[i] = Tmp1[i];  }
//    else {rotateEngineeringVoigtStrain( Tmp2,Tmp1,-singleAngleXYZ[1],2);}
//
//	if(fabs(singleAngleXYZ[0])<1.0e-6)
//	{	for(int i=0; i<numStrains; i++)newStress[i] = Tmp2[i]; 	}
//    else {rotateEngineeringVoigtStrain( newStress,Tmp2,-singleAngleXYZ[0],1);}
// 
//	delete Tmp1;
//	delete Tmp2;
//	return;
//
//	//rotateEngineeringVoigtStrain(newStress,oldStress,angle,rotationAxis);	
//}
//#undef FORMAT
