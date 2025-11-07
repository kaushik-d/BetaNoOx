#include "stdafx.h"

#include "utility/utility.h"
#include "defines/beta_core_defines.h"
#include "dmElastic.hpp"

#include "../material/BasicDegradationModel.hpp"

void rotateEngineeringVoigtStrain(double *rotatedStress, double *stress,
                       const double &theta,const int idir);

void rotateVoigtMaterial(Matrix &Cr, Matrix &cMat,const double &theta, 
	const int idir);
void rotateVoigtStress(double *rotatedStress, double *stress,
                       double theta,int idir);


void rotateEngineeringVoigtStrain(double *rotatedStress, double *stress,
                       const double &theta,const int idir);


#define LINE "--------------------------------------------------"
#define FORMAT setw(15)<<setprecision(4)

//==========================================================================
dmElasticMaterial::dmElasticMaterial()
{
	DegradationModel=0;
}
//==========================================================================
dmElasticMaterial::dmElasticMaterial(const dmElasticMaterial& arg):ElasticMaterial(arg)
{
	outOfPlaneDir = arg.outOfPlaneDir;
	numISV = arg.numISV;

	for(int i=0; i<9; i++)
	{
		strength[i] = arg.strength[i];
	}

	DegradationModel = arg.DegradationModel->clone();
}
//==========================================================================
dmElasticMaterial::~dmElasticMaterial(void)
{
	if(DegradationModel) delete DegradationModel;
}
//==========================================================================
bool dmElasticMaterial::read(istream * inStream)
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

    while(2==2) { 
        if(getLineAndTokenize(inStream,"exitElasticMaterial",tokenList,numberOfTokens)  ==1)
        {
            BETA_OUT<<"End of input for this material"<<LINE<<"\n\n"<<endl;
            return true; 
        }
        token = tokenList[0];
        banner(token, BETA_OUT);	

        if(COMPARE(token,"readModuli")==0 ){
            for(i=0;i<9;i++){ (*inStream)>>E[i];}

            nu21 = E[3] * E[1]/E[0] ;
            nu32 = E[4] * E[2]/E[1] ;
            nu31 = E[5] * E[2]/E[0] ;
	    
            BETA_OUT<<"Input data"<<endl;		 
            oldIOS = BETA_OUT.setf(ios::scientific);
            oldPrecision = BETA_OUT.precision(4);
            BETA_OUT<<"E11,  E22,  E33    = "<< E[0]<<"  "<<E[1]<<"  "<<E[2]<<endl;
            BETA_OUT<<"NU12, NU23, NU13   = "<< E[3]<<"  "<<E[4]<<"  "<<E[5]<<endl;
            BETA_OUT<<"G12,  G23,  G13    = "<< E[6]<<"  "<<E[7]<<"  "<<E[8]<<endl;		 
	
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

            // +xtang 990420: set up outOfPlaneDir=> actually is fiber direction
            if ((E[0] >= E[1]) && (E[0] >= E[2])) {outOfPlaneDir=1;}
            if ((E[1] >= E[0]) && (E[1] >= E[2])) {outOfPlaneDir=2;}
            if ((E[2] >= E[0]) && (E[2] >= E[1])) {outOfPlaneDir=3;}

            if ( (fabs((E[0]-E[1])/E[0]) <1e-2) && (fabs((E[0]-E[2])/E[0]) <1e-2) ){
                outOfPlaneDir=0;// isotropic case
            }
            // -xtang

            (*cMat) = (*S);
            cMat->invert();

            materialState |= haveC;
            //print(BETA_OUT);
            cMat->print("Constitutive matrix",&BETA_OUT);
            (*globalCmat) = (*cMat);  //check this !!
            (*localizedCmat) = (*cMat);  //check this !!
      
            //Should add check for transverse isotropy
        } //end of read moduli
    
        if(COMPARE(token,"readAngles")==0 || COMPARE(token,"readTransformationAngles")==0){
            // This code creates a rotation object that will be applied to all instances of this material after mangles is applied.
            if(COMPARE(token,"readTransformationAngles")==0){
                cout << "WARNING - readTransformationAngles is depreciated and will be removed." << endl;
                cout << "Use readAngles instead." << endl;
            }

            double NumAngles;
            double axis, angle;
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

        if(COMPARE(token,"readStrength")==0){
            (*inStream) >> strength[0] >> strength[1]
            >>strength[2]>>strength[3]>>strength[4]
            >>strength[5];

            // assume identical tensile and compressive strength (depends on readCompressiveStrength being called after readStrength)
            strength[6]=-strength[0];
            strength[7]=-strength[1];
            strength[8]=-strength[2];

            BETA_OUT << "Strength Values:" << endl;
            BETA_OUT<<FORMAT<<strength[0]<<FORMAT<<strength[1]<<FORMAT<<strength[2]
            <<FORMAT<<strength[3]<<FORMAT<<strength[4]<<FORMAT<<strength[5]<<endl;
        }

        // xtang 08042003 read compressive strength
        if(COMPARE(token,"readCompressiveStrength")==0){
            (*inStream) >> strength[6] >> strength[7] >> strength[8];

            for(int i=1;i<=3;++i){
                if(strength[i+5] > 0){
                    cout << "Warning in dmElasticMaterial::read(): " << endl;
                    cout << "Compressive " << i << i << " strength input as tensile value, flipping sign automatically." << endl;
                    strength[i+5] = -strength[i+5];
                }
            }
            BETA_OUT << "Compressive Strength Values:" << endl;
            BETA_OUT<<FORMAT<<strength[6]<<FORMAT<<strength[7]<<FORMAT<<strength[8]<<endl;
        }

    } //end of while

}
//==========================================================================
void dmElasticMaterial::initializeDegradationModel(char** localTokenList, int ParamStartIndex, int numberOfTokens)
{
	if(DegradationModel==0)
		EXIT_BETA("Error! calling initializeDegradationModel() before DegradationModel has been created");

	BETA_OUT << "Initializing DegradationModel for material #" << getGroupNum() << " using these parameters: " << localTokenList[ParamStartIndex] << endl;
	
	DegradationModel->InitializeDegradationModel(outOfPlaneDir, localTokenList, ParamStartIndex, numberOfTokens);
}
//==========================================================================
bool dmElasticMaterial::UpdateISV(BasicISV* &isvToDegrade, double* QuadPointStress)
{
    DamageISV* Disv = (DamageISV*) isvToDegrade;
	return DegradationModel->UpdateISVDegradation(Disv, QuadPointStress, strength);
}
//==========================================================================
void dmElasticMaterial::updateLocalizedCmat(Rotation* OrientationAtPoint, BasicISV* aISV)
{
    DamageISV* isv = (DamageISV*) aISV;
    // Temporarily store inverted
	(*localizedCmat)=(*S);

    // Note that at this moment, localizedCmat is actually the compliance
	DegradationModel->updateSMatrix(localizedCmat, isv);
    localizedCmat->invert(); //Inversion in-place

    OrientationAtPoint->rotate_Cmat_from_local_to_global(*localizedCmat,*globalCmat);
}
//==========================================================================
//These next two methods are obsolete due to the way they interface with the ISV object.
//They treat degradation differently if some material is under compression, in consideration
//that cracks will be closed and thus compressive stiffness will be unaffected.

//The methods will remain here until it is determined they will not be needed, or their
//functionality is incorporated into a new type of DegradationModel.

//void dmElasticMaterial::updateLocalizedCmat( Rotation* OrientationAtPoint, BasicISV* aISV, double *strain)
//{
// //Note that this is set up for one rotation only
//    DamageISV* isv = (DamageISV*) aISV;
//
//// Change this to use a saved transformation matrix ...
//	//same for stresses
//
//	int    ii,jj;
//
//	Matrix *mat;
//    mat  = new Matrix(numStrains,numStrains);
//	(*mat)=(*S);
//
//    if (isv->inTension11) (*mat)[0][0]/=isv->E11dam;
//    if (isv->inTension22) (*mat)[1][1]/=isv->E22dam;
//    if (isv->inTension33) (*mat)[2][2]/=isv->E33dam;
//    (*mat)[3][3]/=isv->G12dam;
//    (*mat)[4][4]/=isv->G23dam;
//    (*mat)[5][5]/=isv->G13dam;
//
//	 mat->invert();
//
//	 double stress[6], sum, newstrain[6];
//	 bool moduliChange=false;
//
//	// xtang 050102
//	// strain: in GCS     
//	double *dd=new double [numStrains];
//	for (int i=0; i<numStrains; i++) {dd[i]=strain[i];}
//    //Strains should be rotated through mangles, then xtAngles
//	rotateEngineeringVoigtStrain(newstrain,dd,-theta,dir); ///xtang
//	if (xtNumAngles > 0) // Caution:  assuming from G to L	
//	{
//		rotateEngineeringVoigtStrain(dd,strain,-xtAngles[0], xtAxis[0]);
//	}
//	delete [] dd;
//
//
//	 for (ii=0; ii<3; ii++) {
//		 for(sum=0, jj=0; jj<3; jj++) sum+=(*mat)[ii][jj]*newstrain[jj]; // strain dosen't include eigen strain!!!
//		 stress[ii]=sum;
//	 }
//
//	 // get the original compliance again
//	(*mat)=(*S);
//    if (isv->inTension11) (*mat)[0][0]/=isv->E11dam;
//    if (isv->inTension22) (*mat)[1][1]/=isv->E22dam;
//    if (isv->inTension33) (*mat)[2][2]/=isv->E33dam;
//    (*mat)[3][3]/=isv->G12dam;
//    (*mat)[4][4]/=isv->G23dam;
//    (*mat)[5][5]/=isv->G13dam;
//
//#define CompChangeCheck(a) \
//		if(isv->inTension ## a ## a) { \
//			if(stress[a-1] <=-1e-5) { \
//				(*mat)[a-1][a-1] *= isv->E ## a  ## a ## dam; \
//				isv->inTension ## a ## a = false; \
//			} \
//		} else { \
//			if(stress[a-1]>=1e-5) {  \
//				(*mat)[a-1][a-1] /= isv->E ## a ## a ## dam;  \
//				isv->inTension ## a ## a = true; \
//			} \
//		} 
//
//	 CompChangeCheck(1)
//	 CompChangeCheck(2)
//	 CompChangeCheck(3)
//
//	 mat->invert();
//    if(xtNumAngles == 0 && theta == 0){
//        (*localizedCmat) = (*mat);
//    } else {
//        //The first rotation will be about the xtAngles.  This is typically a rotation about the z axis
//        if(xtNumAngles>0)
//        {
//            for(int i=0;i<xtNumAngles;++i){
//                rotateVoigtMaterial( (*localizedCmat),(*mat),xtAngles[i], xtAxis[i]);
//                (*mat) = (*localizedCmat);
//            }
//        }
//        //Now rotate due to mangles
//        if(theta!=0 && dir != 0) {
//            rotateVoigtMaterial( (*localizedCmat),(*mat),theta, dir );
//        }
//    }
//	delete mat;
//
//
//}
////==========================================================================
//void dmElasticMaterial::updateLocalizedCmat2( Rotation* OrientationAtPoint, BasicISV* aISV, double *strain)
//{
//	// when fiber failed, turn off the CM for fiber direction
//    DamageISV* isv = (DamageISV*) aISV;
//
// //Note that this is set up for one rotation only
//
//// Change this to use a saved transformation matrix ...
//	//same for stresses
//
//	int    ii,jj;
//
//	Matrix *mat;
//    mat  = new Matrix(numStrains,numStrains);
//	(*mat)=(*S);
//
//    if (isv->inTension11) (*mat)[0][0]/=isv->E11dam;
//    if (isv->inTension22) (*mat)[1][1]/=isv->E22dam;
//    if (isv->inTension33) (*mat)[2][2]/=isv->E33dam;
//    (*mat)[3][3]/=isv->G12dam;
//    (*mat)[4][4]/=isv->G23dam;
//    (*mat)[5][5]/=isv->G13dam;
//
//	 mat->invert();
//
//	 double stress[6], sum, newstrain[6];
//	 bool moduliChange=false;
//
// 	 //xtang 050102
//	// xtang 050102
//	double *dd=new double [numStrains];
//	for (int i=0; i<numStrains; i++) {dd[i]=strain[i];}
//    //Strains should be rotated through mangles, then xtAngles
//	rotateEngineeringVoigtStrain(newstrain,dd,-theta,dir); ///xtang
//	if (xtNumAngles > 0) // Caution:  assuming from G to L	
//	{
//		rotateEngineeringVoigtStrain(dd,strain,-xtAngles[0], xtAxis[0]);
//	}
//	delete [] dd;
//
//	for (ii=0; ii<3; ii++) {
//		 for(sum=0, jj=0; jj<3; jj++) sum+=(*mat)[ii][jj]*newstrain[jj]; // strain dosen't include eigen strain!!!
//		 stress[ii]=sum;
//	}
//
//	 // get the original compliance again
//	(*mat)=(*S);
//
//	// if fiber is failed or in tension, apply reduction factor to its compliance
//
//    if ((outOfPlaneDir==1) || isv->inTension11) (*mat)[0][0]/=isv->E11dam;
//    if ((outOfPlaneDir==2) || isv->inTension22) (*mat)[1][1]/=isv->E22dam;
//    if ((outOfPlaneDir==3) || isv->inTension33) (*mat)[2][2]/=isv->E33dam;
//
//    (*mat)[3][3]/=isv->G12dam;
//    (*mat)[4][4]/=isv->G23dam;
//    (*mat)[5][5]/=isv->G13dam;
//
//#define CompChangeCheck(a) \
//		if(isv->inTension ## a ## a) { \
//			if(stress[a-1] <=-1e-5) { \
//				(*mat)[a-1][a-1] *= isv->E ## a  ## a ## dam; \
//				isv->inTension ## a ## a = false; \
//			} \
//		} else { \
//			if(stress[a-1]>=1e-5) {  \
//				(*mat)[a-1][a-1] /= isv->E ## a ## a ## dam;  \
//				isv->inTension ## a ## a = true; \
//			} \
//		} 
//
//	 if (outOfPlaneDir !=1 )  {CompChangeCheck(1);}
//	 if (outOfPlaneDir !=2 )  {CompChangeCheck(2);}
//	 if (outOfPlaneDir !=3 )  {CompChangeCheck(3);}
//
//	 mat->invert();
//    if(xtNumAngles == 0 && theta == 0){
//        (*localizedCmat) = (*mat);
//    } else {
//        //The first rotation will be about the xtAngles.  This is typically a rotation about the z axis
//        if(xtNumAngles>0)
//        {
//            for(int i=0;i<xtNumAngles;++i){
//                rotateVoigtMaterial( (*localizedCmat),(*mat),xtAngles[i], xtAxis[i]);
//                (*mat) = (*localizedCmat);
//            }
//        }
//        //Now rotate due to mangles
//        if(theta!=0 && dir != 0) {
//            rotateVoigtMaterial( (*localizedCmat),(*mat),theta, dir );
//        }
//    }
//	delete mat;
//}
//==========================================================================
void dmElasticMaterial::allocateISVs(std::vector<BasicISV*> &isv, const int &numISVs)const
{
 if(DegradationModel==0)
		EXIT_BETA("Degradation Model not set. Cannot create ISV unless Degradation Model is set");
	DegradationModel->createISVSpace(isv, numISVs);
}
//==========================================================================
void dmElasticMaterial::deallocateISVs(std::vector<BasicISV*> &isv) const
{
    DegradationModel->deallocateISVs(isv);
}
//==========================================================================
double dmElasticMaterial::checkStrengthRatio(BasicISV* &isv, double *stress)
{
    DamageISV* Disv = (DamageISV*) isv;
	return DegradationModel->CalculateMaxFailureIndex(Disv, stress, strength);
} 
//==========================================================================
void dmElasticMaterial::updateFailureIndex(BasicISV* &isv, const double *stress) const
{
    DamageISV* Disv = (DamageISV*) isv;
    DegradationModel->CalculateFailureModesAndFailureIndices(Disv, stress, strength);
} 
//==========================================================================
/*
double dmElasticMaterial::checkStrengthRatioWithConcentration(ISV &aISV, double *stress)
{
	return aISV.checkMaxRatio(strength, stress, aISV.concValue);
} 
*/
//==========================================================================
/*
void dmElasticMaterial::printDamageState(ISV* isv, int numISV, ostream *ost)
{
	for (int j=0; j<numISV; j++)
	{
		ISV *isvcp=&isv[j];
		char info;
		info=DegradationModel->getDamageFactorSymbol(isvcp->damFactorForE11);(*ost)<<info;
		info=DegradationModel->getDamageFactorSymbol(isvcp->damFactorForE22);(*ost)<<info;
		info=DegradationModel->getDamageFactorSymbol(isvcp->damFactorForE33);(*ost)<<info;
		info=DegradationModel->getDamageFactorSymbol(isvcp->damFactorForG12);(*ost)<<info;
		info=DegradationModel->getDamageFactorSymbol(isvcp->damFactorForG23);(*ost)<<info;
		info=DegradationModel->getDamageFactorSymbol(isvcp->damFactorForG13);(*ost)<<info;
		info=DegradationModel->getDamageFactorSymbol(isvcp->damFactorForNu12);(*ost)<<info;
		info=DegradationModel->getDamageFactorSymbol(isvcp->damFactorForNu23);(*ost)<<info;
		info=DegradationModel->getDamageFactorSymbol(isvcp->damFactorForNu13);(*ost)<<info;
		(*ost)<<endl;
	}
	(*ost)<<endl;
}
*/
//==========================================================================
void dmElasticMaterial::printDamageState(std::vector<BasicISV*> &isv, int numISV, ostream *ost)
{
	string state;
	for (int j=0; j<numISV; j++)
	{
        DamageISV* Disv = (DamageISV*) isv[j];
		DegradationModel->getDamageStateString(Disv, state);
		(*ost)<< state << endl;
	}
	//(*ost)<<endl;
}
//==========================================================================
void dmElasticMaterial::LoadDamageState(std::vector<BasicISV*> &isv, int numISV, istream *is)
{
	string state;
	for (int j=0; j<numISV; j++) 
	{
        DamageISV* Disv = (DamageISV*) isv[j];
		*is>>state;
		DegradationModel->RestoreISV(Disv, state);
	}
}