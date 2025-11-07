#include "stdafx.h"

#include "utility/utility.h"
#include "heatTransferMaterial.hpp"
#include "math/matrix.hpp"
#include "math/tensor.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

#include "elements/InterpAndSolDerivatives.hpp"

#define LINE "--------------------------------------------------"
#define FORMAT setw(15)<<setprecision(4)

//==============================================================
HeatTransferMaterial::HeatTransferMaterial()
{
	cMat=0;
	globalCmat = 0;
	numberOfFluxes = 3;   
	setType("HeatTransferMaterial");
}
//--------------------------------------------------------------------
HeatTransferMaterial::~HeatTransferMaterial(void)
{
	if(cMat) delete cMat;
	if(globalCmat) delete globalCmat;
}
//--------------------------------------------------------------------
HeatTransferMaterial::HeatTransferMaterial(const HeatTransferMaterial& arg)
{
angleZ = arg.angleZ;
numberOfFluxes = arg.numberOfFluxes;
numberOfRotationAngles = arg.numberOfRotationAngles;
currentRotationAngle = arg.currentRotationAngle;

for(int i=0;i<4;i++){
	rotationAxis[i] = arg.rotationAxis[i];
	rotationAngle[i] = arg.rotationAngle[i];
}

cMat = new Matrix(*arg.cMat);
globalCmat = new Matrix(*arg.globalCmat);


}
//--------------------------------------------------------------------
void HeatTransferMaterial::CommonHeatTransferMaterial(void)
{
	cMat       = new Matrix(numberOfFluxes,numberOfFluxes);
	globalCmat = new Matrix(numberOfFluxes,numberOfFluxes);
}
//--------------------------------------------------------------------
void HeatTransferMaterial::initialize()
{
	BETA_OUT<<"I'm in HeatTransferMaterial(int group, char *groupName)"<<endl;
	cMat=0;
	globalCmat = 0;
	numberOfFluxes    = 3;
	setType("HeatTransferMaterial");
	CommonHeatTransferMaterial();
}
//===========================================================
bool HeatTransferMaterial::read(istream * inStream)
{
	int i,j,exitFlag;
	char *token;	
	int  numberOfTokens;
	char *tokenList[20];	


    // Hardwired for orthotropic properties
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            (*cMat)[i][j] =0;
        }
    }

    while(2==2){
        exitFlag=getLineAndTokenize(inStream,
                                    "exitHeatTransferMaterial",
                                    tokenList,
                                    numberOfTokens);
        if(exitFlag==1) return(true);

        token = tokenList[0];

        if(COMPARE(token,"readConductivity")==0 ){
            (*inStream) >> (*cMat)(0,0) >> (*cMat)(1,1) >> (*cMat)(2,2);

            (*globalCmat) = (*cMat);
            cout << "-------------------\tHeatTransferMaterial: " <<
            getGroupName() << "\nGroup Number: " << getGroupNum() << '\n';


            //long oldIOS;
            ios_base:: fmtflags oldIOS;
            int oldPrecision;
            oldIOS = cout.setf(ios::scientific);
            oldPrecision = cout.precision(4);
            cout <<"Conductivity matrix"<<endl;
            cout << "\t" << (*cMat)(0,0)<< "\t" << (*cMat)(0,1)<< "\t" << (*cMat)(0,2)<<endl;
            cout << "\t" << (*cMat)(1,0)<< "\t" << (*cMat)(1,1)<< "\t" << (*cMat)(1,2)<<endl;  
            cout << "\t" << (*cMat)(2,0)<< "\t" << (*cMat)(2,1)<< "\t" << (*cMat)(2,2)<<endl;  
            cout << "-------------------" << endl;
            cout.setf(oldIOS);
            cout.precision(oldPrecision);
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
    }//end of while
    return(false);                                    
}


//==========================================================================

void HeatTransferMaterial::calculateFlux(double *flux,
													  InterpAndSolDerivatives *id_ptr)
{
	InterpAndSolDerivatives &id=*id_ptr;
int i,j;
double *C_i, gradientOfT[3];
gradientOfT[0] = id.p_u1_x1;
gradientOfT[1] = id.p_u1_x2;
gradientOfT[2] = id.p_u1_x3;

//BETA_OUT<<"Flux"<<endl;
//BETA_OUT<<gradientOfT[0]<< "  "<<gradientOfT[1]<<"  "<<gradientOfT[2]<<endl;

//BETA_OUT<<"gradientOfT[0]="<<gradientOfT[0]<<endl;
	for(i=0;i<numberOfFluxes;i++) {
		flux[i]=0;
		 C_i = (*globalCmat)[i];
		for(j=0;j<numberOfFluxes;j++) {
			flux[i] += C_i[j] * gradientOfT[j];			                       
		}
	}

//	(*globalCmat).print("Global cMat",&OUT);
flux[0] = -flux[0]; flux[1] = -flux[1]; flux[2] = -flux[2]; 
//BETA_OUT<<"Fluxes"<<endl;
//BETA_OUT<<flux[0]<<"  "<<flux[1]<<"  "<<flux[2]<<endl;
}		
//==========================================================================


void HeatTransferMaterial::print(ostream &ostrm)
{
	//long oldIOS;
	ios_base:: fmtflags oldIOS;
	int oldPrecision;
		
	cout << "-------------------\tHeatTransferMaterial: " <<
			getGroupName() << "\nGroup Number: " << getGroupNum() << '\n';

	oldIOS = ostrm.setf(ios::scientific);
	oldPrecision = ostrm.precision(4);
	ostrm <<"Conductivity matrix"<<endl;
   ostrm << "\t" << (*cMat)(0,0)<< "\t" << (*cMat)(0,1)<< "\t" << (*cMat)(0,2)<<endl;
	ostrm << "\t" << (*cMat)(1,0)<< "\t" << (*cMat)(1,1)<< "\t" << (*cMat)(1,2)<<endl;  
	ostrm << "\t" << (*cMat)(2,0)<< "\t" << (*cMat)(2,1)<< "\t" << (*cMat)(2,2)<<endl;  
	ostrm << "-------------------" << endl;
	ostrm.setf(oldIOS);
	ostrm.precision(oldPrecision);
}


//==========================================================================
void HeatTransferMaterial::rotateMaterialToGlobal(const Rotation* OrientationAtPoint)
{
    OrientationAtPoint->rotate_3D_2nd_Order_Tensor_from_local_to_global(*cMat,*globalCmat);
}

//==========================================================================
//void HeatTransferMaterial::rotateMaterialToGlobal(double angle, int axis)
//{
//	int i,j,k;
//
//	if(angle==0 || angle ==NULL || axis == 0){
//		(*globalCmat) = (*cMat);
//		return; 
//	}
// 	else if(angle!=0)
//    {
//		RotationTensor a(angle,axis);  //  This is likely slow!!!
//		double aa[3][3];
//
//		// xtang 05-15-00
//		// D_new = T * D_old  * T_transpose
//		// Hardwired to cartesion tensor (3x3) !!!!
//
//		for(i=0;i<3;i++) {
//			for(j=0;j<3;j++) {
//			   aa[i][j]=0.0;
//			   for(k=0;k<3;k++) {
//	              aa[i][j] += a(i+1,k+1)*(*cMat)[k][j];
//			   }
//			}
//		}
//
//		for(i=0;i<3;i++) {
//			for(j=0;j<3;j++) {
//			   (*globalCmat)[i][j]=0.0;
//			   for(k=0;k<3;k++) {
////	              (*globalCmat)[i][j] += aa[i][k]*a(k+1,j+1);//JV042109 bug: this does a T * D_old * T (should be T * D_old * transpose(T) ) 
//	              (*globalCmat)[i][j] += aa[i][k]*a(j+1,k+1); 
//			   }
//			}
//		}
//	}
//}
//==========================================================================
#undef FORMAT
