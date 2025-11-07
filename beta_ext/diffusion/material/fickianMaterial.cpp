#include "stdafx.h"

#include "utility/utility.h"
#include "fickianMaterial.hpp"
#include "math/matrix.hpp"
#include "math/tensor.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

#include "elements/InterpAndSolDerivatives.hpp"

#define LINE "--------------------------------------------------"
#define FORMAT setw(15)<<setprecision(4)

//==============================================================
FickianMaterial::FickianMaterial()
{
	saturationValue=1.0;
}
//==============================================================
FickianMaterial::FickianMaterial(const FickianMaterial& arg):HeatTransferMaterial(arg)
{
saturationValue = arg.saturationValue;
}

//==============================================================
void FickianMaterial::initialize()
{
	setType("FickianMaterial");
	numberOfFluxes    = 3;
	CommonHeatTransferMaterial();
	saturationValue = 1.;
}
//==============================================================
bool FickianMaterial::read(istream * inStream)  // Orthotropic material
{
int i,j,exitFlag;
char *token;	
int  numberOfTokens;
char *tokenList[20];	


for(i=0; i<3; i++){
    for(j=0; j<3; j++) {
        (*cMat)[i][j] =0; 
    }
}

    while(2==2){
        exitFlag=getLineAndTokenize(inStream,
                                    "exitFickianMaterial",
                                    tokenList,
                                    numberOfTokens);
    
        if(exitFlag==1)	{	//Normalize the diffusivity
            for(i=0; i<3; i++) {
                for(j=0; j<3; j++) {
                    (*cMat)[i][j] =(*cMat)[i][j]*(saturationValue);
                }
            }

            (*globalCmat) = (*cMat);

            cMat->print("Normalized cMat",&BETA_OUT);
            return(true);			
        }

        token = tokenList[0];

        if(COMPARE(token,"readDiffusivity")==0 ){
            // Currently limited to orthotropic material
            (*inStream)>>(*cMat)(0,0)>>(*cMat)(1,1)>>(*cMat)(2,2);  // not normalized
            BETA_OUT << "-------------------\tFickianMaterial: " <<
            getGroupName() << "\nGroup Number: " << getGroupNum() << '\n';

            //long oldIOS;
            ios_base:: fmtflags oldIOS;
            int oldPrecision;
            oldIOS = BETA_OUT.setf(ios::scientific);
            oldPrecision = BETA_OUT.precision(4);
            BETA_OUT <<"Diffusivity matrix"<<endl;
            BETA_OUT << "\t" << (*cMat)(0,0)<< "\t" << (*cMat)(0,1)<< "\t" << (*cMat)(0,2)<<endl;
            BETA_OUT << "\t" << (*cMat)(1,0)<< "\t" << (*cMat)(1,1)<< "\t" << (*cMat)(1,2)<<endl;  
            BETA_OUT << "\t" << (*cMat)(2,0)<< "\t" << (*cMat)(2,1)<< "\t" << (*cMat)(2,2)<<endl;  
            BETA_OUT << "-------------------" << endl;
            BETA_OUT.setf(oldIOS);
            BETA_OUT.precision(oldPrecision);
        } //end of readDiffusivity

        if(COMPARE(token,"readSaturationValue")==0 ){
            //default=1
            (*inStream)>>saturationValue;
            BETA_OUT<<"Saturation value = "<< saturationValue<<endl;
        } //end of readSaturationValue



    }//end of while
    return(false);                                    
}


//==========================================================================

void FickianMaterial::calculateFlux(double *flux,
													  InterpAndSolDerivatives *id_ptr)
{
	InterpAndSolDerivatives &id=*id_ptr;
int i,j;
double *C_i, gradientOfT[3];

gradientOfT[0] = id.p_u1_x1;
gradientOfT[1] = id.p_u1_x2;
gradientOfT[2] = id.p_u1_x3;


	for(i=0;i<numberOfFluxes;i++) {
		flux[i]=0;
		 C_i = (*globalCmat)[i];  // normalized!!!
		for(j=0;j<numberOfFluxes;j++) {
			//flux[i] += C_i[j] * gradientOfT[j]/saturationValue;
			flux[i] += C_i[j] * gradientOfT[j];	
		}
	}

flux[0] = -flux[0]; flux[1] = -flux[1]; flux[2] = -flux[2]; 

/* Since the global problem is solved in terms of normalized concentrations, 
 the actual concentration must be defined on an elemen by element basis
...hence it will be appended to the flux output ...  as the fourth column
 To be consistent with the other code, we will calculate the quadraure point value of
the concentration. Later it will be extrapolated to the nodes of the element
*/


/*
Warning: ;

  Code not finished!!!  Need to modify the gradient calculation to "un-normalize" the gradient!

  Easy to implement if I have material update the vol. integral info rather than element
  */
double sum;
double *S= id.S;
double *nodalValues_u1= id.nodalValues_u1;
sum = 0.;
for(i=0; i<id.numberOfInterp; i++) {sum += S[i] * nodalValues_u1[i];}
flux[3] = sum*saturationValue;

}		
//==========================================================================


void FickianMaterial::print(ostream &ostrm)
{
	//long oldIOS;
	ios_base::fmtflags oldIOS;
	int oldPrecision;
		
	BETA_OUT << "-------------------\tFickianMaterial: " <<
			getGroupName() << "\nGroup Number: " << getGroupNum() << '\n';

	oldIOS = ostrm.setf(ios::scientific);
	oldPrecision = ostrm.precision(4);
	ostrm <<"Diffusivity matrix"<<endl;
    ostrm << "\t" << (*cMat)(0,0)<< "\t" << (*cMat)(0,1)<< "\t" << (*cMat)(0,2)<<endl;
	ostrm << "\t" << (*cMat)(1,0)<< "\t" << (*cMat)(1,1)<< "\t" << (*cMat)(1,2)<<endl;  
	ostrm << "\t" << (*cMat)(2,0)<< "\t" << (*cMat)(2,1)<< "\t" << (*cMat)(2,2)<<endl;  
	ostrm << "-------------------" << endl;
	ostrm.setf(oldIOS);
	ostrm.precision(oldPrecision);
}


//==========================================================================
//JV042109 use parent class' function
/*
void FickianMaterial::rotateMaterialToGlobal(double angle, int dir)
{
	int i,j,k;

	if(angle==0 || angle ==NULL || dir == 0){
		(*globalCmat) = (*cMat);
		return; 
	}
 	else if(angle!=0)
    {
		RotationTensor a(angle,dir);  //  This is likely slow!!!
		double aa[3][3];

		// xtang 05-15-00
		// D_new = T * D_old  * T_transpose
		// Hardwired to cartesion tensor (3x3) !!!!

		for(i=0;i<3;i++) {
			for(j=0;j<3;j++) {
			   aa[i][j]=0.0;
			   for(k=0;k<3;k++) {
	              aa[i][j] += a(i+1,k+1)*(*cMat)[k][j];
			   }
			}
		}

		for(i=0;i<3;i++) {
			for(j=0;j<3;j++) {
			   (*globalCmat)[i][j]=0.0;
			   for(k=0;k<3;k++) {
	              (*globalCmat)[i][j] += aa[i][k]*a(j+1,k+1);
			   }
			}
		}
	}
}
*/

//==========================================================================
#undef FORMAT
