#include "stdafx.h"
#include "utility/utility.h"

#include <string.h>

#include "RodElas.hpp"

//======================================================================

RodElas::RodElas(void){ }
RodElas::~RodElas(void){ }

//======================================================================
void RodElas::calculateInternalForce()
{
 internalForce = stress * area;
}
//======================================================================
void RodElas::		      calculateK()
{//=====================================================================
 double factor;
	factor = tangentModulus * area/Lo ;
   (*Ke)[0][0] =  factor;
   (*Ke)[0][1] = -factor;
   (*Ke)[1][0] = -factor;
   (*Ke)[1][1] =  factor;
   //BETA_OUT<<"Stiffness coef. in Ke = "<<factor<<endl;
}
//======================================================================
void RodElas::		      calculateStrain() 
{//====================================================================== 
double dudx   = (xDisp[1]-xDisp[0])/Lo;
strain = dudx ;                        // Engineering
}
//======================================================================
void RodElas::		      calculateStress()
{//=====================================================================
	RodMaterial* material = (RodMaterial*)BasicElement::material;
	material->getMaterialState(strain,stress,tangentModulus);
	if(printStressFlag==1){BETA_OUT<<"stress = "<<stress<<endl;}
}
//======================================================================
void RodElas::		      localizeData(double *uGlobal)
{//=====================================================================
 int i; 
   for(i=0;i<numNodesPerElement;i++) {
      xCoor[i] = node[i].x;            
		}
	Lo = fabs(xCoor[1]- xCoor[0]) ;
	Ln = Lo;
   if(uGlobal != NULL) localizeSolution(uGlobal);
}

//======================================================================
void RodElas::		      localizeSolution(double *uGlobal)
{//=====================================================================
 int i;
 int firstDof;
 //Change this to use extractSolution(..) <= BasicElement
   for(i=0 ; i<numNodesPerElement;i++) {
      firstDof = node[i].getFirstDof();
      xDisp[i] = uGlobal[firstDof];
      }


   Ln = Lo + xDisp[1] - xDisp[0];
}
//======================================================================
void RodElas::		      print(ostream &ostrm)
{//=====================================================================
	ostrm << "------------\n";
	ostrm << "ElementType: RodElas\n";
}
//======================================================================

void RodElas::		      update(char * action, double *uGlobal)
{//=====================================================================
 int calculateStressFlag,calculateInitialForceFlag=0,
     calculateKFlag,     calculateForceFlag;
 
 localizeData(uGlobal); 

 calculateKFlag=calculateForceFlag=calculateStressFlag=0;

 if(strstr(action,"I") !=NULL)calculateInitialForceFlag  = 1;
 if(strstr(action,"K") !=NULL)calculateKFlag = 1;
 if(strstr(action,"F") !=NULL)calculateForceFlag = 1;
 if(strstr(action,"S") !=NULL)calculateStressFlag = 1;
 

 cout<<"action = "<<action<<endl;
 double region;

 RodMaterial* material = (RodMaterial*)BasicElement::material;
 material->getMaterialState(strain,stress,tangentModulus);
 area = (material->getHeight()) * (material->getWidth());

 zeroArrays();

 region = Lo * area;    // This is original area !!

 if(calculateForceFlag==1||analysisType==1)
	{
	 calculateStrain();
	 calculateStress();

	 calculateInternalForce();
	 Fe[0] = -internalForce;
	 Fe[1] =  internalForce;
	}

 if(calculateKFlag==1) { calculateK(); }


//Error check
//..add static member to basicfe to keep track of total volume
	if(region <= 0) banner("Zero or negative volume for element.",BETA_OUT);



if(calculateInitialForceFlag == 1)
{ initialFe[0] = 0.; initialFe[1] = 0.; } 
 }  //end of update()

//======================================================================
void RodElas::		      zeroArrays()
{//=====================================================================
 int i,j;
 for(i=0;i<numDof;i++) {
   Fe[i]=0;
   for(j=i;j<numDof;j++)
   (*Ke)(i,j) = 0;
   }
}//end of zeroArrays
 
//======================================================================
//                            Static Data for RodElas
//======================================================================
// Warning:  Avoid static data in any new code!
 double   RodElas::xCoor[4];
 double   RodElas::xDisp[4];
 double   RodElas::Lo;
 double   RodElas::Ln;
 double   RodElas::strain;
 double   RodElas::area;
 double   RodElas::internalForce;
 int      RodElas::printStressFlag;
 double   RodElas::stress;
 double   RodElas::tangentModulus;
//======================================================================
