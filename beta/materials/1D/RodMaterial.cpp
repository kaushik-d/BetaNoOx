#include "stdafx.h"

#include "RodMaterial.hpp"
#include "math/matrix.hpp"

//=================================================================
//                 RodMaterial
//=================================================================
RodMaterial::RodMaterial()
{
//	BETA_OUT<<"Created rod material object....default constructor"<<endl;
}
//===================================================
RodMaterial::RodMaterial(const RodMaterial& arg):Material(arg)
{
	for(int i=0;i<10;i++) propertyList[i]=arg.propertyList[i];
}
//===================================================
RodMaterial::~RodMaterial()
{
}
//===================================================

bool RodMaterial::read(istream * inStream)
{
	double ea,modulus,height,width, maxTension, maxCompression;
 
   (*inStream)>> modulus >> height >> width >>maxTension>>maxCompression  ;
   ea = modulus * height * width;
   propertyList[0] = modulus;
   propertyList[1] = ea;
   //propertyList[2] = ei;
   propertyList[3] = height;
   propertyList[4] = width;
   propertyList[5] = maxTension;
   propertyList[6] = maxCompression;
   RodMaterial::print(BETA_OUT);

   return true;
}
//===================================================
void RodMaterial::getMaterialState( double strain, double &stress,
												  double &tangentModulus)
{
 stress = propertyList[0] * strain;
 tangentModulus = propertyList[0] ;
}
//===================================================

void RodMaterial::print(ostream & outStream)
{
 outStream << " Rod Material: "<<getGroupName()<<endl;
 outStream << " Group Number:   "<<getGroupNum() <<endl;
 outStream<<"Modulus        = " << getModulus()       <<endl;
 outStream<<"Height         = " << getHeight()        <<endl;
 outStream<<"Width          = " << getWidth()         <<endl;
 outStream<<"maxTension     = " << getMaxTension()    <<endl;
 outStream <<"maxCompression = " << getMaxCompression()<<endl;
}

//=========================================================

void RodMaterial::calculateInitialStress(double *stress,
                                           double *stateVariable)
{
 stress[0] = 0.;
}	
//=========================================================	
void RodMaterial::calculateStress(double *stress, double *strain,
									double *stateVariable)
{
 stress[0] = strain[0] * getEA();   //this is ea * eps= force 
 
 BETA_OUT<<"Stress = "<<stress[0]<<endl;
}		
//=========================================================	
