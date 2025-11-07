#include "stdafx.h"

#include "utility/utility.h"
#include "plastic.hpp"
#include "materials/3D/elastic.hpp"


#include "math/matrix.hpp"
#include "math/tensor.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

extern void tokenize( char *string,char **tokenList, int &number);

#define LINE "--------------------------------------------------"
#define FORMAT setw(15)<<setprecision(4)


void rotateVoigtMaterial(Matrix &Cr, Matrix &cMat,const double &theta, 
	const int idir);


PlasticMaterial::PlasticMaterial()
{
	cMat_ep=0;
}
//====================================================
void PlasticMaterial::initialize()
{
	ElasticMaterial::initialize();
    cMat_ep = new Matrix(numStrains,numStrains);
	numHardParam=1; // default isotropic hardening

}
//====================================================
PlasticMaterial::~PlasticMaterial(void)          
{
	if(cMat_ep) delete cMat_ep;
}

//====================================================
bool PlasticMaterial::read(istream * inStream)
{
	int i;
	char command[80], *token;
	const char comment[]= "//" ;	
	int  numberOfTokens;
	char *tokenList[20];	// XTANG 980329
	int foundMatch;

	ElasticMaterial::read(inStream);

	while(2==2) 
	{ 
		foundMatch=0;				
		top: (*inStream).getline(command, 80,'\n');		
		tokenize( command, tokenList, numberOfTokens);
		if(numberOfTokens==0)goto top;
		token = tokenList[0];		
		if( strncmp(token,comment,2)==0 ) goto top;		
		banner(token, BETA_OUT);	

		if(COMPARE(token,"readOrthotropicParameters")==0 )
		{
			(*inStream) >> F >> G >> H >> L >> M >> N;
	    
			BETA_OUT<<"Orthotroic Parameters"<<endl;
			BETA_OUT<<"F="<< F <<"\tG="<<G<<"\tH="<<H<<"\tL="<< L <<"\tM="<<M<<"\tN="<<N<<endl;
		} 
      

		if(COMPARE(token,"readMasterCurve")==0)
   		{

   	        (*inStream)>> i >>m_E>>m_sigma>>m_n>>m_sigma0;
			if (m_sigma0 == 0.0) {m_sigma0=1.0;} // avoid zero initial sigma
			if (i==1) masterCurve=POWER;
			if (i==2) masterCurve=RICHARDB;
			if (i==3) masterCurve=USER;

   			BETA_OUT<<masterCurve<<"\t"<<m_E<<"\t"<<m_sigma<<"\t"<<m_n<<"\t"<<m_sigma0<<endl;
   	    }

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

		if(COMPARE(token,"exitPlasticMaterial")==0)
		{
	   		BETA_OUT<<"End of input for this material"<<LINE<<"\n\n"<<endl;
	   		return true; 
		}
	} //end of while        
}

GaussISV * PlasticMaterial::createISVSpace(int aType, int aNum)
{
	//int i;
	GaussISV *tmp;

	ISVType=aType;

	tmp= new GaussISV [aNum]; 
	/*
	for (i=0; i<aNum; i++) {
		tmp[i]= GaussISV();  // 
	}
	*/
    return tmp;
}

void PlasticMaterial::deleteISVSpace(GaussISV * isv)
{
	//int i;
	int aNum=27; //hardwired

	/*
	for (i=0; i<aNum; i++) {
		delete [] isv[i];
	}
	*/
    delete [] isv;
}
/*
double PlasticMaterial::HillFunction()
{
	double s[6]; 
	double f;
	int i;

	for (i=0; i<6; i++) {
		if (numHardParam!=1) {
			s[i]=((*gISV).stress[i]+(*gISV).dStress[i])-((*gISV).state[i]+(*gISV).dState[i])-(*gISV).yStress[i];
		} else {
			s[i]=((*gISV).stress[i]+(*gISV).dStress[i])-(*gISV).yStress[i];
		}
	}

	f= 0.5*(F*(s[1]-s[2])*(s[1]-s[2])+
       G*(s[2]-s[0])*(s[2]-s[0])+
	   H*(s[0]-s[1])*(s[0]-s[1])+
	   2*L*s[4]*s[4] +
	   2*M*s[5]*s[5] +
	   2*N*s[3]*s[3]);

	return (f);
}
*/

// 092501
double PlasticMaterial::HillFunction()
{
	double s[6]; 
	double f;
	int i;

	for (i=0; i<6; i++) {
		if (numHardParam!=1) { //092101
			s[i]=((*gISV).stress[i]-(*gISV).originState[i]+(*gISV).dStress[i])-((*gISV).state[i]+(*gISV).dState[i])-(*gISV).yStress[i];
		} else {
			s[i]=((*gISV).stress[i]-(*gISV).originState[i]+(*gISV).dStress[i])-(*gISV).yStress[i];
		}
	}

	f= 0.5*(F*(s[1]-s[2])*(s[1]-s[2])+
       G*(s[2]-s[0])*(s[2]-s[0])+
	   H*(s[0]-s[1])*(s[0]-s[1])+
	   2*L*s[4]*s[4] +
	   2*M*s[5]*s[5] +
	   2*N*s[3]*s[3]);

	return (f);
}


double PlasticMaterial::m_getYieldStress()
{
	double y=m_Stress();
	if (numHardParam==1) {
		if (y > m_sigma0) {return y;}
		else return m_sigma0;  // yield stress based on the total effective plastic strain 
	} else {
		return m_sigma0;// initial yield stress 
	}

}

double PlasticMaterial::getEffectiveStress()
{
	return sqrt(3.0*HillFunction());
}

const double TOL=1e-6;

double PlasticMaterial::yieldFunction(GaussISV * aISV)
{
	double f, sy;

	// make a local copy of aISV
	gISV=aISV;

	// effectiv plastic strain CHECK this??????:  Isotropic hardening only for now
	sy=m_getYieldStress(); 

	f= HillFunction()/(sy*sy/3.0) - 1.0;

	return (f);
}

void PlasticMaterial::calculateLocalStress(double *stress, double *strain)
{
	int i,j;
	double *C_i;

	// stress, strain and cMat are all in the Material Coordinate System!!!
	for(i=0;i<numStrains;i++) {
		stress[i]=0;
		C_i = (*cMat)[i];
		for(j=0;j<numStrains;j++) {
			stress[i] += C_i[j] * strain[j];    
		}
	}
}	


void PlasticMaterial::calculateLocalStrain(double *strain, double *stress)
{
	int i,j;
	double *S_i;

	// stress, strain and cMat are all in the Material Coordinate System!!!
	for(i=0;i<numStrains;i++) {
		strain[i]=0;
		S_i = (*S)[i];
		for(j=0;j<numStrains;j++) {
			strain[i] += S_i[j] * stress[j];    
		}
	}
}	

/*

// version 1.0  120100

int PlasticMaterial::getStressIncrement(GaussISV * aISV)
{
	gISV=aISV;
	double dse[6], 
		   *s0=(*gISV).stress, 
		   *a0=(*gISV).state, 
		   *ds=(*gISV).dStress, 

		   *dstrain=(*gISV).dStrain,         // 052801
		   *dplasticstrain=(*gISV).dPlasticStrain,  // 052801
		   delasticstrain[6],

		   *da=(*gISV).dState, 
		   h[6],  // hardening functions  ?? variable length (iso=1, kinematic=6)
		   sp[6], // plastic stress: c*dg/ds
		   se[6], // c*df/ds
		   nf[6], // normal of yield function f
		   np[6]; // normal of plastic potential g
	
	double cdump[6][6];

	int iter=0;

	double hard, dl, denom, sy; 
	int i,j;

	// The General Baskward-Euler (return mapping) algorithm is used here.
	//    Reference Book:  
	//		Crisfield, M.A. Non-linear finite Element Analysis of Solids and Structures, Vol. 1
	//      John Wiley & Sons, 1991, pp. 152-
	//
	// Explicit algorithms such as Runge-Kutta methods can also be used
	// here.  Thus, several options will be implemented here for optimal 
	// operation.

	// set up initials before Return Mapping Algorithms

	calculateLocalStress(dse, (*gISV).dStrain);  // elastic predictor
	
	for (i=0;i<6;i++) {	ds[i]=dse[i];}
	for (i=0;i<numHardParam;i++) {da[i]=0;}  // numHardParam=1 (iso), =6 (Kinematic)

	(*gISV).dEffectivePlasticStrain=0.0;
	
	s0=(*gISV).stress;
	a0=(*gISV).state;
	ds=(*gISV).dStress;
	(*cMat_ep) = (*cMat); // use elastic properties

	if (yieldFunction(aISV) > 0.0) {  // entering plastic corrector
		(*gISV).inPlasticRegion=true;
		do { // Return Mapping Algorithm
			hard=hardeningModulus();
			dYield(nf);
			dPlastic(np);
			calculateLocalStress(sp,np);
			denom=yieldFunction(aISV);
			sy=m_getYieldStress();
			dl=denom*sy*sy/3/(hard+dotProduct(nf,sp,6));
			for (i=0;i<6;i++) {ds[i]-=dl*sp[i];}     // stress increment
			hardeningFunction(h);
			for (i=0;i<numHardParam;i++) {da[i]+=dl*h[i];}  // state increment
			(*gISV).dEffectivePlasticStrain+=2.0*dl*getEffectiveStress()/3.0; // eq. 3.13
			denom=yieldFunction(aISV);
			iter++;
		} while ((denom>TOL) && (iter<5000));

		// Calculate the tangent stiffness matrix
		dYield(nf);
		dPlastic(np);
		calculateLocalStress(sp,np);
		calculateLocalStress(se,nf);
		denom=hardeningModulus()+dotProduct(nf,sp,6);

		// update C_ep, elastoplastic constitutive matrix
		for(i=0; i<6; i++) {for(j=0; j<6; j++) {cdump[i][j]=sp[i]*se[j];}}
		for (i=0; i<6; i++) {
			for (j=0; j<6; j++) {
				(*cMat_ep)[i][j]=(*cMat)[i][j]-cdump[i][j]/denom;
			}
		}

		// 052801: update  plastic strain
		calculateLocalStrain(delasticstrain,ds);
	    for (i=0;i<6;i++) {	
			dplasticstrain[i]=dstrain[i]-delasticstrain[i];
		}

	} else {
		if ((*gISV).inPlasticRegion) {(*gISV).inPlasticRegion=false;}
	}
    (*gISV).effectiveStress=getEffectiveStress();
	return (iter);
}



*/

//THIS IS THE WORKING ONE DG_2006_03_25
// versio 2.0   052501

int PlasticMaterial::getStressIncrement(GaussISV * aISV)
{
	gISV=aISV;
	double dse[6], 
		   *s0=(*gISV).stress, 
		   *a0=(*gISV).state, 
		   *ds=(*gISV).dStress, 
		   

  		   *dstrain=(*gISV).dStrain,         // 052801
		   *dplasticstrain=(*gISV).dPlasticStrain,  // 052801
		   delasticstrain[6],

		   *da=(*gISV).dState, 

           h[6],  // hardening functions  ?? variable length (iso=1, kinematic=6)
		   sp[6], // plastic stress: c*dg/ds
		   se[6], // c*df/ds
		   nf[6], // normal of yield function f
		   np[6]; // normal of plastic potential g
	
	double cdump[6][6];

	int iter=0;

	double hard, dl, denom, sy; 
	int i,j;

	// The General Baskward-Euler (return mapping) algorithm is used here.
	//    Reference Book:  
	//		Crisfield, M.A. Non-linear finite Element Analysis of Solids and Structures, Vol. 1
	//      John Wiley & Sons, 1991, pp. 152-
	//
	// Explicit algorithms such as Runge-Kutta methods can also be used
	// here.  Thus, several options will be implemented here for optimal 
	// operation.

	// set up initials before Return Mapping Algorithms

	calculateLocalStress(dse, (*gISV).dStrain);  // elastic predictor
	
	for (i=0;i<6;i++) {	ds[i]=dse[i];}
	for (i=0;i<numHardParam;i++) {da[i]=0;}  // numHardParam=1 (iso), =6 (Kinematic)

	(*gISV).dEffectivePlasticStrain=0.0;
	
	s0=(*gISV).stress;
	a0=(*gISV).state;
	ds=(*gISV).dStress;
	(*cMat_ep) = (*cMat); // use elastic properties

	if (yieldFunction(aISV) > 0.0) {  // entering plastic corrector
		(*gISV).inPlasticRegion=true;
		do { // Return Mapping Algorithm
			hard=hardeningModulus();
			dYield(nf);
			dPlastic(np);
			calculateLocalStress(sp,np);
			denom=yieldFunction(aISV);
			sy=m_getYieldStress();
			dl=dotProduct(nf,dse,6)/(hard+dotProduct(nf,sp,6));
			for (i=0;i<6;i++) {
				ds[i]=dse[i]-dl*sp[i];
			}     // stress increment
			hardeningFunction(h);
			for (i=0;i<numHardParam;i++) {da[i]=dl*h[i];}  // state increment
			(*gISV).dEffectivePlasticStrain=2.0*dl*getEffectiveStress()/3.0; // eq. 3.13
			denom=yieldFunction(aISV);
			iter++;
		} while ((denom>TOL) && (iter<5000));

		// Calculate the tangent stiffness matrix
		dYield(nf);
		dPlastic(np);
		calculateLocalStress(sp,np);
		calculateLocalStress(se,nf);
		denom=hardeningModulus()+dotProduct(nf,sp,6);

		// update C_ep, elastoplastic constitutive matrix
		for(i=0; i<6; i++) {for(j=0; j<6; j++) {cdump[i][j]=sp[i]*se[j];}}
		for (i=0; i<6; i++) {
			for (j=0; j<6; j++) {
				(*cMat_ep)[i][j]=(*cMat)[i][j]-cdump[i][j]/denom;
			}
		}

  		// 052801: update  plastic strain
		calculateLocalStrain(delasticstrain,ds);
	    for (i=0;i<6;i++) {	
			dplasticstrain[i]=dstrain[i]-delasticstrain[i];
		}
	  
	} else {
		if ((*gISV).inPlasticRegion) {(*gISV).inPlasticRegion=false;}
	}
    (*gISV).effectiveStress=getEffectiveStress();
	return (iter);
}


/*
// versio 3.0   092501

// including loading-unloading-reloading

int PlasticMaterial::getStressIncrement(GaussISV * aISV)
{
	gISV=aISV;
	double dse[6], 
		   *s0=(*gISV).stress, 
		   *a0=(*gISV).state, 
		   *ds=(*gISV).dStress, 
		   

  		   *dstrain=(*gISV).dStrain,         // 052801
		   *dplasticstrain=(*gISV).dPlasticStrain,  // 052801
		   delasticstrain[6],

		   *da=(*gISV).dState, 

           h[6],  // hardening functions  ?? variable length (iso=1, kinematic=6)
		   sp[6], // plastic stress: c*dg/ds
		   se[6], // c*df/ds
		   nf[6], // normal of yield function f
		   np[6]; // normal of plastic potential g
	
	double cdump[6][6];

	int iter=0;

	double hard, dl, denom, sy; 
	int i,j;

	// The General Baskward-Euler (return mapping) algorithm is used here.
	//    Reference Book:  
	//		Crisfield, M.A. Non-linear finite Element Analysis of Solids and Structures, Vol. 1
	//      John Wiley & Sons, 1991, pp. 152-
	//
	// Explicit algorithms such as Runge-Kutta methods can also be used
	// here.  Thus, several options will be implemented here for optimal 
	// operation.

	// set up initials before Return Mapping Algorithms

	calculateLocalStress(dse, (*gISV).dStrain);  // elastic predictor
	
	for (i=0;i<6;i++) {	ds[i]=dse[i];}

	for (i=0;i<numHardParam;i++) {da[i]=0;}  // numHardParam=1 (iso), =6 (Kinematic),  

	(*gISV).dEffectivePlasticStrain=0.0;
	
	s0=(*gISV).stress;
	a0=(*gISV).state;
	ds=(*gISV).dStress;
	(*cMat_ep) = (*cMat); // use elastic properties


	if (yieldFunction(aISV) > 0.0) {  // entering plastic corrector
		(*gISV).inPlasticRegion=true;
		do { // Return Mapping Algorithm
			hard=hardeningModulus();
			dYield(nf);
			dPlastic(np);
			calculateLocalStress(sp,np);
			denom=yieldFunction(aISV);
			sy=m_getYieldStress();
			dl=dotProduct(nf,dse,6)/(hard+dotProduct(nf,sp,6));
			for (i=0;i<6;i++) {
				ds[i]=dse[i]-dl*sp[i];
			}     // stress increment
			hardeningFunction(h);
			for (i=0;i<numHardParam;i++) {da[i]=dl*h[i];}  // state increment
			(*gISV).dEffectivePlasticStrain=2.0*dl*getEffectiveStress()/3.0; // eq. 3.13
			denom=yieldFunction(aISV);
			iter++;
		} while ((denom>TOL) && (iter<5000));

		// Calculate the tangent stiffness matrix
		dYield(nf);
		dPlastic(np);
		calculateLocalStress(sp,np);
		calculateLocalStress(se,nf);
		denom=hardeningModulus()+dotProduct(nf,sp,6);

		// update C_ep, elastoplastic constitutive matrix
		for(i=0; i<6; i++) {for(j=0; j<6; j++) {cdump[i][j]=sp[i]*se[j];}}
		for (i=0; i<6; i++) {
			for (j=0; j<6; j++) {
				(*cMat_ep)[i][j]=(*cMat)[i][j]-cdump[i][j]/denom;
			}
		}

  		// 052801: update  plastic strain
		calculateLocalStrain(delasticstrain,ds);
	    for (i=0;i<6;i++) {	
			dplasticstrain[i]=dstrain[i]-delasticstrain[i];
		}
	} else {
		if ((*gISV).inPlasticRegion) { // means entering unloading stage
			for (i=0;i<6;i++) (*gISV).originState[i]=(*gISV).stress[i];
			if ((*gISV).weakMemory) { 
				for (i=0;i<6;i++) (*gISV).state[i]=0;
				(*gISV).effectivePlasticStrain=0.0;
			}
		}
		(*gISV).inPlasticRegion=false;
	}
    (*gISV).effectiveStress=getEffectiveStress();
	return (iter);
}
*/



Matrix * PlasticMaterial::getTangentStiffness(GaussISV * aISV)
{
	gISV=aISV;
	double sp[6], // plastic stress: c*dg/ds
		   se[6], // c*df/ds
		   nf[6], // normal of yield function f
		   np[6]; // normal of plastic potential g
	double cdump[6][6], denom;
	int i,j;
	(*cMat_ep) = (*cMat); // use elastic properties

	
	if ((*gISV).inPlasticRegion) {
		// Calculate the tangent stiffness matrix
		dYield(nf);
		dPlastic(np);
		calculateLocalStress(sp,np);
		calculateLocalStress(se,nf);
		denom=hardeningModulus()+dotProduct(nf,sp,6);
		for(i=0; i<6; i++) {for(j=0; j<6; j++) {cdump[i][j]=sp[i]*se[j];}}
		for (i=0; i<6; i++) { for (j=0; j<6; j++) {
				(*cMat_ep)[i][j]=(*cMat)[i][j]-cdump[i][j]/denom;
		}
		}
	}

	return cMat_ep;
}

double PlasticMaterial::m_dsdep()
{
	// There three function form for master curve:
	// POWER, RICHARDB, USER
	//
	// dsdep: is only depends on the total effective Plastic Strain of current 
	// EffectivePlasticStrain + dEffectivePlasticStrain
	// 

	double ep=(*gISV).effectivePlasticStrain+(*gISV).dEffectivePlasticStrain;


	switch (masterCurve) {

	case POWER: {
		if (fabs(ep) <1e-12) {return m_E;}
		else return m_E*pow(ep,m_n-1.0); // initial yield stress
		}
	case RICHARDB: // eq. 3.28
		{
			return m_E/pow(1+ pow(m_E*ep/m_sigma,m_n),(1.0+m_n)/m_n);
		}
	case USER: {cout << "The material has not been defined yet..."<<endl;}
	}
	return 0;



	return 0;
}




double PlasticMaterial::m_Stress()
{
	double ep=(*gISV).effectivePlasticStrain+(*gISV).dEffectivePlasticStrain;

	switch (masterCurve) {
		case POWER: return m_E*pow(ep,m_n)+m_sigma0; // initial yield stress
		case RICHARDB: // eq. 3.26
			{
				double a=m_E*ep/pow(1+ pow(m_E*ep/m_sigma,m_n),1.0/m_n)+m_sigma0;
				return a;
			}
		case USER: {cout << "The material has not been defined yet..."<<endl;}
	}
	return 0;
}


double PlasticMaterial::m_pStrain()
{
	double ep=(*gISV).effectivePlasticStrain+(*gISV).dEffectivePlasticStrain;
	double es=getEffectiveStress();

	switch (masterCurve) {

	case POWER: {
		if (es > m_sigma0) {return pow((es-m_sigma)/m_E, 1.0/m_n); } 
		else return 0.0;
				}
	case RICHARDB: // eq. 3.26
		{
		if (es > m_sigma0) {	return pow( pow(m_E/es, m_n) - pow(m_E/m_sigma,m_n), -1.0/m_n);  }
		else return 0.0;
		}
	case USER: {cout << "The material has not been defined yet..."<<endl;}
	}
	return 0;
}


double PlasticMaterial::hardeningModulus()
{
	double mm=getEffectiveStress();

	return 4.0/9.0*mm*mm*m_dsdep();
}

void PlasticMaterial::hardeningFunction(double * h)
{
	
	double nf[6];
	double ds[6];
	int i;

	double mm=getEffectiveStress();


	if (numHardParam==1) { // isotropic hardening
		h[0]=2.0/3.0*mm;
	} else {  // kinematic hardening
		for (i=0; i<6; i++)	{
			ds[i]=((*gISV).stress[i]-(*gISV).originState[i]  // 092101
				+(*gISV).dStress[i])-((*gISV).state[i]+(*gISV).dState[i])-(*gISV).yStress[i];
		}
		dYield(nf);
		double a=4.0/9.0*mm*mm*m_dsdep()/dotProduct(nf,ds,6);
		for (i=0; i<6; i++)	{h[i]=ds[i]*a;}
	}
	
}

void PlasticMaterial::dYield(double * nf)
{
	double s[6]; 
	int i;

	for (i=0; i<6; i++) {
		if (numHardParam!=1) {
			s[i]=((*gISV).stress[i]-(*gISV).originState[i]+(*gISV).dStress[i])-((*gISV).state[i]+(*gISV).dState[i])-(*gISV).yStress[i];
		} else {
			s[i]=((*gISV).stress[i]-(*gISV).originState[i]+(*gISV).dStress[i])-(*gISV).yStress[i];
		}
	}

	nf[0]=-G*(s[2] - s[0]) + H*(s[0] - s[1]);
  	nf[1]=F*(s[1] - s[2]) - H*(s[0] - s[1]);
	nf[2]=-F*(s[1] - s[2]) + G*(s[2] - s[0]);
	nf[3]=2*N*s[3];
	nf[4]=2*L*s[4];
	nf[5]=2*M*s[5];
}

void PlasticMaterial::dPlastic(double * np)
{

	// Using associate flow rule 
	dYield(np);
}


//==========================================================================
GaussISV::GaussISV():
materialOrientation()
{
	for (int i=0; i<6; i++) {
		stress[i]=
		strain[i]=
		plasticStrain[i]=  // 052801
		
		state[i]=
		originState[i]=  // 092101 for unloading-reloading

		dStress[i]=
		dStrain[i]=
        dPlasticStrain[i]=  // 052801
		dState[i]=  // back stress
		yStress[i]=0;
	}
	dEffectivePlasticStrain=effectiveStress=effectivePlasticStrain=0;
	inPlasticRegion=false;
	weakMemory=true;
}


void GaussISV::updateISV()
{
	for (int i=0; i<6; i++) {
		stress[i] += dStress[i];
		strain[i] += dStrain[i];
        plasticStrain[i]+=dPlasticStrain[i];  // 052801
		state[i]  += dState[i];  // depends on numHardParam

		dStress[i] = 0.0;
		dStrain[i] = 0.0;
        dPlasticStrain[i] = 0.0;  // 052801
		dState[i] = 0.0;
	}

	effectivePlasticStrain += dEffectivePlasticStrain;
    dEffectivePlasticStrain=0.0;
}

int PlasticMaterial::updateISV(GaussISV * aISV, double * de)  // de: total increment solution
{
	gISV=aISV;

	int i;
	for (i=0;i<6;i++) {(*gISV).dStrain[i]=de[i];}
	return (getStressIncrement(aISV)); 
}


int PlasticMaterial::updateISV_2(GaussISV * aISV, double * de)  // de: total increment solution
{
	gISV=aISV;

	int i;
	for (i=0;i<6;i++) {(*gISV).dStress[i]=de[i];}
	return (getStrainIncrement(aISV)); 
}


void PlasticMaterial::updateLocalizedCmat( Rotation* OrientationAtPoint, GaussISV &aISV)
{
	// output cMat_ep !!!!!!
   *localizedCmat = *getTangentStiffness(&aISV); 

   OrientationAtPoint->rotate_Cmat_from_local_to_global(*localizedCmat,*globalCmat);

//   if(theta==0 || theta ==NULL || dir == 0){
//	  (*localizedCmat) = (*cMat_ep);
//	 // return;
//   }
//   else if(theta!=0){
//	  rotateVoigtMaterial( (*localizedCmat),(*cMat_ep),theta, dir );
//   }
//	//xtang 050102
//    // L to G
//
//	if(xtNumAngles>0)
//	{
//		Matrix *mat;
//		mat  = new Matrix(numStrains,numStrains);
//		(*mat)=(*localizedCmat);
////		(*mat)=(*cMat);   // 03122003 for debug only
//		rotateVoigtMaterial( (*localizedCmat),(*mat),xtAngles[0], xtAxis[0]);
//		
//// (*mat).print("before transformation",6,6,&BETA_OUT);  // 03122003 for debug only
//// (*localizedCmat).print("after transformation",6,6,&BETA_OUT);  // 03122003 for debug only
//
//		delete mat;
//	}
//
//
}
/*


#ifdef DEBUG_PlasticMaterial

#include <fstream>
#include <stdio.h>

///////////////////////////////
//   CONSTITUTIVE TESTER     //
//
//   input: stress or strain
//   output: strain or stress and the change of internal state variables.
//
///////////////////////////////

const int MaxStressNum=6;
void main(void)
{
	PlasticMaterial * pm;
	GaussISV *isv;
	double dstrain[MaxStressNum], dstress[MaxStressNum];
	int i,j;
	

	// read material properties;
	int group;
	char label[300];
	ifstream is("\\2work\\ALPHA_CORE_reconstruct_v2\\sample\\lam\\metal3.mat");
  	is.getline( label, 80,'\n'); 
	is >>group;
	is.getline( label, 80,'\n'); 
	pm = new PlasticMaterial(group, label );
	pm->read(&is);

	// specify the load history
	for (i=0;i<MaxStressNum;i++) {
		dstrain[i]=0.0;
		dstress[i]=0.0;
	}

	isv=new GaussISV();
	//dstrain[2]=0.00022;  // uniaxial strain loading
	//dstrain[4]=0.0006;  // uniaxial strain loading

	dstress[2]=2*60e6/300; 
//	dstress[3]=60e6/3000; 

	double rate=1.0/100;
	double ccmat[6][6];
	int i1,i2;
	Matrix *cc, *cc2;
		
	int iter;
	ofstream os("\\2work\\ALPHA_CORE_reconstruct_v2\\sample\\lam\\matTest.dat");
	for (i=0; i<300; i++) { 
          // stress controlled
		   iter=pm->updateISV_2(isv,dstress);

          // strain controlled
		  // iter=pm->updateISV(isv,dstrain);

	      isv->updateISV();

		  os <<iter<<"\t";
		  os <<isv->strain[0]<<"\t"<<isv->stress[0]<<"\t"<<isv->plasticStrain[0]<<"\t";
		  os <<isv->strain[1]<<"\t"<<isv->stress[1]<<"\t"<<isv->plasticStrain[1]<<"\t";
		  os <<isv->strain[2]<<"\t"<<isv->stress[2]<<"\t"<<isv->plasticStrain[2]<<"\t";
		  os <<isv->strain[3]<<"\t"<<isv->stress[3]<<"\t"<<isv->plasticStrain[3]<<"\t";
		  os <<isv->strain[4]<<"\t"<<isv->stress[4]<<"\t"<<isv->plasticStrain[4]<<"\t";
		  os <<isv->strain[5]<<"\t"<<isv->stress[5]<<"\t"<<isv->plasticStrain[5]<<"\t";

		  os <<isv->effectivePlasticStrain<<"\t"<<isv->effectiveStress<<"\t";


		  os<<endl;

    }
}
#endif
*/






/*
// version 1.0

int PlasticMaterial::getStrainIncrement(GaussISV * aISV)
{
	gISV=aISV;
	double dse[6], 
		   *s0=(*gISV).stress, 
		   *a0=(*gISV).state, 
		   *ds=(*gISV).dStrain, 
		   *da=(*gISV).dState, 
		   h[6],  // hardening functions  ?? variable length (iso=1, kinematic=6)
		   sp[6], // plastic stress: c*dg/ds
		   se[6], // c*df/ds
		   nf[6], // normal of yield function f
		   np[6]; // normal of plastic potential g
	
	double cdump[6][6];

	int iter=0;

	double hard, dl, denom, sy; 
	int i,j;

    for (i=0;i<numHardParam;i++) {da[i]=0;}  // numHardParam=1 (iso), =6 (Kinematic)

	(*gISV).dEffectivePlasticStrain=0.0;
	
	calculateLocalStrain(ds,(*gISV).dStress);

	(*cMat_ep) = (*cMat); // use elastic properties

	if (yieldFunction(aISV) > 0.0) {  // entering plastic corrector
		(*gISV).inPlasticRegion=true;
		do { // Return Mapping Algorithm
			hard=hardeningModulus();
			dPlastic(np);
			denom=yieldFunction(aISV);
			sy=m_getYieldStress();
			dl=denom*sy*sy/3/hard;
			for (i=0;i<6;i++) {ds[i]+=dl*np[i];}     // plastic strain increment
			hardeningFunction(h);
			for (i=0;i<numHardParam;i++) {da[i]+=dl*h[i];}  // state increment
			(*gISV).dEffectivePlasticStrain+=2.0*dl*getEffectiveStress()/3.0; // eq. 3.13
			denom=yieldFunction(aISV);
			iter++;
		} while ((denom>TOL) && (iter<5000));

		
		// update C_ep, elastoplastic constitutive matrix
		dYield(nf);
		dPlastic(np);
		calculateLocalStress(sp,np);
		calculateLocalStress(se,nf);
		denom=hardeningModulus()+dotProduct(nf,sp,6);

		for(i=0; i<6; i++) {for(j=0; j<6; j++) {cdump[i][j]=sp[i]*se[j];}}
		for (i=0; i<6; i++) {
			for (j=0; j<6; j++) {
				(*cMat_ep)[i][j]=(*cMat)[i][j]-cdump[i][j]/denom;
			}
		}
	} else {
		if ((*gISV).inPlasticRegion) {(*gISV).inPlasticRegion=false;}
	}
    (*gISV).effectiveStress=getEffectiveStress();
	return (iter);
}

*/



// version 2.0

int PlasticMaterial::getStrainIncrement(GaussISV * aISV)
{
	gISV=aISV;
	double *s0=(*gISV).stress, 
		   *a0=(*gISV).state, 
		   *ds=(*gISV).dStrain, 

		   *dstrain=(*gISV).dStrain,         // 052801
		   *dplasticstrain=(*gISV).dPlasticStrain,  // 052801
		   delasticstrain[6],

		   *da=(*gISV).dState, 
		   h[6],  // hardening functions  ?? variable length (iso=1, kinematic=6)
		   sp[6], // plastic stress: c*dg/ds
		   se[6], // c*df/ds
		   nf[6], // normal of yield function f
		   np[6]; // normal of plastic potential g
	
	double cdump[6][6];

	int iter=0;

	double hard, dl, denom, sy; 
	int i,j;

    for (i=0;i<numHardParam;i++) {da[i]=0;}  // numHardParam=1 (iso), =6 (Kinematic)

	(*gISV).dEffectivePlasticStrain=0.0;
	
	calculateLocalStrain(ds,(*gISV).dStress);

	(*cMat_ep) = (*cMat); // use elastic properties

	if (yieldFunction(aISV) > 0.0) {  // entering plastic corrector
		(*gISV).inPlasticRegion=true;
		do { // Return Mapping Algorithm
			hard=hardeningModulus();
			dPlastic(np);
			denom=yieldFunction(aISV);
			sy=m_getYieldStress();
			dl=denom*sy*sy/3/hard;
			for (i=0;i<6;i++) {ds[i]+=dl*np[i];}     // plastic strain increment
			hardeningFunction(h);
			for (i=0;i<numHardParam;i++) {da[i]+=dl*h[i];}  // state increment
			(*gISV).dEffectivePlasticStrain+=2.0*dl*getEffectiveStress()/3.0; // eq. 3.13
			denom=yieldFunction(aISV);
			iter++;
		} while ((denom>TOL) && (iter<5000));

		
		// update C_ep, elastoplastic constitutive matrix
		dYield(nf);
		dPlastic(np);
		calculateLocalStress(sp,np);
		calculateLocalStress(se,nf);
		denom=hardeningModulus()+dotProduct(nf,sp,6);

		for(i=0; i<6; i++) {for(j=0; j<6; j++) {cdump[i][j]=sp[i]*se[j];}}
		for (i=0; i<6; i++) {
			for (j=0; j<6; j++) {
				(*cMat_ep)[i][j]=(*cMat)[i][j]-cdump[i][j]/denom;
			}
		}

		// 052801: update  plastic strain
		calculateLocalStrain(delasticstrain,(*gISV).dStress);
	    for (i=0;i<6;i++) {	
			dplasticstrain[i]=dstrain[i]-delasticstrain[i];
		}
	
	} else {
		// 092101		
		if ((*gISV).inPlasticRegion) { // means entering unloading stage
			for (i=0;i<6;i++) (*gISV).originState[i]=(*gISV).stress[i];
			if ((*gISV).weakMemory) { 
				for (i=0;i<6;i++) (*gISV).state[i]=0;
				(*gISV).effectivePlasticStrain=0.0;
			}
		}
		(*gISV).inPlasticRegion=false;
	}
    (*gISV).effectiveStress=getEffectiveStress();
	return (iter);
}

