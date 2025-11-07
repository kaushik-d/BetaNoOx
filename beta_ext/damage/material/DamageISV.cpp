#include "stdafx.h"

#include<cmath>
#include "utility/utility.h"
#include "DamageISV.hpp"
//==========================================================================
DamageISV::DamageISV(int numFMs):
numFailureModes(numFMs),
TriggeredFailureModes(numFMs,false)
{
    damFactorForE11=damFactorForE22=damFactorForE33=1.0f;
	damFactorForG12=damFactorForG23=damFactorForG13=1.0f;
	damFactorForNu12=damFactorForNu23=damFactorForNu13=1.0f;

	inTension11 = inTension22 = inTension33 = 1;
	
	propertyFailedE11=propertyFailedE22=propertyFailedE33=0;
	propertyFailedG12=propertyFailedG23=propertyFailedG13=0;
	propertyFailedNu12=propertyFailedNu23=propertyFailedNu13=0;
	
	for (int i=0; i<6; i++) {
		stress[i]=strain[i]=0.0;
	}
	E11dam=E22dam=E33dam=G23dam=G13dam=G12dam=1.0f;


    ratio = new double [numFailureModes];
    lastRatio = new double [numFailureModes];
    for(int i=0;i<numFailureModes;++i){
        ratio[i]=lastRatio[i]=0.0;
    }
}

//==========================================================================
    
DamageISV::~DamageISV()
{
    delete [] ratio;
    delete [] lastRatio;
}

//==========================================================================

void DamageISV::SwapRatioAndLastRatio(){
    double *temp;
    temp = lastRatio;
    lastRatio = ratio;
    ratio = temp;
}

//==========================================================================

void DamageISV::getDeltaToNextFailure(const double &trialDelta, double& currentMinDeltaLoadFactor, const double &failureRatio) const
{
    double slope = 0.0;
    for(int i=0;i<numFailureModes;++i){
		if(!TriggeredFailureModes[i]){//skip if the element has already failed under this mode
		    slope = (ratio[i] - lastRatio[i])/trialDelta;
			//Check for increasing load resulting in increasing normal or increasing positive shear stress
			if(slope>0.0)
            {
				double deltaToFailure=( 1.0 - lastRatio[i])/slope;
				if(deltaToFailure>0.0 && deltaToFailure<currentMinDeltaLoadFactor)
					currentMinDeltaLoadFactor=deltaToFailure;
			}
            //Check for increasing load resulting in increasing magnitude of negative shear stress
			if(slope<0.0 && (i>=3 && i<=5))
            {
				double deltaToFailure=( -1.0 - lastRatio[i])/slope;
				if(deltaToFailure>0.0 && deltaToFailure<currentMinDeltaLoadFactor)
					currentMinDeltaLoadFactor=deltaToFailure;
			}
		}//End if ith failure mode is not triggered
	}//End Loop over Stress Components
}

//==========================================================================

/*
void DamageISV::UpdateDamageDueToOxidation(double conc)
{
	double maxIncrease=20.8;
	
	if(conc<0) conc=0;
	if(conc>1) conc=1;

	double increase=conc*maxIncrease;

	damFactorForE11conc = damFactorForE22conc = damFactorForE33conc = 1 + increase/100;

}
*/
/*
1/D_2.7027026678716470522450246169839 = D_0.37000000476837158
8/3 = D_2.6666666666666666666666666666667 //I had used D_2.6666667f
3/8 = D_0.375
*/
//From C++ 1/D_0.37f = D_2.7027027606964111
//From PowerCalc: 1/D_0.37f = D_2.7027027027027027027027027027027


//==========================================================================

void DamageISV::DegradeStrength(double * originalSig, double *Sig, double conc)
{
	double factor = 1 - (conc * 0.20); //20% degradation if fully concentrated

	int i;
	//JV101607
	for (i=0; i<9; i++) {
		Sig[i] = originalSig[i] * factor;
	}
}

//==========================================================================

double DamageISV::checkMaxRatio(double * originalSig, double *sigma, double conc)
{
	int i;
	
	double Sig[9];
	DegradeStrength(originalSig, Sig, conc);

	for (i=0; i<9; i++) {
        if(i>=3 && i<=5) {ratio[i]=fabs(sigma[i]/Sig[i]);}
        else {ratio[i]=sigma[i]/Sig[i];}		
	}

	double maxRatio=0.0;

	for(i=0; i<numFailureModes; i++){
		if( (TriggeredFailureModes[i]==false)  && (ratio[i] > maxRatio)) 
		{
			maxRatio = ratio[i];
		} 
	}

	return (maxRatio);

}