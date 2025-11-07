#include "stdafx.h"
#include "BasicDegradationModel.hpp"
#include "dmElastic.hpp"

// Degradation Factors for the model to pick up:

BasicDegradationModel::~BasicDegradationModel(){
    for(int i=0;i<NumFailureModes;++i){
        delete [] DegradationFactors[i];
    }
    delete [] DegradationFactors;
    delete [] FailureModes;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BasicDegradationModel::InitializeDegradationModel(int FiberDirection, char** localTokenList, int ParamStartIndex, int numberOfTokens)
{
    string InputFileName=localTokenList[ParamStartIndex];
    //Allocates arrays for DegradationFactor storage and reads degradation factors from specified file.
    ifstream *DegFacFile=0;
    DegFacFile=filemanager->OpenInputFileFromSearchPath(InputFileName);

    if(DegFacFile==0) 
        EXIT_BETA("Cannot open the degradation file");

    *DegFacFile>> NumFailureModes >> NumDegradationFactorsPerNode;
    FailureModes = new bool [NumFailureModes];

    DegradationFactors = new double* [NumFailureModes];
    for(int i=0; i<NumFailureModes; ++i){
        DegradationFactors[i] = new double[NumDegradationFactorsPerNode];
    }

    for(int i=0; i<NumFailureModes; ++i){
        for(int j=0;j<NumDegradationFactorsPerNode;++j){
            *DegFacFile >> DegradationFactors[i][j];
        }
    }

    *DegFacFile >> alpha >> beta >> gamma;
    filemanager->CloseInputStream(DegFacFile);

    // If the material is transversely isotropic and has been defined such that the fiber direction is not x1, 
    // then the degradation factors must be modified so that damage is correctly applied.
    if(!(FiberDirection == 0 || FiberDirection == 1)){
        TransformDegradationFactors(FiberDirection);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool BasicDegradationModel::CalculateFailureModesAndFailureIndices(DamageISV* isvToDegrade, const double* QuadPointStress, const double* strength) const
{
    //Based on current stress state and material properties, determine which failure modes have been met.
    //Also update the failure ratios of the ISV object.
    
    bool NewFailure = false;
    for ( int i=0; i<3; i++) 
    {
        //Tensile stress, compare to tensile strength
        isvToDegrade->ratio[i] = (QuadPointStress[i]/strength[i]);
        if(isvToDegrade->ratio[i]>1.0){FailureModes[i] = NewFailure = true;}
        else {FailureModes[i] = false;}

        //Compressive Stress, compare to compressive strength
        isvToDegrade->ratio[i+6] = (QuadPointStress[i]/strength[i+6]);//Compressive strength should be negative
        if(isvToDegrade->ratio[i+6]>1.0){FailureModes[i+6] = NewFailure = true;}
        else {FailureModes[i+6] = false;}
    }

    //Shear modes, no compressive failure
    for ( int i=3; i<6; i++) 
    {
        isvToDegrade->ratio[i] = (QuadPointStress[i]/strength[i]);
        if(fabs(isvToDegrade->ratio[i])>1){FailureModes[i] = NewFailure = true;}
        else {FailureModes[i] = false;}
    }
    return NewFailure;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BasicDegradationModel::CalculateMaxFailureIndex(DamageISV* isvToDegrade, double* QuadPointStress, double* strength)
{
    //Based on current stress state and material properties, determine which failure modes have been met.
    //Also update the failure ratios of the ISV object.
    bool NewFailure=CalculateFailureModesAndFailureIndices(isvToDegrade, QuadPointStress, strength);

    MaxFailureIndex=0.0;
    MaxFailureIndex_index=0;
    for (int i=0; i<NumFailureModes; i++) 
    {
        if(!isvToDegrade->TriggeredFailureModes[i] && fabs(isvToDegrade->ratio[i]) > MaxFailureIndex) {
            MaxFailureIndex_index=i;
            MaxFailureIndex=isvToDegrade->ratio[i];
        }
    }

    return MaxFailureIndex;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BasicDegradationModel::TransformDegradationFactors(int FiberDirection){
// This function modifies the array DegradationFactors so that it can be applied to a
// transversely isotropic material whose material properties have been described such
// that the fibers run in the x2 direction (x3 currently not implemented).  This is
// necessary because if the fiber direction is x2 and the material fails due to sigma22,
// then the properties in the 2 direction should be modified like the 1 direction properties
// would be for fiber failure in a material with fiber direction x1.  This method will be
// depreciated if material properties are always defined such that x1 is the fiber direction.

// Assumes failure modes are ordered in the following manner:
//     0 - sigma11
//     1 - sigma22
//     2 - sigma33
//     3 - sigma12
//     4 - sigma23
//     5 - sigma13
//     6 - sigma11
//     7 - sigma22
//     8 - sigma33
//     ... and so on
// It also assumes that the degradation factors are applied in the following order:
//     0 - E11
//     1 - E22
//     2 - E33
//     3 - G12
//     4 - G23
//     5 - G13
//     6 - Nu12
//     7 - Nu23
//     8 - Nu13


    double** temp = new double* [NumFailureModes];
    for(int i=0; i<NumFailureModes; ++i){
        temp[i] = new double [NumDegradationFactorsPerNode];
    }
    // Copy DegradationFactors to temp
    for(int i=0;i<NumFailureModes;++i){
        for(int j=0;j<NumDegradationFactorsPerNode;++j){
            temp[i][j] = DegradationFactors[i][j];
        }
    }

    if(FiberDirection == 2){
        //First, switch rows
        for(int i=0;i<NumDegradationFactorsPerNode;++i){
            //Switch Rows 0 and 1 (sigma 11 and sigma 22 failure modes)
            DegradationFactors[0][i] = temp[1][i];
            DegradationFactors[1][i] = temp[0][i];
            if(NumFailureModes > 6){
                DegradationFactors[6][i] = temp[7][i];
                DegradationFactors[7][i] = temp[6][i];
            }

            //Switch Rows 4 and 5 (sigma 23 and sigma 13 failure modes)
            DegradationFactors[4][i] = temp[5][i];
            DegradationFactors[5][i] = temp[4][i];
            if(NumFailureModes > 9){
                DegradationFactors[10][i] = temp[11][i];
                DegradationFactors[11][i] = temp[10][i];
            }
        }

        //Re-copy to temp
        for(int i=0;i<NumFailureModes;++i){
            for(int j=0;j<NumDegradationFactorsPerNode;++j){
                temp[i][j] = DegradationFactors[i][j];
            }
        }

        //Now, switch columns
        for(int i=0;i<NumFailureModes; ++i){
            //Modify nu12 (column 6), since the value being degraded is now really nu21 = nu12*E22/E11
            DegradationFactors[i][6] = DegradationFactors[i][6]*DegradationFactors[i][1]/DegradationFactors[i][0];

            //Switch columns 0 and 1 (E11 and E22)
            DegradationFactors[i][0] = temp[i][1];
            DegradationFactors[i][1] = temp[i][0];

            //Switch columns 4 and 5 (G23 and G13)
            DegradationFactors[i][4] = temp[i][5];
            DegradationFactors[i][5] = temp[i][4];

            //Switch columns 7 and 8 (nu23 and nu13)
            DegradationFactors[i][7] = temp[i][8];
            DegradationFactors[i][8] = temp[i][7];
        }
    }else{
        cout << "Fiber direction " << FiberDirection << "not supported." << endl;
        EXIT_BETA("Exiting from BasicDegradationModel::TransformDegradationFactors");
    }

    for(int i=0; i<NumFailureModes; ++i){
        delete [] temp[i];
    }
    delete [] temp;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool BasicDegradationModel::UpdateISVDegradation(DamageISV* &isvToDegrade, double* QuadPointStress, double* strength)
{
    // Modifies degradation of ISV based on current stress state and material strengths 

    bool NewFailure = CalculateFailureModesAndFailureIndices(isvToDegrade, QuadPointStress, strength);
    if(!NewFailure){return false;} // don't bother executing code if there is no new failure.

    bool propertyFailedE11underAnyMode=false;
    bool propertyFailedE33underAnyMode=false;
    bool propertyFailedE22underAnyMode=false;
    bool propertyFailedG12underAnyMode=false;
    bool propertyFailedG23underAnyMode=false;
    bool propertyFailedG13underAnyMode=false;
    bool propertyFailedNu12underAnyMode=false;
    bool propertyFailedNu23underAnyMode=false;
    bool propertyFailedNu13underAnyMode=false;

    for(int i=0;i<NumFailureModes;++i){ // Loop over failure modes
        if(FailureModes[i]) {
            if (!isvToDegrade->TriggeredFailureModes[i]) {
                isvToDegrade->TriggeredFailureModes[i]=true;
                if(isvToDegrade->damFactorForE11<DegradationFactors[i][0]){
                    isvToDegrade->damFactorForE11=DegradationFactors[i][0]; propertyFailedE11underAnyMode=true;}
                if(isvToDegrade->damFactorForE22<DegradationFactors[i][1]){
                    isvToDegrade->damFactorForE22=DegradationFactors[i][1]; propertyFailedE22underAnyMode=true;}
                if(isvToDegrade->damFactorForE33<DegradationFactors[i][2]){
                    isvToDegrade->damFactorForE33=DegradationFactors[i][2]; propertyFailedE33underAnyMode=true;}
                if(isvToDegrade->damFactorForG12<DegradationFactors[i][3]){
                    isvToDegrade->damFactorForG12=DegradationFactors[i][3]; propertyFailedG12underAnyMode=true;}
                if(isvToDegrade->damFactorForG23<DegradationFactors[i][4]){
                    isvToDegrade->damFactorForG23=DegradationFactors[i][4]; propertyFailedG23underAnyMode=true;}
                if(isvToDegrade->damFactorForG13<DegradationFactors[i][5]){
                    isvToDegrade->damFactorForG13=DegradationFactors[i][5]; propertyFailedG13underAnyMode=true;}
                if(isvToDegrade->damFactorForNu12<DegradationFactors[i][6]){
                    isvToDegrade->damFactorForNu12=DegradationFactors[i][6]; propertyFailedNu12underAnyMode=true;}
                if(isvToDegrade->damFactorForNu23<DegradationFactors[i][7]){
                    isvToDegrade->damFactorForNu23=DegradationFactors[i][7]; propertyFailedNu23underAnyMode=true;}
                if(isvToDegrade->damFactorForNu13<DegradationFactors[i][8]){
                    isvToDegrade->damFactorForNu13=DegradationFactors[i][8]; propertyFailedNu13underAnyMode=true;}
            }
        }
    }// i

    if(propertyFailedE11underAnyMode) isvToDegrade->propertyFailedE11 = 1;
    if(propertyFailedE22underAnyMode) isvToDegrade->propertyFailedE22 = 1;
    if(propertyFailedE33underAnyMode) isvToDegrade->propertyFailedE33 = 1;
    if(propertyFailedG12underAnyMode) isvToDegrade->propertyFailedG12 = 1;
    if(propertyFailedG23underAnyMode) isvToDegrade->propertyFailedG23 = 1;
    if(propertyFailedG13underAnyMode) isvToDegrade->propertyFailedG13 = 1;
    if(propertyFailedNu12underAnyMode) isvToDegrade->propertyFailedNu12 = 1;
    if(propertyFailedNu23underAnyMode) isvToDegrade->propertyFailedNu23 = 1;
    if(propertyFailedNu13underAnyMode) isvToDegrade->propertyFailedNu13 = 1;

    //return if any new failure occurred
    return (propertyFailedE11underAnyMode || propertyFailedE22underAnyMode || propertyFailedE33underAnyMode
         || propertyFailedG12underAnyMode || propertyFailedG23underAnyMode || propertyFailedG13underAnyMode
         || propertyFailedNu12underAnyMode || propertyFailedNu23underAnyMode || propertyFailedNu13underAnyMode  );

}
//======================================================================
char BasicDegradationModel::getDamageFlag(unsigned failureModeStatus, double damFactor) 
{
    //If damage factor =1.0 that means damage did not occur, so return 0
    //Lets passed arguments were s11failureModeStatus and E11damFactor
    //otherwise, find out the value of damage factor and 
        //If Damaged by amount alpha, return A (if s11 failed) or a (if s11 did not failed)
        //If Damaged by amount  beta, return B (if s11 failed) or b (if s11 did not failed)
        //If Damaged by amount gamma, return C (if s11 failed) or c (if s11 did not failed)
    
    if (fabs(damFactor - 1.0) < 1e-5) { return '0';}

    if (fabs(damFactor-alpha) < 1e-5) 
    {
        if (failureModeStatus) { return 'A';} else {return 'a';}
    }
    else if (fabs(damFactor-beta) < 1e-5)
    {
        if (failureModeStatus) { return 'B'; } else {return 'b';}
    }
    else if (fabs(damFactor-gamma) < 1e-5) 
    {
        if (failureModeStatus) { return 'C';} else {return 'c';}
    }
    else {
        EXIT_BETA("Problem exists in getDamageFlag()");
    }
}
//======================================================================
char BasicDegradationModel::getDamageFactorSymbol(double damFactor) 
{
    //If damage factor =1.0 that means damage did not occur, so return 0
    //otherwise:
        //If Damaged by amount alpha, return A 
        //If Damaged by amount  beta, return B 
        //If Damaged by amount gamma, return C 
    
    if (fabs(damFactor - 1.0) < 1e-5) { return '0';}
    else if (fabs(damFactor-alpha) < 1e-5) { return 'A';} 
    else if (fabs(damFactor-beta) < 1e-5)  { return 'B'; }
    else if (fabs(damFactor-gamma) < 1e-5) { return 'C';} 
    else {
        EXIT_BETA("Problem exists in getDamageFactorSymbol()");
    }
}
//======================================================================
void BasicDegradationModel::getDamageStateString(DamageISV* isv, string &state)
{
    char info;
    state.clear();
    for(int i=0;i<NumFailureModes;++i){
        if(isv->TriggeredFailureModes[i]) state.push_back('1');else state.push_back('0');
    }

    info=getDamageFactorSymbol(isv->damFactorForE11);    state.push_back(info);
    info=getDamageFactorSymbol(isv->damFactorForE22);    state.push_back(info);
    info=getDamageFactorSymbol(isv->damFactorForE33);    state.push_back(info);
    info=getDamageFactorSymbol(isv->damFactorForG12);    state.push_back(info);
    info=getDamageFactorSymbol(isv->damFactorForG23);    state.push_back(info);
    info=getDamageFactorSymbol(isv->damFactorForG13);    state.push_back(info);
    info=getDamageFactorSymbol(isv->damFactorForNu12);   state.push_back(info);
    info=getDamageFactorSymbol(isv->damFactorForNu23);   state.push_back(info);
    info=getDamageFactorSymbol(isv->damFactorForNu13);   state.push_back(info);
}
//======================================================================
void BasicDegradationModel::RestoreISV(DamageISV* isv, string &state)
{
    for(int i=0;i<NumFailureModes;++i){
        if(state[i]=='1') isv->TriggeredFailureModes[i]=1; else isv->TriggeredFailureModes[i]=0;
    }

    isv->propertyFailedE11 =setDamageFactor(state[NumFailureModes+0], isv->damFactorForE11);
    isv->propertyFailedE22 =setDamageFactor(state[NumFailureModes+1], isv->damFactorForE22);
    isv->propertyFailedE33 =setDamageFactor(state[NumFailureModes+2], isv->damFactorForE33);
    isv->propertyFailedG12 =setDamageFactor(state[NumFailureModes+3], isv->damFactorForG12);
    isv->propertyFailedG23 =setDamageFactor(state[NumFailureModes+4], isv->damFactorForG23);
    isv->propertyFailedG13 =setDamageFactor(state[NumFailureModes+5], isv->damFactorForG13);
    isv->propertyFailedNu12=setDamageFactor(state[NumFailureModes+6], isv->damFactorForNu12);
    isv->propertyFailedNu23=setDamageFactor(state[NumFailureModes+7], isv->damFactorForNu23);
    isv->propertyFailedNu13=setDamageFactor(state[NumFailureModes+8], isv->damFactorForNu13);
}
//======================================================================
unsigned BasicDegradationModel::setDamageFactor(char c, double &v) 
{
    unsigned failFlag;

    if (c=='A' || c=='a') { 
        v=alpha; 
        if(c=='A') { failFlag=1;} else {failFlag=0;}
    }
    else if (c=='B' || c=='b') {
        v=beta; 
        if(c=='B') { failFlag=1;} else {failFlag=0;}
    }        
    else if (c=='C' || c=='c') { 
        v=gamma; 
        if(c=='C') { failFlag=1;} else {failFlag=0;}
    }
    else if (c=='0') {
        v=1.0;
        failFlag=0;
    }
    return failFlag;
}
//======================================================================
unsigned restoreBlackDamageFlag(char c,double &v) 
{
    unsigned f;
    const double alpha=(double)0.125;
    const double beta=(double)0.37;
    const double gamma1=(double)0.01;

    if (c=='A' || c=='a') { 
        v=alpha; 
        if(c=='A') { f=1;} else {f=0;}
    }
    if (c=='B' || c=='b') {
        v=beta; 
        if(c=='B') { f=1;} else {f=0;}
    }        
    if (c=='C' || c=='c') { 
        v=gamma1; 
        if(c=='C') { f=1;} else {f=0;}
    }
    if (c=='0') {
        v=1.0;
        f=0;
    }
    return f;
}
//======================================================================
/*
ISV * BasicDegradationModel::createISVSpace(int aNum)
{
    ISV *tmp=0;
    tmp=new ISV [aNum];
    for(int i=0; i<aNum; i++){
        tmp[i].numFailureModes=NumFailureModes;
        tmp[i].TriggeredFailureModes.reserve(NumFailureModes);
        tmp[i].TriggeredFailureModes.assign(NumFailureModes,false);
    }
    return tmp;
}
*/
void BasicDegradationModel::createISVSpace(std::vector<BasicISV*> &isv, int numISVs)const
{
    //For this case, allocate ISVs out of a continuous block of memory...
    //isv[0] points to the allocated block...
    if(!isv.empty()){
        deallocateISVs(isv);
    }
    isv.resize(numISVs);
    isv[0] = new DamageISV[numISVs];
    DamageISV* start = (DamageISV*)isv[0];
    for(int i=1;i<numISVs;++i){
        isv[i]=&(start[i]);
    }
}
//======================================================================
void BasicDegradationModel::deallocateISVs(std::vector<BasicISV*> &isv) const 
{
    if(isv[0]){
        delete [] (DamageISV*)isv[0];
        isv.assign(isv.size(),NULL);
    }
}
//======================================================================
void BasicDegradationModel::updateSMatrix(Matrix *mat, DamageISV *isv)
{
    DamageISV &aISV=*isv;

    (*mat)[0][0] = (*mat)[0][0] * aISV.damFactorForE11;
    (*mat)[1][1] = (*mat)[1][1] * aISV.damFactorForE22;
    (*mat)[2][2] = (*mat)[2][2] * aISV.damFactorForE33;
    (*mat)[3][3] = (*mat)[3][3] * aISV.damFactorForG12;
    (*mat)[4][4] = (*mat)[4][4] * aISV.damFactorForG23;
    (*mat)[5][5] = (*mat)[5][5] * aISV.damFactorForG13;
    (*mat)[0][1] = (*mat)[0][1] * aISV.damFactorForE11 / aISV.damFactorForNu12;
    (*mat)[1][0] = (*mat)[0][1];
    (*mat)[0][2] = (*mat)[0][2] * aISV.damFactorForE11 / aISV.damFactorForNu13;
    (*mat)[2][0] = (*mat)[0][2];
    (*mat)[1][2] = (*mat)[1][2] * aISV.damFactorForE22 / aISV.damFactorForNu23;
    (*mat)[2][1] = (*mat)[1][2];
}
//======================================================================

BasicDegradationModel::BasicDegradationModel(const BasicDegradationModel &arg)
{

    NumFailureModes = arg.NumFailureModes;
    NumDegradationFactorsPerNode = arg.NumDegradationFactorsPerNode;

    FailureModes = new bool [NumFailureModes];
    for (int j=0;j<NumFailureModes;j++){
        FailureModes[j]= arg.FailureModes[j]; }

    DegradationFactors = new double * [NumFailureModes];
    for (int i=0;i<NumFailureModes;i++) {
        DegradationFactors[i] = new double [NumDegradationFactorsPerNode];

        for(int j=0;j<NumDegradationFactorsPerNode;j++) {
            DegradationFactors[i][j] = arg.DegradationFactors[i][j];}
    }
        
    alpha = arg.alpha;
    beta = arg.beta;
    gamma = arg.gamma;

}
