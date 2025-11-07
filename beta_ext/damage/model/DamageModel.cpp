#include "stdafx.h"

#include <ctime>
#include <string>
#include <sstream>
#include <omp.h>
#include "DamageModel.hpp"

#include "factory/Factory.hpp"
#include "utility/formWriter.hpp"
#include "math/gaussQuadrature/gaussPointList.hpp"
#include "mesh/MeshUtility.h"

#include "elements/3D/ElasticityElement3D.hpp"
#include "damage/element/dmElasticElement.hpp"

#include "../material/BasicDegradationModel.hpp"
//=========================================================================
extern Factory	*factory;
void deleteFiles(string files);
void setVectorFromNodalValueFile(char* filename, double *vector, int nn, NodeGroup &node);
BasicDegradationModel* CreateDegradationModel(string DegradationModelName, FileManager *filemanager);
//-------------------------------------------------------------------
CreateErrorHandler(DamageModel);
//--------------------------------------------------------------------
DamageModel::DamageModel(void)
{
	dispRequired=false;
	stepIncrements=1;
	NodalLCS_stressRequired=false;//-DG_Oct24_2006
	stepIncrementsNodalLCS_stress=1;//-DG_Oct24_2006
	NodalGCS_stressRequired=false ;//-DG_Oct24_2006
	stepIncrementsNodalGCS_stress=1;//-DG_Oct24_2006
	Quad_stressRequired=false ;//-DG_Oct24_2006
	stepIncrementsQuad_stress=1;//-DG_Oct24_2006
	Damage_Factors_Required=false; //-DG_March29_2007
	stepIncrementsDamage_Factors=1; //-DG_March29_2007

    restartVolAvgStrainFile = "";
    restartVolAvgStressFile = "";

	OutputBeforeIteration=false;
    UseNewtonRaphsonIteration=false;
	DamageInducingFactor=1.000;

	LoadFactor=1.000;
	LoadFactorPrevious=1.000;
	InitialLoadFactor=1.000;
    NewtonRaphsonResidualTolerance = 1.0e-4;
    MaxNewtonRaphsonIterationsBeforeReassemble = 10;

	LoadStepping=DamageModel::NoLoadStepping;

	count=0;

	TransientTimeUnit=minutes;
	mainTimeCount=0.0;
}
//--------------------------------------------------------------------
DamageModel::~DamageModel(void)
{
}
//--------------------------------------------------------------------
int DamageModel::processCommand(char** localTokenList,ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;
	char newFileName[80];

	char token[200];
	int  foundMatch=0;


strcpy(token,localTokenList[0]);
BETA_OUT<<"\nCommand processor = 'DamageModel::processCommand'          "<<token<<endl;


	if( COMPARE(token,"openFile")==0 ) {   //This is for this pass only
		//Should be 3 tokens
		strcpy(newFileName,localTokenList[1]);
		currentStream = openFileForInput(newFileName);
		strcpy(token,localTokenList[2]);
		BETA_OUT<<"Input stream for command ***"<<token<<"*** = "<<newFileName<<endl;
		BETA_OUT<<"Stream(openFile) = "<<currentStream<<endl;
		BETA_OUT<<"Stream(openFile) = "<<(*currentStream)<<endl;  
	}

	//----------- Insert new commands below this line ----------------------
	IS_COMMAND(DoDMAnalysis)
	IS_COMMAND(SetDegradationModelType)		

	if( COMPARE(token,"ReadVolS1S2")==0 ){
		(*currentStream)>>volS1>>volS2;
		BETA_OUT<<" volS1 = "<<volS1<<"   volS2 = "<<volS2<<endl;
		OK
	}

	if( COMPARE(token,"useCompressionModification")==0 )
	{
		 DMElasticElement::compressionModification=true;
		 // xtang 05202003
		 DMElasticElement::jamming=false;
		 BETA_OUT<<" useCompressionModification"<<endl;
		 OK
	}

	if( COMPARE(token,"ReadIterationData")==0 ){
		(*currentStream)>>MaxAllowableLoadSteps>>MaxAllowableIterations;
		BETA_OUT<<" MaxAllowableLoadSteps = "<<MaxAllowableLoadSteps<<"    MaxAllowableIterations = "<<MaxAllowableIterations<<endl;
		OK
	}

	if( COMPARE(token,"InitialDamageFile")==0 ){
		initialDamageFile=localTokenList[1];
		if( doesFileExist(initialDamageFile) == false)
			EXIT_BETA("initialDamageFile ("<< initialDamageFile <<") does not exist!");
		BETA_OUT<<" Initial Damage File: "<<initialDamageFile<<endl;
		OK;
	}

    if( COMPARE(token,"RestartVolAvgStrainFile")==0 ){
		restartVolAvgStrainFile=localTokenList[1];
		if( doesFileExist(restartVolAvgStrainFile) == false)
			EXIT_BETA("restartVolAvgStrainFile ("<< restartVolAvgStrainFile <<") does not exist!");
		BETA_OUT<<" Volume Average Strain File for restart: "<<restartVolAvgStrainFile<<endl;
		OK;
	}

    if( COMPARE(token,"RestartVolAvgStressFile")==0 ){
		restartVolAvgStressFile=localTokenList[1];
		if( doesFileExist(restartVolAvgStressFile) == false)
			EXIT_BETA("restartVolAvgStressFile ("<< restartVolAvgStressFile <<") does not exist!");
		BETA_OUT<<" Volume Average Stress File for restart: "<<restartVolAvgStressFile<<endl;
		OK;
	}

	if( COMPARE(token,"SetDamageInducingFactor")==0 ) {
		DamageInducingFactor=atof(localTokenList[1]);
		BETA_OUT	<< "DamageInducingFactor = " << DOUBLE_FORMAT << DamageInducingFactor << endl;
		OK
	}

	if( COMPARE(token,"SetLoadSteppingScheme")==0 ) {
		if(numToken==2){
			if( COMPARE(localTokenList[1],"Adaptive")==0 )
				LoadStepping=DamageModel::Adaptive;
			else if( COMPARE(localTokenList[1],"Uniform")==0 )
				LoadStepping=DamageModel::Uniform;
			else if( COMPARE(localTokenList[1],"NoLoadStepping")==0 )
				LoadStepping=DamageModel::NoLoadStepping;
			else EXIT_BETA("Cannot recognize LoadSteppingScheme: " << localTokenList[1]);
		}else
			EXIT_BETA("Incorrect syntax for SetLoadSteppingScheme");
		OK
	}

    if( COMPARE(token,"UseNewtonRaphsonIteration")==0 ) {
		UseNewtonRaphsonIteration = true;
		BETA_OUT	<< "UseNewtonRaphsonIteration" << endl;
		OK
	}

    if( COMPARE(token,"SetNewtonRaphsonResidualTolerance")==0 ) {
		NewtonRaphsonResidualTolerance = atof(localTokenList[1]);
		BETA_OUT	<< "NewtonRaphsonResidualTolerance set to " << NewtonRaphsonResidualTolerance << endl;
		OK
	}

    if( COMPARE(token,"SetMaxNewtonRaphsonIterationsBeforeReassemble")==0 ) {
		MaxNewtonRaphsonIterationsBeforeReassemble = atoi(localTokenList[1]);
		BETA_OUT	<< "MaxNewtonRaphsonIterationsBeforeReassemble set to " << MaxNewtonRaphsonIterationsBeforeReassemble << endl;
		OK
	}
	//----------- Insert new commands above this line ----------------------

	if(currentStream!=originalStream) {
		BETA_OUT<<"Close alternate input stream"<<endl;
		(*currentStream).close();
		delete currentStream;
		currentStream=originalStream;
		BETA_OUT<<"Automatically reverts to default stream "<<currentStream
			<<" upon return to readCommands"<<endl;
	}

	if(foundMatch ==0){	
        BETA_OUT<<"Process command by parent\n";
		foundMatch=ElasticityModel::processCommand(localTokenList,currentStream, numToken);
	}

	return foundMatch;
}
//======================================================================
void DamageModel::ProcessSetElementProperty(string command, istream * inStream)
{
	int intCommand=0 ;
	enum {
		SET_ELEMENT_MOISTURE, //xtang 122000
		SET_ELEMENT_ISV,      // +xtang 990420 
		SET_ELEMENT_ISV_new,  // +xtang 991031 
		SET_ELEMENT_ISV_gauss,  // +xtang 991031 
		SET_ELEMENT_ISV_propertyDegrdationModel //Dgoyal_Oct_15_2006
	};

	ChangeToUpper(command)
	if(command=="SETELEMENTMOISTURE")
		{intCommand=SET_ELEMENT_MOISTURE;}
	else if(command=="SETELEMENTISV")
		{intCommand=SET_ELEMENT_ISV;}
	else if(command=="SETELEMENTISV_NEW")
		{intCommand=SET_ELEMENT_ISV_new;}
	else if(command=="SETPROPERTYDEGRADATIONMODELTYPE")
		{intCommand=SET_ELEMENT_ISV_propertyDegrdationModel;}
	else {
		return ElasticityModel::ProcessSetElementProperty(command, inStream);
	}

    int i,first,last,increment;
    double moisture; 

    while (2==2) {				
	    (*inStream)>>first>>last>>increment;
 	        //BETA_OUT<<"First, last, increment="
		    //	<<first<<" "<<last<<" "<<increment<<endl;
	    if(first<0) return;

	    //Error checking...
	    if(first < 0 || first    > mesh->numElements || first > last
				     || increment<1 )
	      {
		    BETA_OUT<<"Fatal error...input is bad"<<endl;
	        BETA_OUT<<"First, last, increment="
			    <<first<<" "<<last<<" "<<increment<<endl;
	        exit(1);
        }

	    BasicElement *e=0;

        switch(intCommand) {
		        //xtang 122000	
		    case SET_ELEMENT_MOISTURE:
		        (*inStream)>>moisture;
		        BETA_OUT << "\tMoisture = " << moisture << " for elements "
		         << first <<" to "<<last << " by " << increment << '\n';
		        for(i=first;i<=last;i+=increment) {
			        e=&mesh->element[i];
			        e->setMoisture(moisture);}
		        break;
	        default:
		        BETA_OUT<<"Match not found"<<endl;
		        exit(1);
        }

    }//end of while (2==2)
} //end of ProcessSetElementProperty
//=======================================================================
//fix this function - the output is different from the alpha version FIXXXXXXXX
void DamageModel::outVolumeStrConstituents(double *disp,int matVol, int numCpt, double ** cpt, ofstream * volumeAverageout, bool byConstituent=false)
{
	// this function works for multi-regions and multi-constituents.  In addition to the selected
	// regions, the volume averaging can also be carried out for the individual constituent in the
	// regions.  

	/* 1. Count constituents in the model (by material ID)
	   2. Set up a lookup table for material IDs
	   3. Allocate the data arrays for the volume averages
	   4. Loop over the model and get the volume averages
	   5. Output the averages.
	*/
	
#define maxStresses 6
#define FORMAT setw(15)<<setprecision(6)
#define FF1 setw(7)<<setprecision(4)<<setiosflags(ios::fixed)

	int numMats=0;     // number of material groups in the model
	int matID[100];  // can hold 100 mat groups at most.

    double *volAvgStress, *volAvgStrain, volume, *SED;  
    double **TotalvolAvgStress, **TotalvolAvgStrain, *Totalvolume;
	double *TotalSED;

    int i, is,  jj, mid;

	double x,y,z;
	double xc[20], yc[20],zc[20];
	
	bool inRegion;
	int numElements=mesh->numElements;
	DMElasticElement *e=0;

	for(i=0;i<100;i++) {matID[i]=0;}

	if (byConstituent) {
		for(i=0;i<numElements;i++) {
			e=(DMElasticElement*)&mesh->element[i];

			is=e->getMaterialNumber();
			jj=getIndex(is,matID,numMats+1);
			if (jj < 0) {numMats++; matID[numMats]=is; }
		}
	}

    volAvgStress = new double [maxStresses];
    volAvgStrain = new double [maxStresses];
    SED = new double [3];

	TotalvolAvgStress = new double * [numMats+1];
	TotalvolAvgStrain = new double * [numMats+1];
    Totalvolume =  new double [numMats+1];
	TotalSED = new double [numMats+1];

	for (i=0; i<=numMats; i++) {
		TotalvolAvgStress[i] = new double [maxStresses];
		TotalvolAvgStrain[i] = new double [maxStresses];

		Totalvolume[i]=0.0;
		TotalSED[i]=0.0;
		for(jj=0;jj<maxStresses;jj++){
			TotalvolAvgStress[i][jj]=0.0;
			TotalvolAvgStrain[i][jj]=0.0;
		}
	}

    for(i=0;i<numElements;i++) 
	{
		e=(DMElasticElement*)&mesh->element[i];
	    is=e->getMaterialNumber();
		mid=getIndex(is,matID,numMats+1);

		if ((matVol !=0) && (matVol !=e->getMaterialNumber())) 
		{
			continue;
		}

		inRegion=true;
		if (numCpt>0) 
		{
			inRegion=false;
			x=y=z=0.0;
			e->extractNodalCoordinates(xc,yc,zc);
			for (jj=0; jj<e->numNodesPerElement; jj++) 
			{
				x+=xc[jj];
				y+=yc[jj];
				z+=zc[jj];
			}
			x/=e->numNodesPerElement;
			y/=e->numNodesPerElement;
			z/=e->numNodesPerElement;

			for (jj=0; jj< numCpt; jj++) 
			{
				inRegion= inRegion || !(
					x > cpt[jj][3] || x < cpt[jj][0] || 
					y > cpt[jj][4] || y < cpt[jj][1] || 
					z > cpt[jj][5] || z < cpt[jj][2]);
			}
		}
		
		if (! inRegion) 
		{
			continue;
		}

        e->calculateDofList();
        volume=0.;
        ((ElasticityElement3D*)e)->getVolAvgValues(disp,volAvgStress,volAvgStrain, SED, &volume);  
		
	    Totalvolume[0]   += volume;
		TotalSED[0] += SED[0];

		if (mid >0) 
		{
			Totalvolume[mid] += volume;
			TotalSED[mid] += SED[0];
		}
        
        for(is=0; is<maxStresses; is++)
		{
           TotalvolAvgStress[0][is] += volAvgStress[is];
           TotalvolAvgStrain[0][is] += volAvgStrain[is];
		   if (mid >0) 
		   {
				TotalvolAvgStress[mid][is] += volAvgStress[is];
				TotalvolAvgStrain[mid][is] += volAvgStrain[is];
		   }
        }
    
	}

	ElasticMaterial * matPoint;
	e=(DMElasticElement*)&mesh->element[0]; // Fix this... what should be done ??

	for (i=0; i<=numMats; i++) 
	{
		*(volumeAverageout) <<matID[i]<<"\t"<<Totalvolume[i]<<"\t";  
		matPoint = (ElasticMaterial * )  (&materialList[matID[i]]);
		/*if (matPoint) 
		{
			*(volumeAverageout) << matPoint->xtNumAngles<<"\t"<<matPoint->xtAngles[0]<<"\t"<<matPoint->xtAxis[0];
		}
		else*/ *(volumeAverageout) << -1<<"\t"<<0<<"\t"<<0;
		*(volumeAverageout)<<"\t"<<TotalSED[i];
		*(volumeAverageout)<<endl; 
		for(is=0; is<maxStresses; is++)  
			*(volumeAverageout)<<FORMAT<<TotalvolAvgStrain[i][is]/Totalvolume[i];
	    *(volumeAverageout)<<endl;
		for(is=0; is<maxStresses; is++)  
			*(volumeAverageout)<<FORMAT<<TotalvolAvgStress[i][is]/Totalvolume[i];
	    *(volumeAverageout)<<endl;
	}
    *(volumeAverageout)<<endl;

    delete [] volAvgStress;  // xt062901
    delete [] volAvgStrain;// xt062901
    delete [] SED;// xt062901
    delete [] Totalvolume;// xt062901.
    delete [] TotalSED;// xt062901.

	for (i=0; i<=numMats; i++) 
	{
		delete [] TotalvolAvgStress[i];// xt062901
		delete [] TotalvolAvgStrain[i];// xt062901
	}

    delete [] TotalvolAvgStress;// xt062901
    delete [] TotalvolAvgStrain;// xt062901

#undef FF1
#undef maxStresses
#undef FORMAT
}
//======================================================================
void DamageModel::SetInitialState()
{
	if(initialDamageFile.empty()==false)
		LoadDamageState(initialDamageFile);

}
//======================================================================
void DamageModel::DoDMAnalysis(istream * inStream)
{
//============================================================
// o for linear elastic material only
// o scaling load step                  
// o solving for q_i+1 at each loadstep, not delta-q_i+1
//============================================================
	checkModel();
	createAssembler();
	allocateForAnalysis();
    equations.setUsePreviousFactorizedMatrixAtRuntime(true);
	InitializeOptionalOutput();

	SetInitialState();

	int countStart=0;

	if(doRestart && !restartInitiated){
		LoadRestartData();
		countStart=restart_iterationNumber;
		restartInitiated=true;
        calculateZeroLoadState();
	}else{
        count = 0;
        //Currently the code depends on thermal loads not causing damage.
        //If thermal loads do cause damage, then it will <probably> result
        //in the LoadFactor jumping up to 1.0 immediately since the model is
        //not allowed to take "negative" loadsteps
        calculateZeroLoadState();
	    countStart=1;
	}
	////////////////////////////////////////////////////////////
    double *currentSolutionInEquations;
	currentSolutionInEquations=equations.getSolution();//get pointer to solution in equations object.
	for(count=countStart; count<=MaxAllowableLoadSteps; count++){    // "Load loop" begins here
        //Increase load so failure occurs (or maximum load is reached)
		double NewLoadFactor=calculateNewLoadFactor();
		if(NewLoadFactor > 1.00000001){
			BETA_OUT << "Next Loadstep exceeds specified maximum load. Terminating load stepping." << endl;
			cout	 << "Next Loadstep exceeds specified maximum load. Terminating load stepping." << endl;
			break;
		}
        //Iterate on load
        SolveLoadStep();
    } // end load step loop
    
    CalculateEnergy(internalForce, currentSolution);  // jae 981023

    if(equations.getUsePreviousFillInOrdering() == true){
		equations.ReleaseMemory();
	}

	FinalizeOptionalOutput();
} // end of DoDMAnalysis
//===================================================================================
void DamageModel::SolveLoadStep(const bool &OutputData)
{
    iIter=0;
    if(UseNewtonRaphsonIteration){ cout << "Load step #" << count;}
    else{cout << "Load step #" << count << " Iteration: " << setw(4) << iIter;}
    SolveLinearStep();

    copyVector(equations.getSolution(),currentSolution,totalNumDof);
    if(OutputBeforeIteration && count != 0 && OutputData) WriteOptionalOutput(PRE_DAMAGE);
    //Volume Average Data always needs to be output pre-damage...
    if(!OutputBeforeIteration && count != 0 && Volume_AveragesRequired) printVolumeAverageData();

    //Some things for Newton-Raphson
    bool recalculateStiffnessMatrix = false;

    do {  // iteration loop begins here
        bool hasNewDamage=UpdateISV(); // Check ISV: for new damage 
        //outputVolumeDistribution();//Nobody using the output from this function right now. uncomment if the need arises
        if(hasNewDamage==false && iIter != 0) break; //no damage, therefore no iteration. Residual is used only for getting delta-q, not for convergence.

        if(!UseNewtonRaphsonIteration) {cout << "\b\b\b\b" << setw(4) << iIter;}
        else {cout << endl << "Iteration: " << setw(4) << iIter;}

        if(UseNewtonRaphsonIteration){
            if(recalculateStiffnessMatrix){
                assembler->assembleK(NULL);
                recalculateStiffnessMatrix = false;
            }
            double residualNorm;
            // Update the Residual
            // CalculateResidual() uses int(B^T * sigma * dV) to assemble the resultant vector 
            // (it doesn't use the global K matrix)
            double lastResidualNorm = residualNorm = CalculateResidual(currentSolution);
            double *deltaSolution;
            //Begin Newton-Raphson iterations using old stiffness matrix
            //(if using Pardiso, this should still be factorized so 
            // the solve should only require a back-substitution)
            // This can't improve performance with iterative solvers
            int MaxNewtonRaphsonIterationsBeforeReassemble = 5;
            int NewtonRaphsonIterations = 0;
            while(true){
                cout << " N-R Iteration: " << setw(4) << NewtonRaphsonIterations << " Residual: " << scientific << setw(10) << residualNorm;

                if(residualNorm < NewtonRaphsonResidualTolerance ) break; //Exit Condition for N-R
                ++NewtonRaphsonIterations;

                //Solve for and apply increment to displacement using the tangent stiffness matrix (which may already be factorized, but might be outdated)
                deltaSolution=equations.incrementalSolutionUsingResidual(currentSolution);
                equations.addIncrementalDisplacements(deltaSolution);
                copyVector(equations.getSolution(),currentSolution,totalNumDof);

                //Calculate the residual using the new solution
                lastResidualNorm = residualNorm;
                residualNorm = CalculateResidual(currentSolution);
                    
                //If your residual is getting worse, then your old stiffness matrix is probably very wrong.
                //Assemble the stiffness matrix for the current damage state and use it.
                if(lastResidualNorm < residualNorm ){
                    assembler->assembleK(NULL);
                }

                //If you exceed MaxNewtonRaphsonIterationsBeforeReassemble, then recalculate the global stiffness matrix next time
                if(NewtonRaphsonIterations > MaxNewtonRaphsonIterationsBeforeReassemble) recalculateStiffnessMatrix = true;
                for(int i = 0;i<54;++i) {cout << "\b";} //Erases previous NR Iteration output data on this line
            }
            //for(int i = 0;i<54;++i) {cout << "\b";} //Erases previous NR Iteration Output Data
        }// End solution using Newton-Raphson iteration

        else {// Solve iteration exactly by assembling and running a linear solve.
            assembler->assembleKandF(NULL);
            equations.solve();
            copyVector(equations.getSolution(),currentSolution,totalNumDof);
        }

        iIter++;
    } while (iIter <= MaxAllowableIterations);  //end of iterate loop
    cout << endl;
    if(iIter > MaxAllowableIterations) EXIT_BETA("Exceeded maximum allowed number of iterations per loadstep!");
    if(OutputData) {WriteOptionalOutput(POST_DAMAGE);}
}
//===================================================================================
void DamageModel::calculateZeroLoadState()
{
    equations.zeroLoadVector(); // Set up applied loadVector
	LoadFactor=LoadFactorPrevious=0.0;
	for(int i=0; i<(int)loads.size(); i++){// Send loads to equation object
        if(loads[i] == NULL) break;
        if (loads[i]->isIncremental()) { 
            loads[i]->setRate(LoadFactor);
        }
        loads[i]->applyTo(equations,&mesh->node,mesh->numNodes, this);
	}

    assembler->assembleKandF(NULL);
    SolveLoadStep(!doRestart); // Only output data if this isn't a restart.

    if(!doRestart) {
        WriteOptionalOutput(POST_DAMAGE);
        calculateDeltaLoadForNextFailure();
    }

	InitialLoadFactor=LoadFactor;
}
//===================================================================================
double DamageModel::calculateNewLoadFactor()
{
    LoadFactorPrevious=LoadFactor;
    double NewLoadFactor;

    if (LoadStepping==DamageModel::Adaptive) { //adaptive
        if(count!=1){ //If this is the first loadstep, the load is already properly set.
            updateFailureIndices();
            calculateDeltaLoadForNextFailure();
        }
    } else if(LoadStepping==DamageModel::Uniform){ //uniform
        NewLoadFactor=LoadFactor+(1.0-InitialLoadFactor)/double(MaxAllowableLoadSteps);
        LoadFactor = NewLoadFactor;
    } else if(LoadStepping==DamageModel::NoLoadStepping){ //NoLoadStepping
        NewLoadFactor=1.0;
        LoadFactor = NewLoadFactor;
    }
    return LoadFactor;
}
//===================================================================================
void DamageModel::calculateDeltaLoadForNextFailure(){
    //Note - assumes the model is in an equilibrium state

	//Solve for a trial load
	double TrialDelta = 0.1;
	double OldLoadFactor = LoadFactor;
	LoadFactor = LoadFactor + TrialDelta;
	ApplyLoadFactorToLoads();

    //AssembleKandF must be used currently so that the constraints get properly accounted for in equations::pseudoLoadVector
	assembler->assembleKandF(NULL);
	double *currentSolutionInEquations=equations.solve();
	copyVector(currentSolutionInEquations,currentSolution,totalNumDof); // U_r after solving

	//Now find the smallest delta LoadFactor that will cause failure at a quad point
	updateFailureIndices(); // Goes through elements and updates failure index (does not degrade properties)
	double minDeltaLoadFactor = 1.0-OldLoadFactor; //Largest that this value could be and there still be meaningful loadsteps to perform

    //Loop over Elements to set minDeltaLoadFactor
	for(int i=0;i<mesh->numElements;++i){
		DMElasticElement* e = (DMElasticElement*) &(mesh->element[i]);
        e->getDeltaToNextFailure(TrialDelta,minDeltaLoadFactor);
    }//End Loop over Elements
	//Now set the new load factor
	LoadFactor = (OldLoadFactor+minDeltaLoadFactor)*DamageInducingFactor;
}
//===================================================================================
void DamageModel::outVolumeDistribution(int numCol, const char *cdump, int count, int iIter)
{
    // Limitations: 10 material system, 50 levels
	// region defined as whole FE model.
    // material su=ystems are numbered as 1, 2, 3 ...
	

	const int numMat=40;
	const int levels=50;
	int col=numCol;

	double sH[numMat][12][50], vH[numMat][12][50];  // mat, col, level
	double smin[numMat][12], smax[numMat][12], sstep[numMat][12]; // mat, col
	
	double s[12], x[3], v, totV[numMat];
	int grp;

	int i,j,k, iele, ip;

	// 0. Open space for
	for (i=0; i < numMat; i++) {
		totV[i]=0.0;
		for (j=0; j< col; j++) {
			smin[i][j]=1.0e300; 
			smax[i][j]=-1.0e300;
		    for (k=0; k< levels; k++) {
				sH[i][j][k]=0.0;
				vH[i][j][k]=0.0;
			}
		}
	}

	eQuadInfo->seekg(0,ios::beg);

	// 2. get min and max of column values for each mat'l 
    for(iele=0;iele<mesh->numElements;iele++) {

		
		int totalNumIPs=((DMElasticElement * )(&mesh->element[iele]))->getTotalNumberOfIPs();

		//////// cout << " Test out" << iele << endl;
		for(ip=0;ip<totalNumIPs;ip++) {
			// get data from "_quadInfo.dat"

			eQuadInfo->read((char *)s, col*sizeof(double));
    		eQuadInfo->read((char *) &grp, sizeof(int));
    		eQuadInfo->read((char *) &v, sizeof(double));
    		eQuadInfo->read((char *) x, sizeof(x));

			totV[grp]+=v;
	
			/*
			for (i=0; i< 7; i++) {cout << s[i] <<" "; }
			cout << grp <<" ";
			cout << v << " ";
			cout << x[0] <<" "<< x[1]<<" "<< x[2] <<endl;
			*/

		    for (i=0; i< col; i++) { 
			  if (s[i] <= smin[grp][i] ) smin[grp][i] =s[i];	
			  if (s[i] >= smax[grp][i] ) smax[grp][i] =s[i];
			}
		}
	}
	
	for (i=0; i < numMat; i++) {
		for (j=0; j< col; j++) 	{  
			sstep[i][j] = (smax[i][j]-smin[i][j])/ double (levels);	
		}
	}


	for (i=0; i < numMat; i++) {
		for (j=0; j< col; j++) 	{  
  			for (k=0; k<levels; k++) {
				sH[i][j][k] = smin[i][j]+sstep[i][j]*(double (k) +0.5); 
			}
		}
	}
	
	eQuadInfo->seekg(0,ios::beg);

	// get histogram 
    for(iele=0;iele<mesh->numElements;iele++) {
		
		int totalNumIPs=((DMElasticElement * )(&mesh->element[iele]))->getTotalNumberOfIPs();

		for(ip=0;ip<totalNumIPs;ip++) {
			// get data from "_quadInfo.dat"

    		eQuadInfo->read((char *)s, col*sizeof(double));
    		eQuadInfo->read((char *) &grp, sizeof(int));
    		eQuadInfo->read((char *) &v, sizeof(double));
    		eQuadInfo->read((char *) x, sizeof(x));


			for (i=0; i< col; i++) { 
                 if (fabs(sstep[grp][i]) <= 1e-30) {
				    k=0;
				 }
                 else {
                    k=int ( fabs((s[i]-smin[grp][i])/sstep[grp][i] - 1.0)  ) ;
					if (k < 0 || k > 49) {k=0;}
				 }
                 vH[grp][i][k]=vH[grp][i][k]+v;
			}
			
		}
	}

	// output result for accumulative    
	double aa;

	// format of "VolDist.dat"
	// step, iter, material

	int nn=4;
	for (i=1; i< nn; i++) {  // 3 materials are allowed
		
		if (totV[i]==0.0) {totV[i]=1.0;}
		else {totV[i]=totV[i]/100.0;}  // for normalizing volume

		(*eVolDist).write((char*) &count, sizeof(int));
		(*eVolDist).write((char*) &iIter, sizeof(int));
		(*eVolDist).write((char*) &i, sizeof(int));

		// sprintf(fn,"%d_%s",i,cdump);
		// ost.open(fn);
  	    for (k=0; k<levels; k++) {
		    for ( j=0; j<col; j++) {
				aa=0.0;
			    for(int l=k; l<levels; l++)	{aa=aa+vH[i][j][l];}
				aa=aa/totV[i];
				// ost << sH[i][j][k] << "  "<<aa/totV[i]<<" ";
		        (*eVolDist).write((char*) &sH[i][j][k], sizeof(double));
		        (*eVolDist).write((char*) &aa, sizeof(double));
			}
			// ost << endl;
		}
		// ost.close();
	}

	
}
//===================================================================================
void DamageModel::printVolumeAverageData()
{
    double stressVolumeSum[6], strainVolumeSum[6],
           elementStressVolume[6], elementStrainVolume[6],
           volumeSum=0.0, elementVolume;
    double dummy[3]; //This is for the SED parameter in getVolAvgValues
    for(int i=0;i<6;++i){
        stressVolumeSum[i]=strainVolumeSum[i]=0.0;
    }
    DMElasticElement* e;
    for(int i=0;i<mesh->numElements;++i){
        e = (DMElasticElement*) &mesh->element[i];
        e->getVolAvgValues(currentSolution,elementStressVolume,elementStrainVolume,dummy,&elementVolume);
        for(int j=0;j<6;++j){
            stressVolumeSum[j]+=elementStressVolume[j];
            strainVolumeSum[j]+=elementStrainVolume[j];
        }
        volumeSum+=elementVolume;
    }
    for(int i=0;i<6;++i){
        *volAvgStressFile << stressVolumeSum[i]/volumeSum << "\t";
        *volAvgStrainFile << strainVolumeSum[i]/volumeSum << "\t";
    }
    *volAvgStressFile << endl;
    *volAvgStrainFile << endl;

}
//===================================================================================
int DamageModel::processOptionalOutput(char** localTokenList,ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;

	char token[200];
	int  foundMatch=0;

    strcpy(token,localTokenList[0]);
    BETA_OUT<<"\nCommand processor = 'DamageModel::processOptionalOutput'          "<<token<<endl;

    //========
	if( COMPARE(token,"Displacements")==0 )
	{
		dispRequired=true;
		stepIncrements=1;
		if(numToken==2){
			stepIncrements=atoi(localTokenList[1]);
			if(stepIncrements < 1)
				EXIT_BETA("Invalid Option in OptionalOutput");
		}
		OK 
	}

    if( COMPARE(token,"GlobalForces")==0 )
	{
		Global_ForcesRequired=true;
		stepIncrementsGlobal_Forces=1;
		if(numToken==2){
			stepIncrementsGlobal_Forces=atoi(localTokenList[1]);
			if(stepIncrementsGlobal_Forces < 1)
				EXIT_BETA("Invalid Option in OptionalOutput");
		}
		OK 
	}

	if( COMPARE(token,"LCS_stress")==0 )
	{
		NodalLCS_stressRequired=true;
		stepIncrementsNodalLCS_stress=1;
		if(numToken==2){
			stepIncrementsNodalLCS_stress=atoi(localTokenList[1]);
			if(stepIncrementsNodalLCS_stress < 1)
				EXIT_BETA("Invalid Option in OptionalOutput");
		}
		OK 
	}

	if( COMPARE(token,"GCS_stress")==0 )
	{
		NodalGCS_stressRequired=true;
		stepIncrementsNodalGCS_stress=1;
		if(numToken==2){
			stepIncrementsNodalGCS_stress=atoi(localTokenList[1]);
			if(stepIncrementsNodalGCS_stress < 1)
				EXIT_BETA("Invalid Option in OptionalOutput");
		}
		OK 
	}

	if( COMPARE(token,"Quad_stress")==0 )
	{
		Quad_stressRequired=true;
		stepIncrementsQuad_stress=1;
		if(numToken==2){
			stepIncrementsQuad_stress=atoi(localTokenList[1]);
			if(stepIncrementsQuad_stress < 1)
				EXIT_BETA("Invalid Option in OptionalOutput");
		}
		OK 
	}

	if( COMPARE(token,"Damage_Factors")==0 )
	{
		Damage_Factors_Required=true;
		stepIncrementsDamage_Factors=1;
		if(numToken==2){
			stepIncrementsDamage_Factors=atoi(localTokenList[1]);
			if(stepIncrementsDamage_Factors < 1)
				EXIT_BETA("Invalid Option in OptionalOutput");
		}
		OK 
	}

    if( COMPARE(token,"Binary_Damage_File")==0 )
	{
		Binary_Damage_File_Required=true;
		OK 
	}

	if( COMPARE(token,"OutputBeforeIteration")==0 )
	{
		if(numToken==2){
			if( COMPARE(localTokenList[1],"true")==0 )
				OutputBeforeIteration=true;
			else if( COMPARE(localTokenList[1],"false")==0 )
				OutputBeforeIteration=false;
			else EXIT_BETA("Incorrect syntax for OutputBeforeIteration");
		}else
			EXIT_BETA("Incorrect syntax for OutputBeforeIteration");
		OK
	}


	if(foundMatch ==0){	
		foundMatch=ElasticityModel::processOptionalOutput(localTokenList,currentStream, numToken);
	}

	return foundMatch;
}
//==========================================================
void DamageModel::WriteOptionalOutput(bool BeforeDamage)
{
#define FORMAT setw(18)<<setprecision(12)

	SET_SCIENTIFIC(BETA_OUT);

	BETA_OUT << count <<"\t"<< FORMAT << mainTimeCount <<"\t"<< LoadFactor <<"\t";
    BETA_OUT<< iIter << "\t";
	if (DMElasticElement::jamming)
		BETA_OUT << "JM" << "\t";
	else 
		BETA_OUT << "NJ" << "\t";
    BETA_OUT << endl;

#undef FORMAT 

	ofstream *ost;
	int i;
	int numElements=mesh->numElements;
	DMElasticElement *e=0;

//	bool lastTimestep=false;
//	if (totalTimeSteps==mainTimeStepCount) lastTimestep=true;

//////////////////////////////////////////////////////////////////
	if(dispRequired &&
	   (count==1 || count % stepIncrements == 0) )	{
		stringstream number;
		number << count;
		string filename;
		if(BeforeDamage)
			filename = number.str() + ".pre.disp" ;
		else
			filename = number.str() + ".disp" ;
		ost = filemanager->OpenOutputStream(filename);
		mesh->printNodalData(currentSolution, "displacements", ost );
		filemanager->CloseOutputStream(ost);
	}
    if(Global_ForcesRequired &&
	   (count==1 || count % stepIncrementsGlobal_Forces == 0) ){
		stringstream number;
		number << count;
		string filename;
		if(BeforeDamage)
			filename = number.str() + ".pre.Global_Forces" ;
		else
			filename = number.str() + ".Global_Forces" ;
		ost = filemanager->OpenOutputStream(filename);
		int    dofForMaxResidual;
		double maxResidual;
		double * mappedResidualVector = new double[totalNumDof];
	
		getGlobalForces(internalForce, currentSolution,mappedResidualVector, maxResidual,dofForMaxResidual );
		mesh->printNodalData(internalForce, "displacements", ost);
		delete [] mappedResidualVector;
		filemanager->CloseOutputStream(ost);
	}

//checking if elemental output is needed
	bool doElementLoop=false;
	bool printNodalLCS_stress=false;
	bool printNodalGCS_stress=false;
	bool printQuad_stress=false;
    bool printQuadPoints=false;

	ofstream *NodalLCS_stress_ost;
	ofstream *NodalGCS_stress_ost;
	ofstream *Quad_stress_ost;
    ofstream *QuadPoints_ost;

	if(NodalLCS_stressRequired &&
		(count==1 || count % stepIncrementsNodalLCS_stress == 0) ){
		doElementLoop=true;
		printNodalLCS_stress=true;
		stringstream number;
		number << count;
		string filename;
		if(BeforeDamage)
			filename = number.str() + ".pre.Nodal.LCS.stress" ;
		else
			filename = number.str() + ".Nodal.LCS.stress" ;
		NodalLCS_stress_ost = filemanager->OpenOutputStream(filename);
		(*NodalLCS_stress_ost)<<"stress"<<endl;
	}
	if(NodalGCS_stressRequired &&
		(count==1 || count % stepIncrementsNodalGCS_stress == 0) ){
		doElementLoop=true;
		printNodalGCS_stress=true;
		stringstream number;
		number << count;
		string filename;
		if(BeforeDamage)
			filename = number.str() + ".pre.Nodal.GCS.stress" ;
		else
			filename = number.str() + ".Nodal.GCS.stress" ;
		NodalGCS_stress_ost = filemanager->OpenOutputStream(filename);
		(*NodalGCS_stress_ost)<<"stress"<<endl;
	}
	if(Quad_stressRequired && 
		(count==1 || count % stepIncrementsQuad_stress == 0) ){
		doElementLoop=true;
		printQuad_stress=true;
		stringstream number;
		number << count;
		string filename;
		if(BeforeDamage)
			filename = number.str() + ".pre.Quad.stress" ;
		else
			filename = number.str() + ".Quad.stress" ;
		Quad_stress_ost = filemanager->OpenOutputStream(filename);
	}
    if(QuadPointsRequired && 
        //Only run once per model (This data doesn't change during analysis)
		(count==1) ){
		doElementLoop=true;
		printQuadPoints=true;
		string filename="QuadPoints";
		QuadPoints_ost = filemanager->OpenOutputStream(filename);
	}
    if(Binary_Damage_File_Required){
		doElementLoop=true;
    }
	if (doElementLoop==false)
		goto SKIP_UPDATE_F_LOOP;


//writing elemental output
	GaussPointList  gaussPointList;
	int totalNumIPs;
	equations.zeroResultantVector();	
	for(i=0; i<numElements; i++)
	{
		e=(DMElasticElement*)(&mesh->element[i]);
		int group = e->getMaterialNumber();		
		e->update("F", currentSolution);
		e->getQuadraturePoints(gaussPointList);
		totalNumIPs=gaussPointList.totalNumIPs;
		if(printNodalGCS_stress){
			e->printNodalStressesInGCS(*NodalGCS_stress_ost);
		}
		if(printNodalLCS_stress){
			e->printNodalStressesInLCS(*NodalLCS_stress_ost);
		}
		if(printQuad_stress){
			e->printQuadStresses(totalNumIPs, *Quad_stress_ost);
		}
        if(printQuadPoints){
            (*QuadPoints_ost) << e->elementNumber << "\t" << group << "\t" << e->getTotalNumberOfIPs() << endl;
            e->printQuadPoints(totalNumIPs, *QuadPoints_ost);
        }
	}
	if(printNodalGCS_stress)	filemanager->CloseOutputStream(NodalGCS_stress_ost);
	if(printNodalLCS_stress)	filemanager->CloseOutputStream(NodalLCS_stress_ost);
	if(printQuad_stress)		filemanager->CloseOutputStream(Quad_stress_ost);

SKIP_UPDATE_F_LOOP:
	//not printing damage state if BeforeDamage==true becuase it will be the same as the damage state from the previous load step
	if(Damage_Factors_Required && BeforeDamage==false &&
	  (count==1 || count % stepIncrementsDamage_Factors == 0) ){
		stringstream number;
		number << count;
		string filename;
		filename = number.str() + ".Damage.Factors" ;
		ost = filemanager->OpenOutputStream(filename);
		printDamageState(ost);
		filemanager->CloseOutputStream(ost);
	}

    if(Volume_AveragesRequired) printVolumeAverageData();

	//binary output is used by the ExtractDV utility program.
	if(!BeforeDamage && Binary_Damage_File_Required) printDamageStateBinary();
}
//==========================================================
void DamageModel::printDamageStateBinary()
{
	eDamage->write((char*) &count, sizeof(int));
    for (int i=0; i<mesh->numElements; i++) {
		DMElasticElement* e=(DMElasticElement*)(&mesh->element[i]);
		e->printDamageStateBinary(eDamage);
    }
}
//==========================================================
void DamageModel::printDamageState(ostream *ost)
{
	for(int i=0; i<mesh->numElements; i++) 
	{
		DMElasticElement*e=(DMElasticElement*)(&mesh->element[i]);
		e->printDamageState(ost);
	}
}
//==========================================================
void DamageModel::LoadDamageState(istream *is)
{
	for(int iele =0;iele<mesh->numElements;iele++) { 
		DMElasticElement*e=(DMElasticElement*)(&mesh->element[iele]);
		e->LoadDamageState(is);
	}
}
//==========================================================
void DamageModel::InitializeOptionalOutput()
{
	eMF=filemanager->OpenBinaryOutputStream("_stiff");
	eQuadInfo =filemanager->OpenBinaryOutputStream("_quadInfo");
	eVolDist=filemanager->OpenBinaryOutputStream("volDist");
    if(Binary_Damage_File_Required){
	    eDamage=filemanager->OpenBinaryOutputStream("Damage.dat");
    }
    if(Volume_AveragesRequired && !doRestart){
        volAvgStressFile=filemanager->OpenOutputStream("VolumeAverageStress");
        volAvgStrainFile=filemanager->OpenOutputStream("VolumeAverageStrain");
    }
    if(Volume_AveragesRequired && doRestart){
        // This is a little bit of a special situation.
        // We want to retain all the data in the file up until the 
        // restart iteration.
        // First, make sure that the volume average files have been declared:
        // If it has, copy everything from that file up to the restart point.
        if(restartVolAvgStrainFile.compare("")==0){
            cout << "No Volume Average Strain file defined for restart - creating a new file" << endl;
            volAvgStrainFile=filemanager->OpenOutputStream("VolumeAverageStrain");
        }else{
            string TempStrainVals("");
            string PreLine,PostLine;
            ifstream VolAvgStrainIn(restartVolAvgStrainFile.c_str());
            getline(VolAvgStrainIn,PostLine);
            TempStrainVals+=PostLine + '\n';
            //Read in two lines per loadstep...
            for(int i=1;i<restart_iterationNumber;++i){
                getline(VolAvgStrainIn,PreLine);
                getline(VolAvgStrainIn,PostLine);
                TempStrainVals+=PreLine + '\n' + PostLine + '\n';
            }
            VolAvgStrainIn.close();
            volAvgStrainFile=filemanager->OpenOutputStream("VolumeAverageStrain");
            *volAvgStrainFile << TempStrainVals;
        }
        if(restartVolAvgStressFile.compare("")==0){
            cout << "No Volume Average Stress file defined for restart - creating a new file" << endl;
            volAvgStrainFile=filemanager->OpenOutputStream("VolumeAverageStress");
        }else{
            string TempStressVals("");
            string PreLine,PostLine;
            ifstream VolAvgStressIn(restartVolAvgStressFile.c_str());
            getline(VolAvgStressIn,PostLine);
            TempStressVals+=PostLine + '\n';
            //Read in two lines per loadstep...
            for(int i=1;i<restart_iterationNumber;++i){
                getline(VolAvgStressIn,PreLine);
                getline(VolAvgStressIn,PostLine);
                TempStressVals+=PreLine + '\n' + PostLine + '\n';
            }
            VolAvgStressIn.close();
            volAvgStressFile=filemanager->OpenOutputStream("VolumeAverageStress");
            *volAvgStressFile << TempStressVals;
        }

    }

	for(int iele =0;iele<mesh->numElements;iele++) { 
		DMElasticElement*e=0;
		e=(DMElasticElement*)(&mesh->element[iele]);
		e->eQuadInfo=eQuadInfo;
	}

	eVolDist->seekp(0,ios::beg);
    if(Binary_Damage_File_Required){
	    eDamage->seekp(0,ios::beg);
    }

	// XTANG 991018: Output the first line of OUT for iteration info
	BETA_OUT <<"Loadstep#\tmainTimeCount\tLoadFactor\tStrengthRatio\tIterations\t";
	BETA_OUT <<"jamming" << endl;

#define FORMAT setw(15)<<setprecision(6)

	SET_SCIENTIFIC(BETA_OUT);

#undef FORMAT
}
//==========================================================
void DamageModel::FinalizeOptionalOutput()
{
	filemanager->CloseOutputStream(eMF);
	filemanager->CloseOutputStream(eQuadInfo);
	filemanager->CloseOutputStream(eVolDist);
	filemanager->CloseOutputStream(eDamage);
    if(Volume_AveragesRequired){
        filemanager->CloseOutputStream(volAvgStressFile);
        filemanager->CloseOutputStream(volAvgStrainFile);
    }

    // delete "_stiff.dat" and "_quadinfo.dat"
	filemanager->ChangeCurrentDirectoryToInputDirectory();
	deleteFiles(filemanager->GetOutputFilenameWithRelativePath("_stiff"));
	deleteFiles(filemanager->GetOutputFilenameWithRelativePath("_quadInfo"));

}
//======================================================================
void DamageModel::LoadRestartData()
{
	LoadDamageState(initialDamageFile);
}
//========================================================================
void DamageModel::LoadDamageState(string filename)
{
	ifstream *is=filemanager->OpenInputStream(filename);
	if(is==0){
		EXIT_BETA("Cannot open initialDamageFile");
	}else
		LoadDamageState(is);
	filemanager->CloseInputStream(is);
}
//========================================================================
void DamageModel::updateFailureIndices()
{
	for(int i=0;i<mesh->numElements;i++) { 
		mesh->element[i].update("R",currentSolution); // R: check strength value only
	}
}
//========================================================================
bool DamageModel::UpdateISV()
{
    if(parallelAssembly_ESA){
        return UpdateISV_Parallel();
    }
	bool hasNewDamage=false;
	eQuadInfo->seekp(0,ios::beg);

	DMElasticElement *e=0;

	for(int i=0; i<mesh->numElements;i++) { 
//		mesh->element[i].update("VQ",currentSolution); //'Q' is for printing eQuadInfo file. uncomment if needed
		mesh->element[i].update("V",currentSolution);
		e=(DMElasticElement*)&mesh->element[i];
		hasNewDamage = hasNewDamage | e->hasNewDamage;		
	}
	return hasNewDamage;
}
//========================================================================
bool DamageModel::UpdateISV_Parallel()
{
    bool hasNewDamage=false;
    omp_set_num_threads(BETA_num_threads_for_assembly);
    #pragma omp parallel shared(hasNewDamage)
    {
        #pragma omp for reduction(|:hasNewDamage) 
        for(int i=0; i<BETA_num_threads_for_assembly; i++)
        {
            threadWorkspace *thread=&(assembler->threadWorkspaceList[i]);
            DMElasticElement *e=0;	
            for(int j=thread->minElementCalculationNumber;j<=thread->maxElementCalculationNumber;j++) {
                e = (DMElasticElement*) &mesh->element[j];
                e->update("V",currentSolution);
                hasNewDamage |= e->hasNewDamage;
            }
        }
    }
    return hasNewDamage;
}
//========================================================================
void DamageModel::outputVolumeDistribution()
{
	stringstream number;
	number << count;
	string filename= number.str() + "vA.dat";

    eQuadInfo->seekg(0,ios::beg);
    outVolumeDistribution(7,filename.c_str(),count,iIter);  // check ISV for this number
}
//========================================================================
void DamageModel::SolveLinearStep()
{
	ApplyLoadFactorToLoads();
	assembler->assembleKandF(NULL);//accounts for thermal loads
    equations.solve();
}
//========================================================================
void DamageModel::ApplyLoadFactorToLoads(){
	equations.zeroLoadVector();
    for(int i=0; i<(int)loads.size(); i++){ // Send loads to equation object
        if(loads[i] == NULL) break;
        if (loads[i]->isIncremental()) { 
            loads[i]->setRate(LoadFactor);
        }
        loads[i]->applyTo(equations,&mesh->node,mesh->numNodes, this);
	}
}
//========================================================================
void DamageModel::SetDegradationModelType(istream * inStream)
{
	int first,last,increment;

	int numMats=materialList.getNumElements();
	string DegradationModelName;
	int ParamStartIndex=0;

	char *localTokenList[20];
	int  numberOfTokens;

	while (true) {
		if( getLineAndTokenize(inStream,"exitSetDegradationModelType",localTokenList, numberOfTokens)==1) {
			goto ExitSetDegradationModelType;
		}
		if(numberOfTokens>=2 && COMPARE(localTokenList[0],"all")==0 ){
			first=0;
			last=numMats-1;
			increment=1;
			ParamStartIndex=1;
		}else if(numberOfTokens>=4){
			first=atoi(localTokenList[0]);
			last=atoi(localTokenList[1]);
			increment=atoi(localTokenList[2]);
			ParamStartIndex=3;
		}else{
			EXIT_BETA("Incorrect syntax for SetDegradationModelType");
		}

		DegradationModelName=localTokenList[ParamStartIndex];
		ParamStartIndex++;

		BETA_OUT <<first<<" "<<last<<" "<<increment<<endl;
		BETA_OUT << "Setting degradation type for materials "
				<< first <<" to "<<last << " by " << increment << '\n';
		BETA_OUT << "Creating Degradation Model: " << DegradationModelName << endl;
		

		if ( (first < 0) || (last >= numMats) )
			EXIT_BETA("You are trying to set the degradation type for a non-existant material!");

		//create the degradation model
		for(int i=first; i<=last; i+=increment) {
			dmElasticMaterial* mat=(dmElasticMaterial*)(&materialList[i]);
			mat->DegradationModel=CreateDegradationModel(DegradationModelName, filemanager);
			mat->DegradationModel->setFileManager(filemanager);
			mat->initializeDegradationModel(localTokenList, ParamStartIndex, numberOfTokens);
		}
	}
ExitSetDegradationModelType:
	BETA_OUT << "Finished SetDegradationModelType" << endl;
	CreateElementISV();
	BETA_OUT << "Finished creating ISV for elements" << endl;
}
//========================================================================
bool DamageModel::checkModel()
{
	BasicModel::checkModel();

	//This function checks to ensure that degradation models have been assigned to the
	//dmElasticMaterial object. If they are not assigned the code will crash in the
	//deallocation stage (without this check).

	int numMats = materialList.getNumElements();

	for(int i=0;i<numMats;i++){
		dmElasticMaterial *mat = (dmElasticMaterial*) &materialList[i];
		char *type = mat->getGroupType();
		int groupNum = mat->getGroupNum();

		if((groupNum != 0) && (type[0]!=0)) { //this check is needed since material group 0 is always part
							     //of the list even if no group 0 is specified
			if(!(mat->DegradationModel)){
				BETA_OUT<<"Degradation model not assigned to Material "<<i<<" ... exiting."<<endl;
				exit(1);}
		}
	}

	return true;
}
