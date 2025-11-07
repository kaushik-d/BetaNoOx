#include "stdafx.h"

#include <ctime>
#include <string>
#include <sstream>
#include "DiffusionModel.hpp"

#include "utility/directoryManipulation.hpp"

#include "factory/Factory.hpp"
#include "utility/formWriter.hpp"
#include "math/gaussQuadrature/gaussPointList.hpp"
#include "mesh/MeshUtility.h"
#include "mesh/NodeSet.hpp"

#include "BCs/SurfaceLoad.hpp"

#include "../element/TransientElementWorkspace.hpp"
#include "../element/HeatTransferElement3DTransient.hpp"
#include "../element/HeatTransferElement1DTransient.hpp"
//=========================================================================
//-------------------------------------------------------------------
extern int verboseFlag;
//-------------------------------------------------------------------
extern Factory	*factory;
void setVectorFromNodalValueFile(char* filename, double *vector, int nn, NodeGroup &node);

//-------------------------------------------------------------------
CreateErrorHandler(DiffusionModel);
//--------------------------------------------------------------------
DiffusionModel::DiffusionModel(void)
{
	AnalysisParameters.reserve(10);
	TimeSegmentList.reserve(10);

	time_integration_alpha=1.0; // backward differentiation scheme (default)
	TransientTimeUnit=minutes;

	ForceNonNegativeSolution=false;

	mainTimeStepCount=0;
	totalTimeSteps=0;
	mainTimeCount=0.0;
	printTimeCount=0.0;
	interpolationFactor = NULL;

	u=0;
    u_s=0;
    ut_s=0;
	delta_u=0;
	u_s1=0;

	delta_delta_u=0;
	delta_u_temp = 0;

	adap_tol = 5.0e-3;
	initial_ts = 1.0e-6;
	usingTimeStepControl = false;
	newTimeStepSize = 0.0;

	needsRefactorize = true;
	iterationTol = 1.0e-8;

	//concentration optional output
	concentration_out=false;
	concentration_out_Step_Interval=1;
	ConcentrationAnimationInput=0;
	PlotModelAvgConcentrationInput=0;
	//flux optional output
	flux_out=false;
	flux_out_Step_Interval=1;
	FluxAnimationInput=0;
	//volume averages
	volavg_out=false;
	volavg_out_Step_Interval=1;
	//actual conc
	actualConcentration_out=false;
	actualConcentration_out_Step_Interval=1;
	//Global Forces
	GlobalForces_out=false;
	GlobalForces_out_Step_Interval=1;
	
	RunFakeAnalysisToGenerateOutputFiles=false;

}
//--------------------------------------------------------------------
DiffusionModel::~DiffusionModel(void)
{
	ReleaseSolutionVectors();
}
//--------------------------------------------------------------------
int DiffusionModel::processCommand(char** localTokenList,ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;
	char newFileName[80];

	char token[200];
	int  foundMatch=0;

	strcpy(token,localTokenList[0]);
	BETA_OUT<<"\nCommand processor = 'DiffusionModel::processCommand'     "<<token<<endl;

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
//	IS_COMMAND(DoTransientAnalysis)  //retired
	IS_COMMAND(DoTransientNonLinearAnalysis)  
	IS_COMMAND(DoEnhancedTransientNonLinearAnalysis)  
	IS_COMMAND(ReadTransientAnalysisParameters)  	

	IS_COMMAND(DoTransientNonLinearAnalysisAdaptive)
	IS_COMMAND(DoTransientNonLinearAnalysisDecomposition)

	if(COMPARE(token,"ForceNonNegativeSolution") ==0){
		if( atoi(localTokenList[1])==0) 
			ForceNonNegativeSolution=false;
		else
			ForceNonNegativeSolution=true;
		BETA_OUT << "ForceNonNegativeSolution " << ForceNonNegativeSolution << endl;
		OK
	}
	
	if(COMPARE(token,"setTimeIntegrationAlpha") ==0) 
	{ 
		time_integration_alpha = atof(localTokenList[1]);
		BETA_OUT << "Time integration scheme - alpha = " << time_integration_alpha << endl;
		OK
	}

	if(COMPARE(token,"setTimeUnit") ==0) 
	{ 
		if( COMPARE(localTokenList[1],"second")==0 || COMPARE(localTokenList[1],"seconds")==0) 
			TransientTimeUnit=seconds;
		else if( COMPARE(localTokenList[1],"minute")==0 || COMPARE(localTokenList[1],"minutes")==0) 
			TransientTimeUnit=minutes;
		else if( COMPARE(localTokenList[1],"hour")==0 || COMPARE(localTokenList[1],"hours")==0) 
			TransientTimeUnit=hours;
		else
			BETA_OUT << "Cannot recognize time unit." << endl;
			
		BETA_OUT << "Transient Analysis Time Unit = " << TransientTimeUnit << " - " << localTokenList[1] << endl;
		OK
	}

	if(COMPARE(token,"setAlphaAdaptiveTolerance") ==0) 
	{ 
		adap_tol = atof(localTokenList[1]);
		BETA_OUT << "Alpha Adaptive Tolerance :" << adap_tol << endl;
		OK
	}

	if(COMPARE(token,"setAlphaAdaptiveInitialTimeStep") ==0) 
	{ 
		initial_ts = atof(localTokenList[1]);
		BETA_OUT << "Alpha Adaptive Initial Time Step :" << initial_ts << endl;
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
//=========================================================================
bool DiffusionModel::createLoads(char* name, istream *inStream)
{
	bool exitLoop;
	Load *newLoad=0;

	cout << "\n\nIn creatloads in diffusion model";

	IS_LOAD(SurfaceLoad)	
	else
		return ElasticityModel::createLoads(name, inStream);
}
//======================================================================
/*
void DiffusionModel::OptionalOutput(istream * inStream)
{
	char *localTokenList[20];
	int  numberOfTokens;

	//reading the options
	while (true) {
		if( getLineAndTokenize(inStream,"exitOptionalOutput",localTokenList, numberOfTokens)==1) {
			break;
		}
		BETA_OUT<<"Optional Output = "<< localTokenList[0] <<endl;
		if(COMPARE(localTokenList[0],"stress")==0)
			OutputOptions |= localstress;
		else if(COMPARE(localTokenList[0],"strain")==0)
			OutputOptions |= localstrain;
		else if(COMPARE(localTokenList[0],"quadstress")==0)
			OutputOptions |= quadstress;
		else if(COMPARE(localTokenList[0],"quadstrain")==0)
			OutputOptions |= quadstrain;
		else if(COMPARE(localTokenList[0],"volumeAverage")==0)
			OutputOptions |= volumeAverage;
		else{
			cout << localTokenList[0] << endl;
			FatalError("Incorrect Optional Output Specification!");
		}
	}//end of while

	//setting up the files
	HeatTransferElement3Dtransient *e=(HeatTransferElement3Dtransient*)(&mesh->element[0]);
	//e->openOutputFiles(&filemanager, filemanager.OutputOptions);
	//ElasticityElement3D::openOutputFiles(&filemanager, OutputOptions);
	openOutputFiles(filemanager, OutputOptions);

	//output to files
	outputToFiles(OutputOptions);

	//closing files
	//e->closeOutputFiles(&filemanager, OutputOptions);
	closeOutputFiles(filemanager, OutputOptions);

}
//======================================================================
void DiffusionModel::outputToFiles(int flags)
{
	int i;
	double * theSolution=0;

	int totalNumIPs;
	GaussPointList  gaussPointList;

	theSolution=equations.getSolution();

	HeatTransferElement3Dtransient *e=0;

	equations.zeroResultantVector();	
	for(i=0; i<mesh->numElements; i++)
	{
		e=(HeatTransferElement3Dtransient*)(&mesh->element[i]);

		int group = e->getMaterialNumber();		
		e->update("F", theSolution);
		e->getQuadraturePoints(gaussPointList);
		totalNumIPs=gaussPointList.totalNumIPs;

		//if(flags & localstress){
		//	e->extrapolateStresses(totalNumIPs, *localstressout, *localstrainout);
		//}
	
		//if(flags & quadstress){
		//	e->printQuadStresses(totalNumIPs, *(quadstressout) );
		//}

		//if(flags & quadstrain){
		//	e->printQuadStrains(totalNumIPs, *(quadstrainout) );
		//}
	}//end of loop on elements

	if(flags & volumeAverage)
	{
		BETA_OUT << "Printing Volume Averages..." << endl;
		mesh->element[0].initializeSummary(); // since getglobal forces is not called, i have to initialize the summary here 
		outVolumeStrConstituents(theSolution,0, 0, 0, true);
		mesh->element[0].outputSummary(); //no sense in printing summary since it just got initialized
	}
}
*/
//======================================================
/* // retired function
void DiffusionModel::assembleK_F_s(double *u_s, double *ut_s, double a1,double a2, double timeStep)
{
	equations.zeroMatrix();
	//double FluxForce[60000];
	//for (int i=0; i<totalNumDof; i++) FluxForce[i]=0.0;

	for(int i=0;i<mesh->numElements;i++) {

		mesh->element[i].calculateDofList();
//		mesh->element[i].update("K", u_s);//commented out to test 'total' form of time stepping algorithm with 1d diffusion element
		mesh->element[i].update("KM", u_s);

        int *dofList = mesh->element[i].getDofList();
		int numDof= mesh->element[i].getNumDof();
		

		((HeatTransferElement3Dtransient *) &mesh->element[i])->formK_s(a1,a2, timeStep);

		((HeatTransferElement3Dtransient *) &mesh->element[i])->formF_s(a1,a2,u_s,ut_s, timeStep);
		
		Ke = mesh->element[i].getStiffnessMatrix();
        Fe = mesh->element[i].getForceVector();

		equations.addMatrix((*Ke),mesh->element[i].getDofList(),
				                  mesh->element[i].getNumDof(),
				                  Fe);

	//	for(int ii=0; ii<numDof; ii++)
	//	{ FluxForce[dofList[ii]] += Fe[ii];	}

	}
	// printNodalData(FluxForce, "Force  ", dispOut );
}
*/
//======================================================
void DiffusionModel::assembleK_F_s_alpha_family(double *u_s, double *delta_u, double *u_s1, Array<double> &AnalysisParameters, double timeStep)
{
	equations.zeroMatrix();

	for(int i=0;i<mesh->numElements;i++) {

		mesh->element[i].calculateDofList();
		mesh->element[i].update("KM", u_s1); //latest solution

        int *dofList = mesh->element[i].getDofList();
		int numDof= mesh->element[i].getNumDof();

		((HeatTransferElement1Dtransient *) &mesh->element[i])->formK_s(AnalysisParameters, timeStep);

		((HeatTransferElement1Dtransient *) &mesh->element[i])->formF_s_new(AnalysisParameters, u_s, delta_u, u_s1, timeStep); //latest solution
		
		Matrix *Mbar_e = ((HeatTransferElement1Dtransient *) &mesh->element[i])->getKe_orig();
		//Matrix *Ke = mesh->element[i].getStiffnessMatrix();
        double *Fbar_e = mesh->element[i].getForceVector();

		equations.addMatrix((*Mbar_e),mesh->element[i].getDofList(),
				                  mesh->element[i].getNumDof(),
				                  Fbar_e);

	}
}
//======================================================
void DiffusionModel::assembleF_s_alpha_family(double *u_s, double *delta_u, double *u_s1, Array<double> &AnalysisParameters, double timeStep)
{
	for(int i=0;i<mesh->numElements;i++) {

		mesh->element[i].calculateDofList();
		mesh->element[i].update("F", u_s1); //latest solution

        int *dofList = mesh->element[i].getDofList();
		int numDof= mesh->element[i].getNumDof();

		((HeatTransferElement1Dtransient *) &mesh->element[i])->formF_s_new(AnalysisParameters, u_s, delta_u, u_s1, timeStep); //latest solution
		
        double *Fe = mesh->element[i].getForceVector();

		equations.addToLoadVector(mesh->element[i].calculateDofList(),
				              mesh->element[i].getNumDof(),
				              Fe);

	}
}
//======================================================
void DiffusionModel::assembleK_alpha_family(double *u_s, double *delta_u, double *u_s1, Array<double> &AnalysisParameters, double timeStep)
{
	equations.zeroMatrix();

	for(int i=0;i<mesh->numElements;i++) {

		mesh->element[i].calculateDofList();
		mesh->element[i].update("KM", u_s1); //latest solution

        int *dofList = mesh->element[i].getDofList();
		int numDof= mesh->element[i].getNumDof();

		((HeatTransferElement1Dtransient *) &mesh->element[i])->formK_s(AnalysisParameters, timeStep);

		//((HeatTransferElement1Dtransient *) &mesh->element[i])->formR_s_new(a1, a2, u_s, delta_u, u_s1, timeStep); //latest solution
		Matrix* Mbar_e = ((HeatTransferElement1Dtransient *) &mesh->element[i])->getKe_orig();
		//Matrix *Ke = mesh->element[i].getStiffnessMatrix();
        //Fe = mesh->element[i].getForceVector();

		equations.addMatrix((*Mbar_e),mesh->element[i].getDofList(),
				                  mesh->element[i].getNumDof(),
				                  NULL);

	}
}
//======================================================
void DiffusionModel::DoTransientAnalysis(istream * inStream)
{
	BETA_OUT << "if you want to run XT's transient moisture analysis, use beta source prior to 2/7/2008" << endl;	 
} // end of DoTransientAnalysis
//===========================================================
#define FF1 setw(7)<<setprecision(4)<<setiosflags(ios::fixed)

void DiffusionModel::printNodalDataByElement(double *data, char * label, BasicMesh *mesh, ostream *outStream )
{
	int numNodesPerElement=mesh->element[0].numNodesPerElement;
	InterpAndSolDerivatives       id(numNodesPerElement,numNodesPerElement); id.allocate();
	double *elemDisp[3];
	
	//------------------------------------------------
	elemDisp[0] = id.nodalValues_u1; 
	elemDisp[1] = id.nodalValues_u2;
	elemDisp[2] = id.nodalValues_u3;     
    double saturationValue;

	(*outStream) <<label<<endl;

	int iele, inode,j;
	for (iele=0; iele<mesh->numElements; iele++) {
		BasicElement &ae=mesh->element[iele];
		
		id.numberOfInterp = ae.numNodesPerElement;
		ae.extractSolution( data, elemDisp);
        saturationValue = ae.getNodalScalingFactor();
		// here we should notice that dof at each node may be different
		(*outStream) << iele <<" "<<ae.getMaterialNumber()<<endl;
		for(inode=0; inode<ae.numNodesPerElement; inode++){
			for( j=0; j<ae.getNumDofatNode(inode); j++ ) {
				(*outStream)<<FF1<<elemDisp[j][inode]*saturationValue;   
			}
            (*outStream)<<endl;
		}
	}
}
//======================================================
void DiffusionModel::PrintFluxes(ostream &SO, double *sol)
{
	GaussPointList  gaussPointList;
	int totalNumIPs;

	SO <<"stress"<<endl;
	//bool FlagNodalStressInLCS = false;

	HeatTransferElement3Dtransient* e=0;
	for(int i=0;i<mesh->numElements;i++) {
		e=(HeatTransferElement3Dtransient*)(&mesh->element[i]);
		e->update("F", sol);
		e->getQuadraturePoints(gaussPointList);
		totalNumIPs=gaussPointList.totalNumIPs;
		e->printNodalFluxesInGCS(SO);
	}
}
//======================================================
void DiffusionModel::PrintQuadFluxes(ostream &SO, double *sol)
{
	GaussPointList  gaussPointList;
	int totalNumIPs;

	//SO <<"stress"<<endl;
	bool FlagNodalStressInLCS = false;

	HeatTransferElement3Dtransient* e=0;
	for(int i=0;i<mesh->numElements;i++) {
		e=(HeatTransferElement3Dtransient*)(&mesh->element[i]);
		e->update("F", sol);
		e->getQuadraturePoints(gaussPointList);
		totalNumIPs=gaussPointList.totalNumIPs;
		e->PrintQuadFluxes(SO,totalNumIPs);
	}
}
//======================================================
void DiffusionModel::PrintQuadPoints(ostream &SO)
{
	GaussPointList  gaussPointList;
	int totalNumIPs;

	HeatTransferElement3Dtransient* e=0;
	for(int i=0;i<mesh->numElements;i++) {
		e=(HeatTransferElement3Dtransient*)(&mesh->element[i]);
		e->getQuadraturePoints(gaussPointList);
		totalNumIPs=gaussPointList.totalNumIPs;
		SO << e->elementNumber << "\t" << e->getMaterialNumber() << "\t" << e->getTotalNumberOfIPs() << endl;
		e->printQuadPoints(totalNumIPs,SO);
	}
}
//===================================================================================
void DiffusionModel::PrintActualConcentrations(ostream &SO, double *sol)
{
	SO <<"stress"<<endl;
	bool FlagNodalStressInLCS = false;

	HeatTransferElement3Dtransient* e=0;
	for(int i=0;i<mesh->numElements;i++) {
		e=(HeatTransferElement3Dtransient*)(&mesh->element[i]);
		SO << e->elementNumber << '\t' << e->getMaterialNumber() << endl;
		e->PrintActualConcentrations(sol, SO);
	}
}
//===================================================================================
void DiffusionModel::CreateWorkspace(ElementWorkspace * &elemWorkspace)
{
	if(elemWorkspace==0) 
		elemWorkspace=new TransientElementWorkspace;
}
//=======================================================================
void DiffusionModel::allocateWorkspace(ElementWorkspace * &elemWorkspace)
{
	int WorkspaceMaxNumberDof	=elemWorkspace->WorkspaceMaxNumberDof;

	TransientElementWorkspace* workspace=(TransientElementWorkspace*)elemWorkspace;
	workspace->TransientTimeUnit=TransientTimeUnit;

	workspace->F1			=new double [WorkspaceMaxNumberDof];
	workspace->F2			=new double [WorkspaceMaxNumberDof];
	workspace->M_deltaq		=new double [WorkspaceMaxNumberDof];

	BasicModel::allocateWorkspace(elemWorkspace);
}
//========================================================================
void DiffusionModel::AllocateSolutionVectors()
{
	u = equations.getSolution();  // unmapped
    u_s = new double[totalNumDof];//solution from previous timestep
	ut_s= new double[totalNumDof];
	delta_u= new double[totalNumDof]; 
	u_s1= new double[totalNumDof]; 
	delta_u_temp = new double [totalNumDof];
	delta_delta_u = new double [totalNumDof];

	for (int i=0;i<totalNumDof;i++) {u_s[i]=u[i]=ut_s[i]=0.0; u_s1[i]=delta_u[i]=delta_delta_u[i]=delta_u_temp[i]=0.0;}

}
//========================================================================
void DiffusionModel::ReleaseSolutionVectors()
{
	if(u_s)		{delete [] u_s;		u_s=0;}
	if(ut_s)	{delete [] ut_s;	ut_s=0;}
	if(u_s1)	{delete [] u_s1;	u_s1=0;}
	if(delta_u)	{delete [] delta_u;	delta_u=0;}
}
//========================================================================
void DiffusionModel::SetupAnalysisParameters()
{
	double *a=0;

	//adding parameters for a1 and a2:
	//a1=AnalysisParameters[0] and
	//a2=AnalysisParameters[1]

	//this function only allocates memory for the analysis parameters.
	//the actual value will be written to the variable later on.
	a=new double(0.0);
	AnalysisParameters.add(*a);
	a=new double(0.0);
	AnalysisParameters.add(*a);
}
//===========================================================
bool DiffusionModel::ReadTransientAnalysisParameters(istream * inStream)
{
	SetupAnalysisParameters();
	createAssembler();
	allocateForAnalysis();
	AllocateSolutionVectors();
	ReadTransientAnalysisCommands(inStream);
	return true;
}
//===========================================================
bool DiffusionModel::ReadTransientAnalysisCommands(istream * inStream)
{
	char *localTokenList[20];
	int  numberOfTokens;

	double dump;
	int i;
	double iTimeCount=0.0; // for calculating sqrt stepping
	totalTimeSteps=0;


	while(true) {//loop reading each time-integration 'segment' command or initial solution command

_AGAIN:
		if( getLineAndTokenize(inStream,"exitReadTransientAnalysisParameters",localTokenList, numberOfTokens)==1) {	
			return true;
		}

		// 07072003 add code specifying initial conditions in heat transfer problems
		// currently only uniform initial condition is allowed
		if (COMPARE(localTokenList[0],"initialConditions")==0) {
			dump=atof(localTokenList[1]);
	        for (i=0;i<totalNumDof;i++) u_s[i]=dump;
			goto _AGAIN;

		} else if (COMPARE(localTokenList[0],"SetInitialSolution")==0) {
			//format: SetInitialSolution nodenum dir value 
			int dumpNode=0,dumpNodeDOF=0;
			dumpNode=atoi(localTokenList[1]);
			dumpNodeDOF=atoi(localTokenList[2]); // dof 'direction' starting with 1
			if( dumpNodeDOF > mesh->node[dumpNode].getNumDof() )
				FatalError("Error when specifying parameters for SetInitialSolution");
			int iDOF = mesh->node[dumpNode].getFirstDof() + dumpNodeDOF - 1;
			dump=atof(localTokenList[3]);
			u_s[iDOF]=dump;
			goto _AGAIN;

		} else if (COMPARE(localTokenList[0],"SetSolutionForNodeGroup")==0) {
			//format: SetSolutionForNodeGroup nodegroupname dir value 
			int dumpNode=0,dumpNodeDOF=0;
			char dumpGroupName[256];
			strcpy(dumpGroupName,localTokenList[1]);
			NodeSet nlist=mesh->getNodeSet(dumpGroupName);
			int numNodes=nlist.getNumNodes();
			Node *n=0;
	
			dumpNodeDOF=atoi(localTokenList[2]); // dof 'direction' starting with 1
			if( dumpNodeDOF > mesh->node[dumpNode].getNumDof() )
				FatalError("Error when specifying parameters for SetInitialSolution");
			dump=atof(localTokenList[3]);

			for(i=0;i<numNodes;i++){
				n=nlist[i];
				if(n->numDof <= 0 || (dumpNodeDOF > n->getNumDof() ) ){
					BETA_OUT << "Node has no dofs. cannot apply solution to this node - node num:" << n->getNodeNum() << endl;
				}else{
					int iDOF = n->getFirstDof() + dumpNodeDOF - 1;
					u_s[iDOF]=dump;
				}
			}
			goto _AGAIN;

		} else {

			if (numberOfTokens < 3) {
				BETA_OUT << " Error in defining control parameters for the trasient analysis! \n";	 
				BETA_OUT << " Format: timeSteppingMethod    timeStepSize     numTimeSteps \n";	 
				exit(1);
			}

 			int timeSteppingMethod=atoi(localTokenList[0]);
			double inputTimeStepSize=atof(localTokenList[1]);
			int numTimeSteps=atoi(localTokenList[2]);
			totalTimeSteps+=numTimeSteps;

			Array<double> *timesegment=0;
			timesegment=new Array<double>;
			(*timesegment).reserve(numTimeSteps);
			double *ts=0;
			
			timeStepSize=0.0;
			for (int iCount=0; iCount<numTimeSteps; iCount++) {
				if(timeSteppingMethod==0) { // linear stepping
					timeStepSize=inputTimeStepSize;
				} else if (timeSteppingMethod==1) {// sqrt stepping for moisture diffusion 
					timeStepSize=(sqrt(iTimeCount)+inputTimeStepSize)*(sqrt(iTimeCount)+inputTimeStepSize)-iTimeCount;
					iTimeCount += timeStepSize;
				}
				ts=new double(timeStepSize);
				(*timesegment).add(*ts);
			}
			TimeSegmentList.add(*timesegment);
		}
	}
}
//========================================================================
char *getTimeUnit(BETA_Time_Unit TimeUnit)
{
	if(TimeUnit==hours) return (" hours");
	else if(TimeUnit==seconds) return (" seconds");
	else if(TimeUnit==minutes) return (" minutes");
	else return (" Undefined time unit");
}
//========================================================================
void DiffusionModel::WriteTimeStamp(double outTime)
{
	cout << "#" << mainTimeStepCount <<" Current Time= " << outTime << getTimeUnit(TransientTimeUnit) << ", ts size= " << timeStepSize << endl;
	if(verboseFlag == Max)	{
	BETA_OUT << "#" << mainTimeStepCount <<" Current Time= " << DECIMAL_FORMAT << outTime << getTimeUnit(TransientTimeUnit) << ", ts size= " << DECIMAL_FORMAT << timeStepSize << endl;
	}
}
//========================================================================
void DiffusionModel::getTimeStepLoopParametersForRestart(int restartTimeStepCount, int &iTimeSegmentStart, int &iCountStart)
{
	double rMainTimeCount=0;
	int rMainTimeStepCount=0;


	int numTimeSegments=TimeSegmentList.getNumElements();
	for (int iTimeSegment=0; iTimeSegment<numTimeSegments; iTimeSegment++) {//loop over timesegments
		int numTimeSteps=TimeSegmentList[iTimeSegment].getNumElements();
		for (int iCount=0; iCount<numTimeSteps; iCount++) {//loop over time steps
			timeStepSize=TimeSegmentList[iTimeSegment][iCount];
			rMainTimeCount+=timeStepSize;
			rMainTimeStepCount++;
			if(rMainTimeStepCount==restartTimeStepCount){
				iTimeSegmentStart=iTimeSegment;
				iCountStart=iCount;
				return;
			}
			mainTimeCount+=timeStepSize;
			mainTimeStepCount++;

			if(usingTimeStepControl){
				printTimeCount += timeStepSize;}
		}
	}
}
//========================================================================
void DiffusionModel::LoadRestartData()
{
	setVectorFromNodalValueFile(restart_dispfilename, u_s, mesh->numNodes, mesh->node);
}
//========================================================================
void DiffusionModel::DoTransientNonLinearAnalysis(istream * inStream)
{
	double * deltaSolution=0;
	equations.setUsePreviousFillInOrdering(true);
	bool DoIterate=false;
	bool isTimeStepComingFromCFD = true; // Hard coded now, Change to input from input script

	checkModel();
	InitializeOptionalOutput();
	
	mainTimeCount=0.0;
	mainTimeStepCount=0;

	int iTimeSegmentStart=0;
	int iCountStart=0;
	//restart capability for DoTransientNonLinearAnalysis(istream * inStream)
	if(doRestart && !restartInitiated){
		LoadRestartData();
		int restartTimeStepCount=restart_iterationNumber;
		getTimeStepLoopParametersForRestart(restartTimeStepCount, iTimeSegmentStart, iCountStart);
		restartInitiated=true;
	}

	BETA_OUT << "solution:" << endl;
	BETA_OUT << "\n Timestep# time solution " << endl;
	displayTime("Time to prepare transient analysis", BETA_OUT);

	int numTimeSegments=TimeSegmentList.getNumElements();
	for (int iTimeSegment=iTimeSegmentStart; iTimeSegment<numTimeSegments; iTimeSegment++) {//loop over timesegments
		int numTimeSteps=TimeSegmentList[iTimeSegment].getNumElements();
		for (int iCount=iCountStart; iCount<numTimeSteps; iCount++) {//loop over time steps

			timeStepSize=TimeSegmentList[iTimeSegment][iCount];

			AnalysisParameters[0] = time_integration_alpha*timeStepSize;
			AnalysisParameters[1] = (1-time_integration_alpha)*timeStepSize;
			
			mainTimeCount+=timeStepSize;
			mainTimeStepCount++;

            		if(RunFakeAnalysisToGenerateOutputFiles) goto _WriteOptionalOutput;
            		//displayTime("startup single step", BETA_OUT);
			// everytime clear rooms in equation xt: 102200
			
			equations.zeroLoadVector();
			equations.zeroMatrix();
			
			// Set up applied loadVector; Send loads to equation object
			for(int i=0; i<(int)loads.size(); i++){
	            		if (loads[i]->isIncremental())
	                		loads[i]->setRate(0.01);
				else
					loads[i]->setRate(1.0);  // uniform loading
				loads[i]->applyTo(equations,&mesh->node,mesh->numNodes, this);
			}
			
			if(isTimeStepComingFromCFD){
				// change the timeStepSize and recalculate AnalysisParameters
				// Currnetly this has to be done after applying SurfaceLoads
				mainTimeCount-=timeStepSize;
				timeStepSize = getTimeStepFromSharedMemory();
				AnalysisParameters[0] = time_integration_alpha*timeStepSize;
				AnalysisParameters[1] = (1-time_integration_alpha)*timeStepSize;
				mainTimeCount+=timeStepSize;
			}

			WriteTimeStamp(mainTimeCount);

			//apply additional equations: important for dummy load !
			equations.applyAdditionalEquations(); 

			//linear solution:
			sumVector(u_s1,	u_s, delta_u,totalNumDof);//why is this needed? appears that delta_u is always going to filled with zeros
			assembleK_F_s_alpha_family(u_s, delta_u, u_s1, AnalysisParameters, timeStepSize);

			deltaSolution=equations.solve(u_s1);
            		
			//displayTime("end single step", BETA_OUT);
			//equations.addIncrementalDisplacements(deltaSolution);

			//temporary forcing to zero.... ONLY for debug !!!!
			//concentration goes below zero
			if(ForceNonNegativeSolution){
				int FoundNegativeSolution=ConvertNegativeNumbersToZero(deltaSolution, totalNumDof);
				if(verboseFlag>=Min && FoundNegativeSolution != -1)
					cout << "Found negative solution. Converted to zero!" << endl;
			}


			copyVector(deltaSolution, delta_u, totalNumDof); //copy deltaU to local array
			sumVector(u_s1,		u_s,		delta_u,totalNumDof);

			TimeStepUpdate();//update Phi

            		//BETA_OUT << "\n" << iCount << " " << count << " " << FORMAT << u_s1[1];

			if(DoIterate) TransientAnalysisIterate();
            		
			_WriteOptionalOutput:
			WriteOptionalOutput();
                          
                        // Write solution to the shared memory
                        writeSolToShMemory(u_s1);

			//BETA_OUT << iCount << " " << count << " " << u_s1[1] << endl;

			initializeVector( delta_u, 0.0, totalNumDof);
			copyVector( u_s1, u_s, totalNumDof);
            		
			//displayTime("Time to finish timestep", BETA_OUT);
		}//end TimeStep loop
		BETA_OUT << endl;
	}//end TimeSegment loop
	BETA_OUT << endl;

	if(equations.getUsePreviousFillInOrdering() == true){
		equations.ReleaseMemory();
	}

	FinalizeOptionalOutput();
} // end of DoTransientNonLinearAnalysis
//
//========================================================================
//
void DiffusionModel::DoTransientNonLinearAnalysisAdaptive(istream * inStream)
{
	double * deltaSolution=0;

	usingTimeStepControl = true;

	equations.setUsePreviousFillInOrdering(true);
	bool DoIterate=false;

	checkModel();
	InitializeOptionalOutput();
	
	printTimeCount = 0.0;
	double printTS = 0.0;

	mainTimeCount=0.0;
	mainTimeStepCount=0;

	timeStepSize = initial_ts;

	int iTimeSegmentStart=0;
	int iCountStart=0;
	//restart capability for DoTransientNonLinearAnalysis(istream * inStream)
	if(doRestart && !restartInitiated){
		LoadRestartData();
		int restartTimeStepCount=restart_iterationNumber;
		getTimeStepLoopParametersForRestart(restartTimeStepCount, iTimeSegmentStart, iCountStart);
		restartInitiated=true;
	}

	BETA_OUT << "solution:" << endl;
	BETA_OUT << "\n Timestep# time solution " << endl;
	displayTime("Time to prepare transient analysis", BETA_OUT);

	int numTimeSegments=TimeSegmentList.getNumElements();
	for (int iTimeSegment=iTimeSegmentStart; iTimeSegment<numTimeSegments; iTimeSegment++) {//loop over timesegments
		int numTimeSteps=TimeSegmentList[iTimeSegment].getNumElements();
		for (int iCount=iCountStart; iCount<numTimeSteps; iCount++) {//loop over time steps
			printTS=TimeSegmentList[iTimeSegment][iCount];
		
			printTimeCount+=printTS;
			mainTimeStepCount++;
			
			WriteTimeStamp(printTimeCount);

if(RunFakeAnalysisToGenerateOutputFiles) goto _WriteOptionalOutput;

while (mainTimeCount<printTimeCount) {

			bool TS_SUCCESS = false;

			while(TS_SUCCESS == false) {

			AnalysisParameters[0] = time_integration_alpha*timeStepSize;
			AnalysisParameters[1] = (1-time_integration_alpha)*timeStepSize;

				// everytime clear rooms in equation xt: 102200
				equations.zeroLoadVector();
				equations.zeroMatrix();
			
				// Set up applied loadVector; Send loads to equation object
				for(int i=0; i<(int)loads.size(); i++){
				    if (loads[i]->isIncremental())
				        loads[i]->setRate(0.01);
					else
						loads[i]->setRate(1.0);  // uniform loading
					loads[i]->applyTo(equations,&mesh->node,mesh->numNodes, this);
				}

				// apply additional equations: important for dummy load !
				equations.applyAdditionalEquations(); 

				//linear solution:
				sumVector(u_s1,		u_s,		delta_u,totalNumDof);//why is this needed? appears that delta_u is always going to filled with zeros
				assembleK_F_s_alpha_family(u_s, delta_u, u_s1, AnalysisParameters, timeStepSize); 

	//			displayTime("Time to setup equations", BETA_OUT);

				deltaSolution=equations.solve(u_s1);

    //			displayTime("Time for solve", BETA_OUT);

				//equations.addIncrementalDisplacements(deltaSolution);

				//temporary forcing to zero.... ONLY for debug !!!!
				//concentration goes below zero
				if(ForceNonNegativeSolution){
					int FoundNegativeSolution=ConvertNegativeNumbersToZero(deltaSolution, totalNumDof);
					if(verboseFlag>=Min && FoundNegativeSolution != -1)
						cout << "Found negative solution. Converted to zero!" << endl;
				}

				//check delta U against adaptive tolerance
				TS_SUCCESS = checkSolutionAndModifyTimeStep(deltaSolution,timeStepSize);

			} // end of loop to check for successful timestep

			//mainTimeCount += timeStepSize;

			copyVector(deltaSolution, delta_u, totalNumDof); //copy deltaU to local array
			sumVector(u_s1,		u_s,		delta_u,totalNumDof);

			TimeStepUpdate();//update Phi

			timeStepSize = newTimeStepSize;

//			BETA_OUT << "\n" << iCount << " " << count << " " << FORMAT << u_s1[1];

			initializeVector( delta_u, 0.0, totalNumDof);
			copyVector( u_s1, u_s, totalNumDof);
}//end of loop to check if t<tprint

			//need to insert option to interpolate here
			interpolateSolutionForPrintTime(printTimeCount,deltaSolution);

			if(DoIterate) TransientAnalysisIterate();
_WriteOptionalOutput:
			WriteOptionalOutput();

			//BETA_OUT << iCount << " " << count << " " << u_s1[1] << endl;


//displayTime("Time to finish timestep", BETA_OUT);
		}//end TimeStep loop
		BETA_OUT << endl;
	}//end TimeSegment loop
	BETA_OUT << endl;

	if(equations.getUsePreviousFillInOrdering() == true){
		equations.ReleaseMemory();
	}

	FinalizeOptionalOutput();
} // end of DoTransientNonLinearAnalysis
//========================================================================
bool DiffusionModel::checkSolutionAndModifyTimeStep(double *deltaSolution,double &delta_t)
{
	//get the maximum absolute value of delta U
	double maxLTE = 0.0;
	for(int i=0;i<totalNumDof;i++)
	{
		double val = fabs(deltaSolution[i]);
		if(val > maxLTE)
		{
			maxLTE = val;
		}
	}

	//calculate new time step value
	//=====================================
	double minTS = initial_ts;
	if(delta_t < minTS)
	{
		delta_t = minTS;
	}
	//===================================== 
	double delta_t_old = delta_t;
	newTimeStepSize = 0.7 * delta_t * adap_tol / maxLTE;
	
	bool flag = false;

	if(maxLTE < adap_tol || delta_t_old == minTS) {
		flag = true;
		mainTimeCount += delta_t_old;
		BETA_OUT<<delta_t_old<<"     "<<maxLTE<<"          "<<"t = "<<mainTimeCount<<endl; //added bco 121609
	}
	else{
		delta_t = newTimeStepSize;	}

	return flag;
}
//========================================================================
void DiffusionModel::interpolateSolutionForPrintTime(double printTime, double *deltaSolution)
{
	interpolationFactor = (printTime - (mainTimeCount - timeStepSize))/timeStepSize;
	for(int i=0;i<totalNumDof; i++)
	{
	ut_s[i] = deltaSolution[i] * interpolationFactor + (u_s1[i] - deltaSolution[i]);

	}
}
//========================================================================
void DiffusionModel::DoTransientNonLinearAnalysisDecomposition(istream * inStream)
{
	double * deltaSolution=0;
	equations.setUsePreviousFillInOrdering(true);
	bool DoIterate=false;

	checkModel();
	InitializeOptionalOutput();
	
	mainTimeCount=0.0;
	mainTimeStepCount=0;

	int iTimeSegmentStart=0;
	int iCountStart=0;

	int localTimeStepCount=0;

	//restart capability for DoTransientNonLinearAnalysis(istream * inStream)
	if(doRestart && !restartInitiated){
		LoadRestartData();
		int restartTimeStepCount=restart_iterationNumber;
		getTimeStepLoopParametersForRestart(restartTimeStepCount, iTimeSegmentStart, iCountStart);
		restartInitiated=true;
	}

	BETA_OUT << "solution:" << endl;
	BETA_OUT << "\n Timestep# time solution " << endl;
	displayTime("Time to prepare transient analysis", BETA_OUT);

	int numTimeSegments=TimeSegmentList.getNumElements();
	for (int iTimeSegment=iTimeSegmentStart; iTimeSegment<numTimeSegments; iTimeSegment++) {//loop over timesegments
		int numTimeSteps=TimeSegmentList[iTimeSegment].getNumElements();
		for (int iCount=iCountStart; iCount<numTimeSteps; iCount++) {//loop over time steps

			timeStepSize=TimeSegmentList[iTimeSegment][iCount];

			AnalysisParameters[0] = time_integration_alpha*timeStepSize;
			AnalysisParameters[1] = (1-time_integration_alpha)*timeStepSize;
			
			mainTimeCount+=timeStepSize;
			mainTimeStepCount++;
			
			WriteTimeStamp(mainTimeCount);

if(RunFakeAnalysisToGenerateOutputFiles) goto _WriteOptionalOutput;

			if(needsRefactorize){
				//cout<<"SOLVING"<<endl;
				// everytime clear rooms in equation xt: 102200
				displayTime("Start Zero Matrix & Apply Loads/AE", BETA_OUT);

				equations.zeroLoadVector();
				equations.zeroMatrix();
			
				// Set up applied loadVector; Send loads to equation object
				for(int i=0; i<(int)loads.size(); i++){
				   if (loads[i]->isIncremental())
						loads[i]->setRate(0.01);
					else
						loads[i]->setRate(1.0);  // uniform loading
					loads[i]->applyTo(equations,&mesh->node,mesh->numNodes, this);
				}

				// apply additional equations: important for dummy load !
				equations.applyAdditionalEquations(); 

				//linear solution:
				displayTime("Start Sum Vector", BETA_OUT);
				sumVector(u_s1,		u_s,		delta_u,totalNumDof);//why is this needed? appears that delta_u is always going to filled with zeros
				displayTime("Start Assembly", BETA_OUT);
				assembleK_F_s_alpha_family(u_s, delta_u, u_s1, AnalysisParameters, timeStepSize); 

				needsRefactorize = false;
				equations.setUsePreviousFactorizedMatrixAtRuntime(true);

				displayTime("Start Solve", BETA_OUT);
				deltaSolution=equations.solve(u_s1);
				
				displayTime("Start Copy Vector", BETA_OUT);
				copyVector(deltaSolution, delta_u_temp, totalNumDof);
				displayTime("End full solve portion", BETA_OUT);
			}

			else {
				//cout<<"ITERATING"<<endl;
				displayTime("Start Iteration Portion", BETA_OUT);
				bool convergenceCheck = false;
				int iterationCount = 0;

				while(convergenceCheck == false) {
					CalculateResidual(u_s,delta_u_temp,u_s1,AnalysisParameters,timeStepSize);
					convergenceCheck = checkResidual();

					if(convergenceCheck || iterationCount >=3){
						break;
					}
					
					delta_delta_u = equations.solve(u_s1);
					sumVector(delta_u_temp,delta_u_temp,delta_delta_u,totalNumDof);

					++iterationCount;

				}
				BETA_OUT<<"# of Iterations:"<<iterationCount<<endl;
				checkIfNeedsRefactorization(localTimeStepCount,iterationCount);
				
			displayTime("End Iteration Portion", BETA_OUT);
			}

			localTimeStepCount++;

				//temporary forcing to zero.... ONLY for debug !!!!
				//concentration goes below zero
				if(ForceNonNegativeSolution){
					int FoundNegativeSolution=ConvertNegativeNumbersToZero(deltaSolution, totalNumDof);
					if(verboseFlag>=Min && FoundNegativeSolution != -1)
						cout << "Found negative solution. Converted to zero!" << endl;
				}


				copyVector(delta_u_temp, delta_u, totalNumDof); //copy deltaU to local array
				sumVector(u_s1,		u_s,		delta_u,totalNumDof);

				TimeStepUpdate();//update Phi
                               

				if(DoIterate) TransientAnalysisIterate();
	_WriteOptionalOutput:
				WriteOptionalOutput();

				initializeVector( delta_u, 0.0, totalNumDof);
				copyVector( u_s1, u_s, totalNumDof);
				
			}//end TimeStep loop
			BETA_OUT << endl;
		}//end TimeSegment loop
	BETA_OUT << endl;

	if(equations.getUsePreviousFillInOrdering() == true){
		equations.ReleaseMemory();
	}

	FinalizeOptionalOutput();
} // end of DoTransientNonLinearAnalysis
//========================================================================
int DiffusionModel::writeSolToShMemory(double * sol)
{
    key_t Infokey = 1234;
    key_t Tempkey = 9101;
    
    int InfoID;
    int * InfoPt;
    int TempID;
    double * TempPt;
    int i;

    // Locate the info segment, only difference from the server call
    // of the shmget function is the 3rd argument 

    InfoID = shmget(Infokey, 2*sizeof(int), 0666);

    // Check if info segment is located correctly
    
    if (InfoID  < 0) {
        printf("Info segment create error");
        return 1;
    }

    // Attach the info segment to the address space.

    InfoPt = (int *) shmat(InfoID, NULL, 0);

    // Check if info segment is attached correctly

    if ( InfoPt == (int *) -1) {
        printf("Info segment attach error");
        return 1;
    }

    // Locate the Temp segment, only difference from the server call
    // of the shmget function is the 3rd argument 

    TempID = shmget(Tempkey, InfoPt[1]*sizeof(double), 0666);

    // Check if Temp segment is located correctly
    
    if (TempID  < 0) {
        printf("Temp segment create error");
        return 1;
    }

    // Attach the Temp segment to the address space.

    TempPt = (double *) shmat(TempID, NULL, 0);

    // Check if temp segment is attached correctly

    if ( TempPt == (double *) -1) {
        printf("Temp segment attach error");
        return 1;
    }

	// get the total number of nodes from the shared memory segment
	int numNodes = InfoPt[1];

	// write the nodal temperatures to the shared memory
	for (int i = 0; i < numNodes; i++)
	{
		TempPt[i] = sol[i];
	}
	
	// temerature has been written, now change the Info status
	InfoPt[0] = 2;


    // Detach segments

    shmdt(InfoPt);
    shmdt(TempPt);

  return 0;
}
//========================================================================
//
double DiffusionModel::getTimeStepFromSharedMemory()
{

    key_t Infokey = 1234;
    key_t Timekey = 1213;
    
    int InfoID;
    int * InfoPt;
    int TimeID;
    double * TimePt;
    int i;

    // Locate the info segment, only difference from the server call
    // of the shmget function is the 3rd argument 

    do
    {
    	InfoID = shmget(Infokey, 2*sizeof(int), 0666);

    	// Check if info segment is located correctly
    
    	if (InfoID  < 0) {
        	cout << "Info segment locate error in getTimeStepFromSharedMemory"<<endl;
        	cout << "Trying to locate again in 1 second"<<endl;
		sleep(1);
	}
    }
    while (InfoID < 0);	

    // Attach the info segment to the address space.

    InfoPt = (int *) shmat(InfoID, NULL, 0);

    // Check if info segment is attached correctly

    if ( InfoPt == (int *) -1) {
        printf("Info segment attach error");
        return 1;
    }

    while (InfoPt[0] != 1)
    {
	cout <<"Time step is not set from CFD yet. Chacking again in 1 second";
	sleep(1);
    }

    // Locate the Time segment, only difference from the server call
    // of the shmget function is the 3rd argument 

    TimeID = shmget(Timekey, 1*sizeof(double), 0666);

    // Check if time segment is located correctly
    
    if (TimeID  < 0) {
        printf("Time segment locate error");
        return 1;
    }

    // Attach the Time segment to the address space.

    TimePt = (double *) shmat(TimeID, NULL, 0);

    // Check if Time segment is attached correctly

    if ( TimePt == (double *) -1) {
        printf("Time segment attach error");
        return 1;
    }

    // Get temperature from the shared memory

    double timeStep = TimePt[0];

    // Detach segments

    shmdt(InfoPt);
    shmdt(TimePt);

  return timeStep;
}
//========================================================================
//
bool DiffusionModel::checkResidual()
{
	double *residual = equations.getLoadVector();
	double maxResidual = 0.0;
	double val = 0.0;
	bool flag = false;

	//could possibly multithread this check?? worth the start up costs?
	for(int i=0;i<totalNumDof;i++)
	{
		val = fabs(residual[i]);
		if(val > maxResidual)
		{
			maxResidual = val;
		}
	}

	if(maxResidual > iterationTol){
		flag = false;}
	else {
		flag = true;}

	return flag;
}
//========================================================================
void DiffusionModel::checkIfNeedsRefactorization(int &localTimeStepCount, int iterationCount)
{

	if(localTimeStepCount >= 20){
		BETA_OUT<<"Allowable Local Time Step Count Exceeded - Refactorizing LHS"<<endl;
		needsRefactorize = true;
	}
	if(iterationCount >= 3){
		BETA_OUT<<"Allowable Iteration Count Exceeded - Refactorizing LHS"<<endl;
		needsRefactorize = true;
	}

	if(needsRefactorize)
	{
		equations.ReleaseMemory();
		equations.setUsePreviousFactorizedMatrixAtRuntime(false);
		localTimeStepCount = 0;
	}
	
}
//========================================================================
void DiffusionModel::DoEnhancedTransientNonLinearAnalysis(istream * inStream)
{
	//equations.setUsePreviousFillInOrdering(true);
	bool DoIterate=false;
	bool UseOldK=false;
	double maxAllowableEnhancedResidual=1e-10;

	checkModel();

	mainTimeCount=0.0;
	mainTimeStepCount=0;

	InitializeOptionalOutput();
	double * deltaSolution;

	BETA_OUT << "solution:" << endl;
	BETA_OUT << "\n Timestep# time solution " << endl;

	int numTimeSegments=TimeSegmentList.getNumElements();
	for (int iTimeSegment=0; iTimeSegment<numTimeSegments; iTimeSegment++) {//loop over timesegments
		int numTimeSteps=TimeSegmentList[iTimeSegment].getNumElements();
		for (int iCount=0; iCount<numTimeSteps; iCount++) {//loop over time steps
			timeStepSize=TimeSegmentList[iTimeSegment][iCount];

			AnalysisParameters[0] = time_integration_alpha*timeStepSize;
			AnalysisParameters[1] = (1-time_integration_alpha)*timeStepSize;
			
			mainTimeCount+=timeStepSize;
			mainTimeStepCount++;
			cout << "Current Time= " << mainTimeCount << " seconds." << endl;

			// everytime clear rooms in equation xt: 102200
			equations.zeroLoadVector();
			if(!UseOldK)
				equations.zeroMatrix();
			
			// Set up applied loadVector; Send loads to equation object
			for(int i=0; i<(int)loads.size(); i++){
	            if (loads[i]->isIncremental())
	                loads[i]->setRate(0.01);
				else
					loads[i]->setRate(1.0);  // uniform loading
				loads[i]->applyTo(equations,&mesh->node,mesh->numNodes, this);
			}

			// apply additional equations: important for dummy load !
			equations.applyAdditionalEquations(); 

			//enhanced solution algorithm
			//first obtain linear solution:
			sumVector(u_s1,		u_s,		delta_u,totalNumDof);
			if(UseOldK){
				assembleF_s_alpha_family(u_s, delta_u, u_s1, AnalysisParameters, timeStepSize); 
//				equations.setDoDecomposeFlag(false);
				deltaSolution=equations.solve(u_s1);
			}
			else{
				assembleK_F_s_alpha_family(u_s, delta_u, u_s1, AnalysisParameters, timeStepSize); 
				deltaSolution=equations.solve(u_s1);
				UseOldK=true;
			}

			copyVector(deltaSolution, delta_u, totalNumDof); //copy deltaU to local array
			sumVector(u_s1,		u_s,		delta_u,totalNumDof);

			//second step is to iterate to reduce the residual
			int iIter=0;
			int nIter = 2;  //fix this later
			//for calculating residual
			double maxResidual;
			int    dofForMaxResidual;
			double * mappedResidualVector=0;

			getEnhancedResidual(internalForce,
							AnalysisParameters, timeStepSize,
							u_s, delta_u, u_s1, 
							mappedResidualVector,maxResidual,dofForMaxResidual );

			while ( (maxResidual > maxAllowableEnhancedResidual) && (iIter <= nIter) ){

				deltaSolution=equations.incrementalSolutionUsingResidual(u_s1);

//				equations.addIncrementalDisplacements(deltaSolution);
				sumVector(delta_u,	delta_u,	deltaSolution,totalNumDof);
				sumVector(u_s1,		u_s,		delta_u,totalNumDof);

				mesh->printNodalData(u_s1, "displacements", &BETA_OUT );
				getEnhancedResidual(internalForce,
								AnalysisParameters, timeStepSize,
								u_s, delta_u, u_s1, 
								mappedResidualVector,maxResidual,dofForMaxResidual );
		    
				iIter++;
			}//end of iterate loop
			cout << "Interations taken: " << iIter << endl;
			if(iIter > nIter){ FatalError("Exceeded allowed number of iterations");}
			//end of second step of enhanced solution

			TimeStepUpdate();//update Phi

//			BETA_OUT << "\n" << iCount << " " << count << " " << FORMAT << u_s1[1];

			if(DoIterate) TransientAnalysisIterate();

			WriteOptionalOutput();

			//BETA_OUT << iCount << " " << count << " " << u_s1[1] << endl;

			initializeVector( delta_u, 0.0, totalNumDof);
			copyVector( u_s1, u_s, totalNumDof);
		}//end TimeStep loop
		BETA_OUT << endl;
	}//end TimeSegment loop
	BETA_OUT << endl;

	FinalizeOptionalOutput();
} // end of DoTransientNonLinearAnalysis
//========================================================================
void DiffusionModel::TransientAnalysisIterate()
{
	int iIter=0;
	int nIter = 50;  //fix this later
    double * deltaSolution;

	//for calculating residual
    double maxResidual;
    int    dofForMaxResidual;
    double * mappedResidualVector=0; // = new double[totalNumDof]; //do not modify mappedResidualVector
	//only for 'reading', because mappedResidualVector points to equation object's data

	getGlobalForcesAndResidual(internalForce,
					AnalysisParameters, timeStepSize,
					u_s, delta_u, u_s1, 
					mappedResidualVector,maxResidual,dofForMaxResidual );

//	BETA_OUT << " (" << mappedResidualVector[1] << ") ";

	while ( (maxResidual > maxAllowableResidual) && (iIter <= nIter) ){

		assembleK_alpha_family(u_s, delta_u, u_s1, AnalysisParameters, timeStepSize); 

		deltaSolution=equations.incrementalSolutionUsingResidual(u_s1);

		equations.addIncrementalDisplacements(deltaSolution);
		sumVector(delta_u,	delta_u,	deltaSolution,totalNumDof);
		sumVector(u_s1,		u_s,		delta_u,totalNumDof);

//				BETA_OUT << " " << u_s1[1];
		mesh->printNodalData(u_s1, "displacements", &BETA_OUT );

		getGlobalForcesAndResidual(internalForce,
						AnalysisParameters, timeStepSize,
						u_s, delta_u, u_s1, 
						mappedResidualVector,maxResidual,dofForMaxResidual );

//				BETA_OUT << " (" << mappedResidualVector[1] << ") ";
    
		iIter++;
	}//end of iterate loop

	cout << "Interations taken: " << iIter << endl;
	if(iIter > nIter){ FatalError("Exceeded allowed number of iterations");}
}
//========================================================================
void DiffusionModel::getGlobalForcesAndResidual( double *internalForce,
							    Array<double> &AnalysisParameters, double timeStep,
							    double *u_s, double *delta_u, double *u_s1, 
								double *mappedResidualVector, double &maxResidual, int &dofForMaxResidual)
{
	double a1=AnalysisParameters[0];
	double a2=AnalysisParameters[1];
	
	mesh->element[0].initializeSummary();

	int i;
	double * elementForce;

	for(i=0;i<totalNumDof;i++) internalForce[i]=0;
	equations.zeroResultantVector();	
	BasicElement *e=0;
	for(i=0; i<mesh->numElements; i++)
	{
		e=&mesh->element[i];
		e->initializeArrays();
		e->update("FM", u_s1); // calculate M as well
		((HeatTransferElement3Dtransient *) e)->formResultantVector(AnalysisParameters, u_s, delta_u, u_s1, timeStep); //latest solution
		if(a2 != 0.0){
			e->update("F", u_s); 
			((HeatTransferElement3Dtransient *) e)->addA2TermToResultantVector(AnalysisParameters, u_s, delta_u, u_s1, timeStep); //previous solution
		}

		int numDof   = e->getNumDof();
		int *dofList = e->calculateDofList();
		elementForce = e->getForceVector();

		equations.assembleResultantVector(dofList,numDof,
			                              elementForce);
		for(int ii=0; ii<numDof; ii++)// InternalForce[] is not mapped
		{ 
			internalForce[dofList[ii]] += elementForce[ii];	
		}
	}//end of loop on elements

	double residualNorm = equations.calculateResidualNorm();
	//BETA_OUT<<"Residual norm= "<<residualNorm
	//    <<endl;
	double *ResidualVectorInEquations = equations.getResidualVector(); // mapped
	for(i=0;i<totalNumDof;i++) {
		mappedResidualVector[i]=ResidualVectorInEquations[i];
	}
	maxResidual       = 0.;
	dofForMaxResidual = -1;

	for(i=0;i<totalNumDof;i++) {
		if( fabs(mappedResidualVector[i] )> maxResidual ){
			maxResidual =  fabs(mappedResidualVector[i]);
			dofForMaxResidual=i; 
		}
	}
}//end of getGlobalForces
//===============================================================
void DiffusionModel::getEnhancedResidual( double *internalForce,
							    Array<double> &AnalysisParameters, double timeStep,
							    double *u_s, double *delta_u, double *u_s1, 
								double *mappedResidualVector, double &maxResidual, int &dofForMaxResidual)
{
	mesh->element[0].initializeSummary();

	int i;
	double * elementForce;

	for(i=0;i<totalNumDof;i++) internalForce[i]=0;
	equations.zeroResultantVector();	
	BasicElement *e=0;
	for(i=0; i<mesh->numElements; i++)
	{
		e=&mesh->element[i];
		int group = e->getMaterialNumber();		

		e->initializeArrays();
		e->update("N", delta_u); // calculates Mbar times delta_q... 
							  //this includes M*delta_q and K*delta_q
		((HeatTransferElement1Dtransient *) e)->formEnhancedResultantVector(AnalysisParameters, u_s, delta_u, u_s1, timeStep);

		int numDof   = e->getNumDof();
		int *dofList = e->calculateDofList();
		elementForce = e->getForceVector();

		equations.assembleResultantVector(dofList,numDof,
			                              elementForce);
		for(int ii=0; ii<numDof; ii++)// InternalForce[] is not mapped
		{ 
			internalForce[dofList[ii]] += elementForce[ii];	
		}
	}//end of loop on elements

	double residualNorm = equations.calculateResidualNorm();
	//BETA_OUT<<"Residual norm= "<<residualNorm
	//    <<endl;
	double *ResidualVectorInEquations = equations.getResidualVector(); // mapped
	for(i=0;i<totalNumDof;i++) {
		mappedResidualVector[i]=ResidualVectorInEquations[i];
	}
	maxResidual       = 0.;
	dofForMaxResidual = -1;

	for(i=0;i<totalNumDof;i++) {
		if( fabs(mappedResidualVector[i] )> maxResidual )
			{ maxResidual =  fabs(mappedResidualVector[i]);
			dofForMaxResidual=i; }
	}
}//end of getEnhancedResidual
//===============================================================
void DiffusionModel::getGlobalFlowVector( double *internalForce,
							    Array<double> &AnalysisParameters, double timeStep,
							    double *u_s, double *delta_u, double *u_s1)
{
	double a1=AnalysisParameters[0];
	double a2=AnalysisParameters[1];

	mesh->element[0].initializeSummary();

	int i;
	double * elementForce;

	for(i=0;i<totalNumDof;i++) internalForce[i]=0;
	equations.zeroResultantVector();	
	HeatTransferElement3Dtransient *e=0;
	for(i=0; i<mesh->numElements; i++)
	{
		e=(HeatTransferElement3Dtransient *)&mesh->element[i];
		e->initializeArrays();
		e->update("FM", u_s1); // calculate M as well
		e->formFlowVector(AnalysisParameters, u_s, delta_u, u_s1, timeStep); //latest solution
		if(a2 != 0.0){
			e->update("F", u_s); 
			e->addA2TermToResultantVector(AnalysisParameters, u_s, delta_u, u_s1, timeStep); //previous solution
		}

		int numDof   = e->getNumDof();
		int *dofList = e->calculateDofList();
		elementForce = e->getForceVector();

		equations.assembleResultantVector(dofList,numDof,
			                              elementForce);
		for(int ii=0; ii<numDof; ii++)// InternalForce[] is not mapped
		{ 
			internalForce[dofList[ii]] += elementForce[ii];	
		}
	}//end of loop on elements
}//end of getGlobalFlowVector
//===============================================================

/*
//======================================================================
void DamageModel::DoDMAnalysis(istream * inStream)
{
	time_t startTime, endTime;
    
	int i, iele, ele, matID;
	char fm[7];

    double maxResidual;
    int    dofForMaxResidual;
    double * mappedResidualVector = new double[totalNumDof];
    
	double aRateInit;

    double *currentSolutionInEquations;
    double * deltaSolution;

    int iIter;
    int count=0;
	int k=0;
	int kk=0;
	int kkk=0;
	int kkkk=0;
	int kkkkk=0;

    ofstream  os, os1;
    char cdump[500];

	double aRate=1.0/double(nStep);  // this is a trial load step...to be scaled later...

    FormWriter fr(BETA_OUT);

	if(allocationState != 1){allocateForStaticAnalysis(elementWorkspace);}
    
    // initializing solution vectors
    solution = new double[totalNumDof];

    equations.zeroLoadVector();
    equations.zeroMatrix(); //xtang 990525

    // "Load loop" begins here

    // now goto nonlinear routine by using residual

	// xtang 990720
	DMElasticElement *e=(DMElasticElement*)(&mesh->element[0]);
	DMElasticElement::openOutputFiles(filemanager, OutputOptions);

	DMElasticElement::eVolDist->seekp(0,ios::beg);
	DMElasticElement::eDamage->seekp(0,ios::beg);


	// xtang 991016: add searching scaling factor for the trial load step
	// (1) Find scaling factor for the 1st load step

	// XTANG 991018: Output the first line of OUT for iteration info
    //OUT <<count<<" 0.0   0  0.0   0.0   0.0 "<<endl;
	BETA_OUT <<count<<"0.0"<<"\t"<<"0"<<"\t"<<"0.0"<<"\t"<<"0.0"<<"\t"<<"0.0"<<endl;

//testing
	ofstream solfile;
//end testing

	checkModel();

    do {    
		cout << "********\nload step number :" << count <<"\n*************\n" << endl;
		//displayTime("LoadStep ");

        // Set up applied loadVector; 
        // Send loads to equation object
        equations.zeroLoadVector();

		for( i=0; i<(int)loads.size(); i++){
            if(loads[i] == NULL) break;
            if (loads[i]->isIncremental()) { 
                loads[i]->setRate(aRate);
            }
            loads[i]->applyTo(equations,&mesh->node,mesh->numNodes);
		}

      
        // apply additional equations: important for dummy load !
        equations.applyAdditionalEquations();

        // switch Primary BC
        currentSolutionInEquations = equations.getSolution();
        copyVector(currentSolutionInEquations,currentSolution,totalNumDof);  // U_r+1 after soving

        iIter=0;
        do {  // iteration begins here
	

            assembleKandF(currentSolution); 
            getGlobalForces(internalForce,currentSolution, 
                mappedResidualVector,maxResidual,dofForMaxResidual );
            
            if(maxResidual< maxAllowableResidual) { break;}

    startTime = time(NULL);

            deltaSolution=equations.incrementalSolutionUsingResidual(currentSolution);
            
    endTime = time(NULL);
            
            equations.addIncrementalDisplacements(deltaSolution);
            
            copyVector(currentSolution,solution,totalNumDof); // U_r after solving

            currentSolutionInEquations = equations.getSolution();
            copyVector(currentSolutionInEquations,currentSolution,totalNumDof);  // U_r+1 after soving
     
			if (count>0) {
    			DMElasticElement::eQuadInfo->seekp(0,ios::beg);
    
				for(int i=0;i<numElements;i++) { 
					mesh->element[i].update("VQ",currentSolution);
				}

                sprintf(cdump,"%dvA.dat",count);
                DMElasticElement::eQuadInfo->seekg(0,ios::beg);
                outVolumeDistribution(7,cdump,count,iIter);  // check ISV for this number
				iIter++;
			}
			else if (count==0) {
                DMElasticElement::maxStrengthRatio=0.0;
			    for(i=0;i<numElements;i++) { 
                    mesh->element[i].update("R",currentSolution); // R: check strength value only
				}
				// we don't have to warry about the maxStrengthRatio=1.0 for now 
		        {aRate/=DMElasticElement::maxStrengthRatio;}

								
				aRateInit=aRate;
				for(i=0; i<totalNumDof; i++) {currentSolution[i]*=(aRateInit*double(nStep));}
				break;
			}

            // epsilon=differenceAB(currentSolution, solution,totalNumDof);

        } while ((iIter <= nIter) );  //end of iterate loop

		count++;

		// xtang 991018: output info for current load step
        BETA_OUT <<count<<"\t"<<aRate <<"\t";
        BETA_OUT<<iIter<<"\t"<<maxResidual<<"\t";
		if (DMElasticElement::jamming) { BETA_OUT << "JM\t"; }
		else BETA_OUT <<"NJ\t";
		
		
		//+DG_Oct11_2006
		if(dispRequired)
		{
			if(count==1 || count==stepIncrements*k)
			{
				char cdumpDG[500];
				char outBase[] = {'d','i','s','p','\0'};
				ofstream *osdg;
				sprintf(cdumpDG,"%d.%s",count,outBase);  // like "1.disp"
				osdg = filemanager->OpenOutputStream(cdumpDG);
				mesh->printNodalData(currentSolution, "Nodal Displacements", osdg );  
				osdg->close();
				if(count==1 && stepIncrements==1){k=k+1;}
				k++;
			}
		}

    } while ((count <= nStep) && (aRate <= 1.0));
                
    delete [] mappedResidualVector;
    delete [] solution;

} // end of DoDMAnalysis
*/
//======================================================================
void DiffusionModel::WriteOptionalOutput()
{
	ofstream *ost;
	char aName[200];

	bool lastTimestep=false;
	if (totalTimeSteps==mainTimeStepCount) lastTimestep=true;

	if ((concentration_out) && 
		( (mainTimeStepCount % concentration_out_Step_Interval == 0) || lastTimestep) ) {
			sprintf(aName,"concentration%d.disp",mainTimeStepCount);
			ost = filemanager->OpenOutputStream(aName);

			double outTime = 0.0;
			double *uPrint = NULL;
			if(usingTimeStepControl){
				outTime = printTimeCount;
				uPrint = ut_s;
			}
			else{
				outTime = mainTimeCount;
				uPrint = u_s1;
			}

			if(verboseFlag > Basic)
				BETA_OUT << "Writing optional output- concentration timestepcount: " << mainTimeStepCount 
					<< ", time: " << outTime << ", filename: " << aName << endl;
			mesh->printNodalData(uPrint, "displacements", ost);
			filemanager->CloseOutputStream(ost);
			*ConcentrationAnimationInput << mainTimeStepCount << "\t" << outTime << "\t" << aName << endl;

			sprintf(aName,"VolAvg%d.txt",mainTimeStepCount);
			char currentDirectory[_MAX_PATH];
			_getcwd( currentDirectory, _MAX_PATH );
			*PlotModelAvgConcentrationInput << currentDirectory <<"\\" << aName<< "\t" << mainTimeCount << "\t" << "1.0" << endl;
	}	
	if ((flux_out) && 
		((mainTimeStepCount % flux_out_Step_Interval == 0) || lastTimestep) ) {
			sprintf(aName,"flux%d.disp",mainTimeStepCount);
			ost = filemanager->OpenOutputStream(aName);
			if(verboseFlag > Basic)
				BETA_OUT << "Writing optional output- flux timestepcount: " << mainTimeStepCount 
					<< ", time: " << mainTimeCount << ", filename: " << aName << endl;
			PrintFluxes(*ost, u_s1);
			filemanager->CloseOutputStream(ost);
			*FluxAnimationInput << mainTimeStepCount << "\t" << mainTimeCount << "\t" << aName << endl;
	}
		if ((quad_flux_out) && 
		((mainTimeStepCount % quad_flux_out_Step_Interval == 0) || lastTimestep) ) {
			sprintf(aName,"flux%d.quad",mainTimeStepCount);
			ost = filemanager->OpenOutputStream(aName);
			if(verboseFlag > Basic)
				BETA_OUT << "Writing optional output- flux timestepcount: " << mainTimeStepCount 
					<< ", time: " << mainTimeCount << ", filename: " << aName << endl;
			PrintQuadFluxes(*ost, u_s1);
			filemanager->CloseOutputStream(ost);
			*QuadFluxAnimationInput << mainTimeStepCount << "\t" << mainTimeCount << "\t" << aName << endl;
	}	


	if ((actualConcentration_out) && 
		((mainTimeStepCount % actualConcentration_out_Step_Interval == 0) || lastTimestep) ) {
			sprintf(aName,"ActualConc%d.disp",mainTimeStepCount);
			ost = filemanager->OpenOutputStream(aName);
			if(verboseFlag > Basic)
				BETA_OUT << "Writing optional output- flux timestepcount: " << mainTimeStepCount 
					<< ", time: " << mainTimeCount << ", filename: " << aName << endl;
			PrintActualConcentrations(*ost, u_s1);
			filemanager->CloseOutputStream(ost);
			*ActualConcentrationAnimationInput << mainTimeStepCount << "\t" << mainTimeCount << "\t" << aName << endl;
	}	

	if ((volavg_out) && 
		((mainTimeStepCount % volavg_out_Step_Interval == 0) || lastTimestep) ) {
		sprintf(aName,"VolAvg%d.txt",mainTimeStepCount);
		ofstream *volumeAverageout = filemanager->OpenOutputStream(aName);
		if(verboseFlag > Basic)
			BETA_OUT << "Printing Volume Averages..." << endl;
		mesh->element[0].initializeSummary(); // since getglobal forces is not called, i have to initialize the summary here 
		outVolumeStrConstituents(u_s1, 0, 0, 0, volumeAverageout, true);
		mesh->element[0].outputSummary(); //no sense in printing summary since it just got initialized
		filemanager->CloseOutputStream(volumeAverageout);
	}

	if ((GlobalForces_out) && 
		((mainTimeStepCount % GlobalForces_out_Step_Interval == 0) || lastTimestep) ) {
		sprintf(aName,"GlobalForces%d.txt",mainTimeStepCount);
		ost = filemanager->OpenOutputStream(aName);
		if(verboseFlag > Basic)
			BETA_OUT << "Printing Global Forces..." << endl;
		int    dofForMaxResidual;
		double maxResidual;
		double * mappedResidualVector = new double[totalNumDof];
		getGlobalForces(internalForce, u_s1, mappedResidualVector,
		                maxResidual,dofForMaxResidual );
		ConvertNumericalZeroToZero(internalForce, totalNumDof, 1e-12);
		mesh->printNodalData(internalForce, "displacements", ost);
		delete [] mappedResidualVector;
		
		filemanager->CloseOutputStream(ost);
	}
	if(QuadPointsRequired && mainTimeStepCount == 1){
        sprintf(aName,"quadPoints",mainTimeStepCount);
		ost = filemanager->OpenOutputStream(aName);
		if(verboseFlag > Basic)
		BETA_OUT << "Writing optional output - quadPoints: " << mainTimeStepCount 
		         << ", filename: " << aName << endl;
		PrintQuadPoints(*ost);
		filemanager->CloseOutputStream(ost);
			
    }



}
//==========================================================
void DiffusionModel::InitializeOptionalOutput()
{
	char aName[200];

	if (concentration_out)  {
			sprintf(aName,"ConcentrationAnimationInput.txt");
			ConcentrationAnimationInput = filemanager->OpenOutputStream(aName);
			*ConcentrationAnimationInput << "mainTimeStepCount \t mainTimeCount \t FileName " << endl;

			sprintf(aName,"PlotModelAvgConcWithTime.txt");
			PlotModelAvgConcentrationInput = filemanager->OpenOutputStream(aName);
			*PlotModelAvgConcentrationInput << "1 6\n0 0" << endl;
	}	
	if (flux_out)  {
			sprintf(aName,"FluxAnimationInput.txt");
			FluxAnimationInput = filemanager->OpenOutputStream(aName);
			*FluxAnimationInput << "mainTimeStepCount \t mainTimeCount \t FileName " << endl;
	}
	if (quad_flux_out)  {
			sprintf(aName,"QuadFluxAnimationInput.txt");
			QuadFluxAnimationInput = filemanager->OpenOutputStream(aName);
			*QuadFluxAnimationInput << "mainTimeStepCount \t mainTimeCount \t FileName " << endl;
	}	
	if (actualConcentration_out)  {
			sprintf(aName,"ActualConcentrationAnimationInput.txt");
			ActualConcentrationAnimationInput = filemanager->OpenOutputStream(aName);
			*ActualConcentrationAnimationInput << "mainTimeStepCount \t mainTimeCount \t FileName " << endl;
	}	
}
//==========================================================
void DiffusionModel::FinalizeOptionalOutput()
{
	if (concentration_out)  {
			filemanager->CloseOutputStream((ofstream*)ConcentrationAnimationInput);
			ConcentrationAnimationInput=0;
			filemanager->CloseOutputStream((ofstream*)PlotModelAvgConcentrationInput);
			PlotModelAvgConcentrationInput=0;
	}	
	if (flux_out)  {
			filemanager->CloseOutputStream((ofstream*)FluxAnimationInput);
			FluxAnimationInput=0;
	}
	if (quad_flux_out)  {
			filemanager->CloseOutputStream((ofstream*)QuadFluxAnimationInput);
			QuadFluxAnimationInput=0;
	}
	if (actualConcentration_out)  {
			filemanager->CloseOutputStream((ofstream*)ActualConcentrationAnimationInput);
			ActualConcentrationAnimationInput=0;
	}	
}
//======================================================================
int DiffusionModel::processOptionalOutput(char** localTokenList,ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;

	char token[200];
	int  foundMatch=0;

strcpy(token,localTokenList[0]);
BETA_OUT<<"\nCommand processor = 'DiffusionModel::processOptionalOutput'          "<<token<<endl;

//========
	if(COMPARE(token,"concentration")==0){
		concentration_out=true;
		concentration_out_Step_Interval=1;
		if(numToken==2)
			concentration_out_Step_Interval=atoi(localTokenList[1]);
		OK 
	}			

	if(COMPARE(token,"flux")==0){
		flux_out=true;
		flux_out_Step_Interval=1;
		if(numToken==2)
			flux_out_Step_Interval=atoi(localTokenList[1]);
		OK 
	}

	if(COMPARE(token,"quadflux")==0){
		quad_flux_out=true;
		quad_flux_out_Step_Interval=1;
		if(numToken==2)
			quad_flux_out_Step_Interval=atoi(localTokenList[1]);
		OK 
	}	

	if(COMPARE(token,"ActualConcentration")==0){
		actualConcentration_out=true;
		actualConcentration_out_Step_Interval=1;
		if(numToken==2)
			actualConcentration_out_Step_Interval=atoi(localTokenList[1]);
		OK 
	}			

	if(COMPARE(token,"volumeaverage")==0){
		volavg_out=true;
		volavg_out_Step_Interval=1;
		if(numToken==2)
			volavg_out_Step_Interval=atoi(localTokenList[1]);
		OK 
	}			

	if(COMPARE(token,"GlobalForces")==0){
		GlobalForces_out=true;
		GlobalForces_out_Step_Interval=1;
		if(numToken==2)
			GlobalForces_out_Step_Interval=atoi(localTokenList[1]);
		OK 
	}			

	if(COMPARE(token,"RunFakeAnalysisToGenerateOutputFiles")==0){
		RunFakeAnalysisToGenerateOutputFiles=true;
		OK 
	}			
	
	if(foundMatch ==0){	
		foundMatch=ElasticityModel::processOptionalOutput(localTokenList,currentStream, numToken);
	}

	return foundMatch;
}
//======================================================================
void DiffusionModel::DoAnalysis(istream * inStream)
{
	int i;

	FormWriter fr(BETA_OUT);

	allocateForAnalysis();

	u = equations.getSolution();  // unmapped
    u_s = new double[totalNumDof];//solution from previous timestep
	u_s1= new double[totalNumDof]; 
	for (i=0;i<totalNumDof;i++) {u_s[i]=u[i]=u_s1[i]=0.0;}

	// Send loads to equation object
	equations.zeroLoadVector();
	
	for( i=0; i<(int)loads.size(); i++){
		loads[i]->applyTo(equations,&mesh->node,mesh->numNodes, this);
	}

	checkModel();

	InitializeOptionalOutput();
		
	BETA_OUT<<"Assemble stiffness matrix"<<endl;
	displayTime("", BETA_OUT);
	
    assembler->assembleInitialF();
	assembler->assembleK(NULL);
//	assembleK(currentSolution);

	//equations.printMatrix(fr);
	//equations.printLoadVector(fr);

	displayTime("Time to assembleKandF", BETA_OUT);
	BETA_OUT<<"Start solving equations (" << totalNumDof << ") :"<<endl;
	displayTime("", BETA_OUT);

	cout << "solving" << endl;
	solution = equations.solve();
	displayTime("Time to solve", BETA_OUT);

	for(i=0; i<totalNumDof; i++){u_s1[i] += solution[i];}
		
	ofstream *os=filemanager->OpenOutputStream("disp");
	mesh->printNodalData(u_s1, "displacements", os );
	filemanager->CloseOutputStream(os);
	
	haveSolution = true;

	WriteOptionalOutput();
	FinalizeOptionalOutput();

	if( BasicElement::getAnalysisType() == 1) iterate();

    delete [] u_s; u_s=0;
	delete [] u_s1; u_s1=0;
} // end of DoAnalysis
//=========================================================================
void DiffusionModel::outVolumeStrConstituents(double *disp,int matVol, int numCpt, double ** cpt, ostream *volumeAverageout, bool byConstituent=false)
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
	HeatTransferElement3Dtransient *e=0;

	for(i=0;i<100;i++) {matID[i]=0;}

	if (byConstituent) {
		for(i=0;i<numElements;i++) {
			e=(HeatTransferElement3Dtransient *)&mesh->element[i];

			is=e->getMaterialNumber();
			jj=getIndex(is,matID,numMats+1);
			if (jj < 0) {numMats++; matID[numMats]=is; }
		}
	}

	int maxStresses=((HeatTransferElement3Dtransient *)&mesh->element[0])->getNumStrains();
    volAvgStress = new double [maxStresses];
    volAvgStrain = new double [maxStresses];
    SED = new double [3];

	for(jj=0;jj<maxStresses;jj++){
		volAvgStress[jj]=0.0;
		volAvgStrain[jj]=0.0;
	}

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
		e=(HeatTransferElement3Dtransient *)&mesh->element[i];
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
        volume=0.0;
		SED[0]=0.0;
        ((HeatTransferElement3Dtransient*)e)->getVolAvgValues(disp,volAvgStress,volAvgStrain, SED, &volume);  
		
	    Totalvolume[0]   += volume;
		TotalSED[0] += SED[0];
		//TotalSED is actually storing the total concentration

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

	HeatTransferMaterial * matPoint;
	e=(HeatTransferElement3Dtransient*)(&mesh->element[0]); // Fix this... what should be done ??

	for (i=0; i<=numMats; i++) 
	{
		*(volumeAverageout) <<matID[i]<<"\t"<<Totalvolume[i]<<"\t";  
		matPoint = (HeatTransferMaterial * ) &(materialList[matID[i]]);				
		if (matPoint) 
		{
//			*(volumeAverageout) << matPoint->xtNumAngles<<"\t"<<matPoint->xtAngles[0]<<"\t"<<matPoint->xtAxis[0];
			*(volumeAverageout) << -1<<"\t"<<0<<"\t"<<0;
		}
		else 
			*(volumeAverageout) << -1<<"\t"<<0<<"\t"<<0;
		*(volumeAverageout)<<"\t"<<TotalSED[i]/Totalvolume[i]; // volume averaged concentration
		*(volumeAverageout)<<endl; 
		for(is=0; is<maxStresses; is++)  
			*(volumeAverageout)<<DOUBLE_FORMAT<<TotalvolAvgStrain[i][is]/Totalvolume[i]<< ' ';
	    *(volumeAverageout)<<endl;
		for(is=0; is<maxStresses; is++)  
			*(volumeAverageout)<<DOUBLE_FORMAT<<TotalvolAvgStress[i][is]/Totalvolume[i]<< ' ';
	    *(volumeAverageout)<<endl;
	}
    *(volumeAverageout)<<endl;

//creating input file for ross' matlab code:
	ofstream *ost;
	char aName[200];
	sprintf(aName,"VolAvgPlotInput%d.txt",mainTimeStepCount);
	ost = filemanager->OpenOutputStream(aName);
	*ost << "group# volume avg.conc avg.gradient [1,2...] avg.flux [1,2...]" << endl;
	for (i=1; i<=numMats; i++) { //ignoring '0' material which is the complete mesh region
		*ost << matID[i] << '\t' << DOUBLE_FORMAT << Totalvolume[i] << '\t' << TotalSED[i]/Totalvolume[i] << '\t';
		for(is=0; is<maxStresses; is++)  
			*(ost)<<DOUBLE_FORMAT<<TotalvolAvgStrain[i][is]/Totalvolume[i] << '\t';
		for(is=0; is<maxStresses; is++)  
			*(ost)<<DOUBLE_FORMAT<<TotalvolAvgStress[i][is]/Totalvolume[i] << '\t';
		*(ost)<<endl;
	}
	filemanager->CloseOutputStream(ost);
//end of creating input file for ross' matlab code

	//JV081007
	//updating the element's vol avg data variable. 
	//so that element[0].outputSummary() prints the volume averaging info in 
	//the output file as well. this is just a convenience fix.
	//should decide a more formal procedure to do this later.
	ElementWorkspace* bag=e->getElementWorkspace();
	bag->totalVolume = Totalvolume[0];
	for(i=0; i<maxStresses; i++) {
		bag->strainVolume_total[i] = TotalvolAvgStrain[0][i];
		bag->stressVolume_total[i] = TotalvolAvgStress[0][i];
		bag->localStressVolume_total[i] =0.;  
		bag->energyVolume_total[i] =0; 
	}


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
//#undef maxStresses
//#undef FORMAT
}
//======================================================================
bool DiffusionModel::checkModel()
{
	ElasticityModel::checkModel(); //calling parent's checkModel first

	if(allocationState==0) EXIT_BETA("AllocationState flag not set.\nThis probably means 'ReadTransientAnalysisParameters' has not been called.\nFix the model script.");

	return true;
}
//==========================================================================
