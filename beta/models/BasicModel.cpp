#include "stdafx.h"

#include "BasicModel.hpp"
#include "factory/Factory.hpp"
#include "utility/formWriter.hpp"
#include "mesh/NodeSet.hpp"

#include "elements/IntegrationDataType.hpp"

#include "assemblers/SerialAssembler.hpp"
#include "assemblers/ESA_Assembler.hpp"
#include "assemblers/OFA_Assembler.hpp"

#include "elements/ExtrapolationMatrix.hpp"

//-------------------------------------------------------------------
#include <omp.h>
#include <numeric>
//-------------------------------------------------------------------
extern Factory	*factory;
extern int verboseFlag;
//-------------------------------------------------------------------

void setVectorFromMeshDisplacements(BasicMesh *mesh, double *vector);
ElementGroup* getElementGroup(char* name, Array<ElementGroup> *GroupList);
ElementGroup* CreateElementGroupFromElematFile(string elematfilename, BasicMesh *mesh, FileManager *filemanager);
void setVectorFromNodalValueFile(char* filename, double *vector, int nn, NodeGroup &node);


//nonlinear relaxation
#define None        0
#define Crisfield 100
#define Felippa   101

extern double findNodeTOLERANCE;

//-------------------------------------------------------------------
CreateErrorHandler(BasicModel);
//--------------------------------------------------------------------
BasicModel::BasicModel(void)
{
	FileManagerType = COMMON_OUTPUT_SCHEME;
	initializeFlags();
	filemanager			= 0;
	maxAllowableResidual= 1.e-10 ;
	solution			= 0;
	allocationState		= 0;
	
//	elementWorkspace=0;

//	Fe=0;
//	initialFe=0;
//	Ke=0;

	totalNumDof		= 0;
	appliedForces	= 0;
	currentSolution	= 0;
	internalForce	= 0;
    zeroVector		= 0;

	mesh			= 0;
	assembler = NULL;
	materialList.reserve(10);

	initializeRestartFlags();
}
//--------------------------------------------------------------------
BasicModel::~BasicModel(void)
{
	if(mesh) delete mesh;
	if(assembler) delete assembler;
//	if(elementWorkspace) delete elementWorkspace;

	int size;
	size=(int)constraints.size();
	if( !constraints.empty() ){
		for(int i=0;i<size;i++)
			delete  constraints[i];
		constraints.clear();
	}

	size=(int)loads.size();
	if( !loads.empty() ){
		for(int i=0;i<size;i++)
			delete  loads[i];
		loads.clear();
	}

	if(appliedForces)	delete [] appliedForces;
	if(currentSolution) delete [] currentSolution;
	if(internalForce)	delete [] internalForce;
	if(zeroVector)		delete [] zeroVector;
	if(filemanager) {
		delete filemanager;
		filemanager=0;
	}
}
//--------------------------------------------------------------------
void BasicModel::initialize()
{
	mesh			= new BasicMesh;
}
//--------------------------------------------------------------------
void BasicModel::initializeFlags()
{
	haveSetElementType  = haveSetElementMaterials = false ;
	haveReadMesh        = haveReadMaterials       = haveReadLoads = false;
	haveSetNumDofPerNode= haveSolution = false;
	StoreAllMe = false;
	

	UseOpenMP=false;
	BETA_num_threads=1;

	parallelAssembly_OFA=false;
	parallelAssembly_ESA=false;
	ParallelAssembly=false;
	BETA_num_threads_for_assembly=1;
}
//--------------------------------------------------------------------
void BasicModel::setFileManagerType(char * option)
{
	if(COMPARE(option, "COMMON_OUTPUT_SCHEME") ==0 )
		FileManagerType = COMMON_OUTPUT_SCHEME;
	else if(COMPARE(option, "QUALIFIED_OUTPUT_SCHEME") ==0 )
		FileManagerType = QUALIFIED_OUTPUT_SCHEME;
	else{
		EXIT_BETA("Cannot recognize output scheme (filemanager type)");
	}
}
//--------------------------------------------------------------------
bool BasicModel::checkModel()
{
	int numMats=materialList.getNumElements();
	//checking if material properties have been assigned to each element
	for(int i=0;i<mesh->numElements;i++){
		int matnum=mesh->element[i].getMaterialNumber();
		if(matnum < 0)
			EXIT_BETA("Material property not assigned for element number:" << i);
		if(matnum >= numMats )
			EXIT_BETA("Material #" << matnum << " assigned for element number:" << i << " does not exist. check material library.");
		if(COMPARE(materialList[matnum].getGroupType(),"")==0)
			EXIT_BETA("Material #" << matnum << " assigned for element number:" << i << " does not exist (or has not been constructed properly). check material library.");
	}
	return true;
}
//--------------------------------------------------------------------
bool BasicModel::readCommands(ifstream *defaultStream)
{
	char *localTokenList[20];
	int  numberOfTokens;
	int foundMatch;

	ifstream *currentStream;
	    BETA_OUT<<"Current stream in BasicModel::readCommands ="<<(*defaultStream)<<endl;
	currentStream = defaultStream;

	while(2==2) {
		foundMatch=0;		
		currentStream = defaultStream;		
		if( getLineAndTokenize(currentStream,"end",localTokenList, numberOfTokens)==1) {
			goto _ExitChoose;
		}
	
		foundMatch=processCommand(localTokenList,currentStream, numberOfTokens); // decendent class first

		if(foundMatch ==0) {	exit(1);}
		
	} //end of while

	_ExitChoose:
	BETA_OUT << "Exiting reader"<<endl; 
	defaultStream->close();
	return true;
}
//===================================================================================
int BasicModel::processCommand(char** localTokenList,ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;
	char newFileName[80];

	char token[200];
	char inputString[300];
	int  foundMatch=0;

strcpy(token,localTokenList[0]);
BETA_OUT<<"\nCommand processor = 'BasicModel::processCommand'          "<<token<<endl;

	if( COMPARE(token,"readFromCommandFile")==0 ){
		strcpy(newFileName,localTokenList[1]);
		BETA_OUT<<"Read commands in file:  "<<newFileName<<endl;
		ifstream  newStream(newFileName);
		readCommands(&newStream);  
		OK
		return foundMatch;
	}

	if( COMPARE(token,"systemCall")==0 ){
		BETA_OUT<<"systemCall:  "<<endl;
		(*currentStream).getline(inputString, 300,'\n');
        system(inputString);

		OK
		return foundMatch;
	}

	if( COMPARE(token,"openFile")==0 ) {   //This is for this pass only
		//Should be 3 tokens
		strcpy(newFileName,localTokenList[1]);
		currentStream = openFileForInput(newFileName);
		strcpy(token,localTokenList[2]);
		BETA_OUT<<"File for new stream = "<<newFileName<<endl;
		BETA_OUT<<"New stream  = "<<currentStream<<endl;
	}

	//========
	if( COMPARE(token,"CreateElements")==0 ) { // could detour to user-defined routine
		banner("CreateElements", BETA_OUT);
		factory->ElementFactory(currentStream, &mesh->element,mesh->numElements, filemanager->outStream);
		haveSetElementType = true; 
		OK
	}

	if( COMPARE(token,"MapQuadPointsToMesh")==0) {
		char data[80];
		char mapFileName[80];
		int numEl;
		int numQPPerEl;
		strcpy(data,localTokenList[1]);
		numEl = atoi(data);
		strcpy(data,localTokenList[2]);
		numQPPerEl = atoi(data);
		strcpy(newFileName,localTokenList[3]);
		strcpy(mapFileName,localTokenList[4]);
		cout<<"MapQuadPointToMesh"<<endl;
		MapQuadPointsToMesh(numEl,numQPPerEl,newFileName,mapFileName);
		OK
	}

		if( COMPARE(token,"MapQuadFile")==0) {
		char data[80];
		char mapFileName[80];
		char originQuadFileName[80];
		char destQuadFileName[80];
		int numEl;
		strcpy(data,localTokenList[1]);
		numEl = atoi(data);
		strcpy(mapFileName,localTokenList[2]);
		strcpy(originQuadFileName,localTokenList[3]);
		strcpy(destQuadFileName,localTokenList[4]);
		cout<<"MapQuadFile"<<endl;
		MapQuadFile(numEl,mapFileName,originQuadFileName,destQuadFileName);
		OK
	}

	IS_COMMAND(DoAnalysis)		
	IS_COMMAND(DoPostProcess) //bco 12/19/08
	IS_COMMAND(OptionalOutput)
	IS_COMMAND(ReadConstraints)
	IS_COMMAND(ReadLoads)
	IS_COMMAND(ReadMesh)
	IS_COMMAND(ReadOptionalOutput)  
    IS_COMMAND(ReadElementOrientations)
	IS_COMMAND(CreateNodeGroup)
	IS_COMMAND(CreateElementGroup)
	IS_COMMAND_2(BasicElement::setAnalysisType)
	IS_COMMAND(SetElementProperty)
	IS_COMMAND(SetNumDofPerNode)

	if( COMPARE(token,"readSpecialElementCommand")==0 ){//this needs to be fixed... maybe add another command to handle passing commands to particular elements
		char command[80];
		(*currentStream)>>command;
		BETA_OUT<<"Command= "<<command<<endl;

		if(numToken==2){ // then the command is for a specific element group, otherwise the first element in the mesh handles the command
			ElementGroup *elist=getElementGroup(localTokenList[1], &mesh->ElementGroupsList);
			if(elist==0) FatalError("readSpecialElementCommand: ElementGroup not found!\n");
			(*elist)[0].readSpecialCommand(*currentStream,elist,elist->getNumElements(), command);

		}else{
		mesh->element[0].readSpecialCommand(*currentStream,&mesh->element,mesh->numElements, command);
		}
		OK 
	} 

	if( COMPARE(token,"maxResidual")==0 ) {  // Move to nonlinear class
		(*currentStream)>>maxAllowableResidual;
		BETA_OUT<<"Set max allowable residual = "<< maxAllowableResidual<<endl;OK
	}

	if( COMPARE(token,"ReadMaterials")==0 ) { // could detoure to user-defined routine
		haveReadMaterials=factory->MaterialFactory(currentStream, &materialList, filemanager);
		OK
	}   

	if( COMPARE(token,"ReadTitle")==0 )	{
		(*currentStream).getline(title, 80,'\n');
		BETA_OUT<<"Title: "<<title << endl;
		OK
	}
		
	if(COMPARE(token,"setVerboseFlag") ==0) 
	{ 
		verboseFlag=(VerboseLevel)(atoi(localTokenList[1])) ;
		OK
	}
		
	if(COMPARE(token,"StoreAllMe") ==0) 
	{ 
		if( atoi(localTokenList[1])==0 )
			StoreAllMe=false;
		else
			StoreAllMe=true;
		BETA_OUT << "StoreAllMe " << StoreAllMe << endl;
		OK
	}

	if( COMPARE(token,"ReadDispFile")==0 ) {
		if(!mesh->readDisplacements(localTokenList[1]))
			FatalError("Error in ReadDispFile");
		OK
	}

	if( COMPARE(token,"ReadDispFileIntoCurrentSolution")==0) {
		if(!(readDisplacementsIntoCurrentSolution(localTokenList[1])))
			FatalError("Error in ReadDispFile");
		OK
	}

	if( COMPARE(token,"SetFindNodeTOLERANCE")==0 ) {
        mesh->setNodeSearchTolerance(atof(localTokenList[1]));
		findNodeTOLERANCE=atof(localTokenList[1]);//Eventually this needs to go away.
		BETA_OUT	<< "NEW findNodeTOLERANCE = " << DOUBLE_FORMAT << findNodeTOLERANCE << endl;
		cout		<< "NEW findNodeTOLERANCE = " << DOUBLE_FORMAT << findNodeTOLERANCE << endl;
		OK
	}

	if( COMPARE(token,"SolverSettings")==0 ) { 
		equations.readCommands(currentStream);
		OK
	}   

	if(COMPARE(token,"CalculateTotalGlobalForcesForNodeGroup") ==0) 
	{ 
		NodeSet nlist=mesh->getNodeSet(localTokenList[1]);
		BETA_OUT << "NodeSet: " << localTokenList[1] << endl;
		CalculateTotalGlobalForcesForNodeGroup(solution, nlist);
		OK
	}

	if(COMPARE(token,"AddToFileSearchPath") ==0) 
	{ 
		filemanager->AddDirectoryToSearchPath(localTokenList[1]);
		BETA_OUT << "Added this directory to file search path :" << localTokenList[1] << endl;
		OK
	}

	if(COMPARE(token,"UseMultiCore") ==0) 
	{ 
		BETA_OUT << "UseMultiCore (specified in input file): " << localTokenList[1] << endl;
		int numCoresRequested=1;
		if(COMPARE(localTokenList[1],"max") ==0)
			numCoresRequested=omp_get_max_threads();
		else
			numCoresRequested=atoi(localTokenList[1]);
		BETA_OUT << "MaxPossibleThreads (from omp_get_max_threads())=" << omp_get_max_threads() << endl;

		if( numCoresRequested==0){ //|| numCoresRequested==1 ){
			UseOpenMP=false;
			BETA_num_threads=1;
		}else{
			UseOpenMP=true;
			int MaxPossibleThreads=omp_get_max_threads();
			BETA_OUT << "This machine has " << MaxPossibleThreads << " core(s) that can be potentially used." << endl;
			if(numCoresRequested != MaxPossibleThreads){
				BETA_OUT << "The requested number of cores is not the same as the max possible cores." << endl;
			}
			if(numCoresRequested > MaxPossibleThreads) 
				EXIT_BETA("you are asking for more cores than max possible... exiting!");
			BETA_num_threads=numCoresRequested;
		}
		/*if(UseOpenMP && BETA_num_threads==1) {
			BETA_OUT << "Specified UseMultiCore but program is allowed only one thread, therefore turning off parallelization." << endl;
			UseOpenMP=false;
		}*/
		BETA_OUT << "UseOpenMP flag set to :" << UseOpenMP << endl;
		BETA_OUT << "BETA will use " << BETA_num_threads << " thread(s)." << endl;
		OK
	}

	if(COMPARE(token,"ParallelAssembly") ==0) 
	{
			
		if(COMPARE(localTokenList[1],"OFA") ==0){
			parallelAssembly_OFA = true;
			ParallelAssembly = true;
		}
		else if(COMPARE(localTokenList[1],"ESA") ==0){
			parallelAssembly_ESA = true;
			ParallelAssembly = true;
		}
		else{
			EXIT_BETA("Requested parallel assembly algorithm not recognized");
		}
		
		if(ParallelAssembly && UseOpenMP==false)
			EXIT_BETA("The UseOpenMP flag is currently false. Use the multicore command and check if multiple cores are available before attempting parallel assembly");
		if(ParallelAssembly && BETA_num_threads < BETA_num_threads_for_assembly)
			EXIT_BETA("Not enough cores to perform parallel assembly");

		BETA_OUT << "ParallelAssembly : " << localTokenList[1] << endl;
		OK
	}

    if(COMPARE(token,"ReadMPCNeutralFile") ==0) 
	{
		BETA_OUT << "Reading MPCs from Neutral MPC file :" << localTokenList[1] << endl;
        equations.readMPCNeutralFile(localTokenList[1]);
		OK
	}

	//restart feature commands
	IS_COMMAND(RestartParameters)
	IS_COMMAND(RestartBackupParameters)
	if( COMPARE(token,"restartDispfile")==0 ) {
		strcpy(restart_dispfilename,localTokenList[1]);
		if( doesFileExist(restart_dispfilename) == false){
			cout << "restart_dispfilename does not exist!\nrestart failed!" << endl;
			exit(1);
		}
		OK;
	}
	if( COMPARE(token,"restartIterationNumber")==0 ) {
		restart_iterationNumber=atoi(localTokenList[1]);
		OK;
	}
	if( COMPARE(token,"restart_backup_interval")==0 ) {
		restart_backup_interval=atoi(localTokenList[1]);
		OK;
	}
	//end restart feature commands



	if(currentStream!=originalStream) {
		BETA_OUT<<"Close alternate input stream"<<endl;
		(*currentStream).close();
		delete currentStream;
		currentStream=originalStream;
		BETA_OUT<<"Revert to default stream "<<*currentStream
			    <<" upon return to BasicModel::readCommands"<<endl;
	}

	if(foundMatch ==0){	
		BETA_OUT << "Command: " << token << endl;
		BETA_OUT << "Did not recognize command: " <<token<<endl;
		exit(1);
	}

	return foundMatch;
}
//========================================================================
//========================================================================
#define F1 setw(5)
#define F2 setw(13)<<setprecision(4)
//========================================================================
void BasicModel::ReadConstraints(istream * inStream)
{
	ConstraintFactory(inStream);
}
//=========================================================================
void BasicModel::ReadLoads(istream * inStream)
{
	LoadFactory(inStream);
	haveReadLoads=true;
}
//=========================================================================
void BasicModel::SetNumDofPerNode(istream * inStream)
{
Node::setNumDofPerNode(inStream, &BETA_OUT,mesh->numNodes, &mesh->node, totalNumDof);
return;
}
//=========================================================================
void readQuadPoints(int numQP, vector<Node> &qpList,istream *is);

void BasicModel::MapQuadPointsToMesh(int numElMappedMesh,int numQPPerElMappedMesh,char* fileName, char*mapFileName)
{
	istream *qpStream;
	int numQP = numElMappedMesh*numQPPerElMappedMesh;
	qpStream = openFileForInput(fileName);
	vector<Node> qpList(numQP);
	readQuadPoints(numElMappedMesh,qpList,qpStream);

	int numEl = mesh->numElements;

	vector<vector<Node> > elementQPListContainer(numEl);

	omp_set_num_threads(BETA_num_threads);
	#pragma omp parallel               
	{
	#pragma omp for		
	for(int i=0; i<numEl; i++)
	{
		BasicElement *e;
		e = &mesh->element[i];
		int count =0;
	
		for(int j=0; j<numQP;j++){
			double pointCoord[3];
			double localCoord[3]={0.0,0.0,0.0};
			bool pointInside = false;
			pointCoord[0] = qpList[j].x;
			pointCoord[1] = qpList[j].y;
			pointCoord[2] = qpList[j].z;
			pointInside = e->SearchForPointInsideElement(pointCoord,localCoord);

			int pointMatNum = qpList[j].getNumDof(); //mat num is being stored in numdof
			bool sameMat = false;
			if(pointMatNum == e->getMaterialGroup()){
				sameMat = true;
			}

			if(pointInside){
				if(!sameMat){
					cout<<"Error: Point is not same mat group as element it is inside"<<endl;
					exit(1);
				}
				elementQPListContainer[i].push_back(qpList[j]);
				
				elementQPListContainer[i][count].x =  localCoord[0];
				elementQPListContainer[i][count].y =  localCoord[1];
				elementQPListContainer[i][count].z =  localCoord[2];
				count ++;

			}

		}

	}
	} //end of parallel bracket


	//output extrapolation factors and index of QPs to file

	ostream *mapFile;
	mapFile = filemanager->OpenOutputStream(mapFileName);
	SET_SCIENTIFIC(*mapFile);

	for(int i=0;i<numEl;i++){
		BasicElement *e;
		GaussPointList gaussPointList;
		e = &mesh->element[i];
		e->getQuadraturePoints(gaussPointList);
		int nDims = mesh->numDims;

		ExtrapolationMatrix extrapMat(nDims,elementQPListContainer[i],gaussPointList,2);

		int numQPForEl = elementQPListContainer[i].size();
		*mapFile<<numQPForEl<<endl;
		for(int j=0;j<numQPForEl;j++){
			*mapFile<<elementQPListContainer[i][j].getNodeNum()<<"\t";	
		}
		*mapFile<<endl;

		int numRows = extrapMat.Size();
		for(int j=0;j<numRows;j++){
			for(int k=0; k<numQPForEl; k++){
				*mapFile<<DOUBLE_FORMAT<<extrapMat[j][k]<<"\t";
			}
			*mapFile<<endl;
		}
	}
	filemanager->CloseOutputStream(mapFile);
	
return;
}
//=======================================================================
void BasicModel::MapQuadFile(int numEl, char *mapFileName, char *originQuadFileName, char *destQuadFileName)
{

//read in origin quad file
	ifstream *originStream;
	originStream = openFileForInput(originQuadFileName);
	int dum=0;
	int mat=0;
	int numQP=0;
	int numCol=0;
	vector<vector<double> > quadVals;
	for(int i=0; i< numEl; i++){
		*originStream >> dum >> dum >> numQP >> numCol;
		for(int j=0; j<numQP;j++){
			double value=0.0;
			vector <double> tempQP;
			for(int k=0; k<numCol; k++){
				*originStream >> value;
				tempQP.push_back(value);
				
			}
			quadVals.push_back(tempQP);
		}
	}

	filemanager->CloseInputStream(originStream);

//read in map for each element
	ifstream *mapFile;
	mapFile = openFileForInput(mapFileName);

	ostream *destFile;
	destFile = filemanager->OpenOutputStream(destQuadFileName);
	SET_SCIENTIFIC(*destFile);

	int numColMap;
	for(int i=0; i<mesh->numElements; i++)
	{
		BasicElement *e;
		GaussPointList gaussPointList;
		e = &mesh->element[i];
		e->getQuadraturePoints(gaussPointList);
		int nQP = gaussPointList.totalNumIPs;

		*destFile<<i<<'\t'<<e->getMaterialGroup()<<'\t'<<nQP<<'\t'<<numCol<<endl;
		*mapFile >> numColMap;
		vector<int> mapIndex;
		mapIndex.reserve(numColMap);
		
		int index=0;
		double extrapVal =0.0;
		for(int j=0;j<numColMap;j++)
		{
			*mapFile >> index;
			mapIndex.push_back(index);
		}
		
		for(int j=0;j<nQP;j++)
		{
			vector<double> extrapRow;
			extrapRow.reserve(numColMap);
			for(int k=0;k<numColMap;k++)
			{
				*mapFile>>extrapVal;
				extrapRow.push_back(extrapVal);
			}
			//construct local vector of values to map over
			for(int m=0;m<numCol;m++)
			{
				vector<double> tempLocalCol;
				for(int n=0;n<numColMap;n++)
				{
					tempLocalCol.push_back(quadVals[mapIndex[n]][m]);
				}
				double extrapolatedVal = inner_product(tempLocalCol.begin(),tempLocalCol.end(),extrapRow.begin(),0.0);
				*destFile<<DOUBLE_FORMAT<<extrapolatedVal<<'\t';
			}
			*destFile<<endl;
		}
	}
	filemanager->CloseOutputStream(destFile);

	return;
}
//=======================================================================
// DECIDE WHAT TO DO WITH THIS !
void BasicModel::SetElementProperty(istream * inStream)  //This will phase out
{

// Maybe this should be a static element method !
// ..will have to pass the list of material pointers
//........Do not do yet... if property classes are used (like
// material classes are currently used), logic would change.

	//JV040507 - modified to let derived models implement their own method to set element properties
	// so, as of now, each model will implement its own ProcessSetElementProperty(string command)

	char command[80];
	(*inStream)>>command;
	banner(command, BETA_OUT);

	ProcessSetElementProperty(command, inStream);
}

void BasicModel::ProcessSetElementProperty(string command, istream * inStream)
{
	int intCommand=0 ;
	enum {SET_ELEMENT_MATERIAL,SET_ELEMENT_TEMPERATURE,SET_QUADRATURE_ORDER};
BETA_OUT<<"\nCommand processor = 'BasicModel::ProcessSetElementProperty'          "<<command<<endl;

ChangeToUpper(command)

if(command == "SETELEMENTMATERIAL" || command == "SELECTELEMENTMATERIAL" ||
   command == "READMATERIALGROUP()" )
                          {intCommand=SET_ELEMENT_MATERIAL;}
else if(command == "SETELEMENTTEMPERATURE")
                          {intCommand=SET_ELEMENT_TEMPERATURE;}
else if(command == "SETINTEGRATIONORDER")
                          {intCommand=SET_QUADRATURE_ORDER;}
else {
	FatalError("Cannot recognize this SetElementProperty command:" + command);
}

int i,first,last,increment;
double temperature;
Material *matPoint;	
int matGroupNum,quadOrder;

char *localTokenList[20];
int  numberOfTokens;

int dataIndex=0;

while (2==2) {				
		if( getLineAndTokenize(inStream,"end",localTokenList, numberOfTokens)==1) {
			return;
		}
		if(numberOfTokens >= 2 && COMPARE(localTokenList[0],"all")==0 ){
			first=0;
			last=mesh->numElements-1;
			increment=1;
			dataIndex=1;
		}else if(numberOfTokens >= 4){
			first=atoi(localTokenList[0]);
			last=atoi(localTokenList[1]);
			increment=atoi(localTokenList[2]);
			dataIndex=3;
		}else{
			first=atoi(localTokenList[0]);
		}
		if(first<0) return;

	//Error checking...
	if(first < 0 || first    > mesh->numElements || first > last
				 || last     > mesh->numElements || increment<1 ){
		EXIT_BETA("Fatal error...input is bad\nFirst, last, increment="<<first<<" "<<last<<" "<<increment);
	}

	BasicElement *e=0;

switch(intCommand) {
	case SET_ELEMENT_MATERIAL:
		matGroupNum=atoi(localTokenList[dataIndex]);
		BETA_OUT << "\tGroup = " << matGroupNum << " for elements "
			  << first <<" to "<<last << " by " << increment << '\n';				
		matPoint = &(materialList[matGroupNum]);				
		for(i=first;i<=last;i+=increment) {
			e=&mesh->element[i];
			e->setMaterial(matPoint);
            e->setMaterialGroup(matGroupNum);}
		break;
	case SET_ELEMENT_TEMPERATURE:
		temperature=atof(localTokenList[dataIndex]);
		BETA_OUT << "\tTemperature = " << temperature << " for elements "
			 << first <<" to "<<last << " by " << increment << '\n';
        for(i=first;i<=last;i+=increment) {
			e=&mesh->element[i];
			e->setTemperature(temperature);}
		break;
	case SET_QUADRATURE_ORDER:    //obsolete
		quadOrder=atoi(localTokenList[dataIndex]);
		BETA_OUT << "\tQuadrature order = " << quadOrder << " for elements "
			  << first <<" to "<<last << " by " << increment << '\n';								
		for(i=first;i<=last;i+=increment) {
			e=&mesh->element[i];
			e->setIntegrationOrder(quadOrder);}
		break;
	default:
		BETA_OUT<<"Match not found"<<endl;
		exit(1);
}

}//end of while (2==2)
} //end of setElementProperty
//--------------------------------------------------------------------

/* // taken out when added parallel version
void BasicModel::allocateForAnalysis(ElementWorkspace * &elemWorkspace, bool allocateEquations)
{
	if(allocationState==1) return;
	CreateAndAllocateWorkspace(elemWorkspace);
	assignElementWorkspacePointers(elemWorkspace);

	if(allocateEquations){
		equations.setupStorage(totalNumDof);
		allocateEquationsMatrix();
		displayTime("Time to allocate equations", BETA_OUT);
	}
}
*/
//=======================================================================
void BasicModel::allocateForAnalysis(bool allocateEquations)
{
	if(allocationState==1) return;

	
	FileManager *filemanager1=0;
	filemanager1 = new FileManager;
	filemanager1->initialize(filemanager->inputFilename);

	allocateArrays();

	setAssemblerPointers();
	assembler->allocateAndInitialize(this);
	
	//This call is required to allocate the equations. Allocation of equations requires
	//certain pointers in an element to be set (ie. dofList).
	//For ParallelAssembly elements pointers will be re-assigned once the elements are
	//divided among the threads.
	assignElementWorkspacePointers(assembler->threadWorkspaceList[0].eWorkspace);

	if(allocateEquations){
	BETA_num_threads_for_assembly = BETA_num_threads;
	equations.numberOfOpenMPThreads=BETA_num_threads_for_assembly;
	equations.setupStorage(totalNumDof);
	allocateEquationsMatrix();
	displayTime("Time to allocate equations", BETA_OUT);
	}

	assembler->makeAssignments(totalNumDof);
    CreateElementISV();

if(allocateEquations) allocationState = 1; //allocationstate is 1 (meaning ready for a complete analysis) only when all arrays, workspaces and equations have been allocated
}

//=======================================================================


void BasicModel::allocateArrays()
{
	if(!appliedForces)		Assert(appliedForces   = new double[totalNumDof]);
	if(!currentSolution)	Assert(currentSolution = new double[totalNumDof]);
	if(!internalForce)		Assert(internalForce   = new double[totalNumDof]);
    if(!zeroVector)			Assert(zeroVector      = new double[totalNumDof]);  

	for(int i=0; i<totalNumDof; i++) {
		currentSolution[i]	= 0.;
		appliedForces[i]	= 0.;
		internalForce[i]	= 0.;
		zeroVector[i]		= 0.;
	}
}
//=======================================================================
void BasicModel::CreateAndAllocateWorkspace(ElementWorkspace * &elemWorkspace)
{
	displayTime("", BETA_OUT);
	CreateWorkspace(elemWorkspace);
	SetElementWorkspaceLimits(elemWorkspace);
	BETA_OUT << "ElementWorkspace array size limits" << endl;
	//elemWorkspace->PrintLimits(&BETA_OUT);
	allocateWorkspace(elemWorkspace);
	displayTime("Time to allocate element workspace", BETA_OUT);
}
//=======================================================================

void BasicModel::CreateWorkspace(ElementWorkspace * &elemWorkspace)
{
	if(elemWorkspace==0) 
		elemWorkspace=new ElementWorkspace;
	else{
		delete elemWorkspace;
		elemWorkspace=new ElementWorkspace;
	}
}
//=======================================================================
void BasicModel::CreateElementISV()
{
    //Allocates the ISVs for any element that has not had an ISV allocated yet.
	for(int i=0; i<mesh->numElements; i++) {
        if(!(&mesh->element[i])->checkForAllocatedISVs()){
            (&mesh->element[i])->allocateISVs();
        }
	}
}
//=======================================================================
void BasicModel::SetElementWorkspaceLimits(ElementWorkspace * elemWorkspace)
{
	elemWorkspace->WorkspaceMaxQuadPoints		=64;
	elemWorkspace->WorkspaceMaxStresses			=10;
	elemWorkspace->WorkspaceMaxFluxes			=3;
	elemWorkspace->WorkspaceMaxNumberDof		=96;	
	elemWorkspace->WorkspaceMaxNodes			=32;
	elemWorkspace->WorkspaceMaxNumberOfStrains	=10;
}
//=======================================================================
void BasicModel::allocateWorkspace(ElementWorkspace * &elemWorkspace)
{
//Warning:  should change to always initialize memory to zero after allocation
// It is probably preferable to keep the memory contiguous => need to allocate
// big block and then set pointers to pieces

	if(elemWorkspace->allocated) return;

	elemWorkspace->equations=&equations;
	elemWorkspace->mesh=mesh;

	int i;
	
	if(!haveReadMesh)      FatalError("Mesh not read.\n");
	if(!haveReadMaterials) FatalError("Materials not read.\n");
//	if(allocateEquations && !haveReadLoads)     FatalError("Loads not read.\n");
	if(!haveReadLoads)     FatalError("Loads not read.\n");
	
	int WorkspaceMaxQuadPoints		=elemWorkspace->WorkspaceMaxQuadPoints;
	int WorkspaceMaxStresses		=elemWorkspace->WorkspaceMaxStresses;
	int WorkspaceMaxFluxes			=elemWorkspace->WorkspaceMaxFluxes;
	int WorkspaceMaxNumberDof		=elemWorkspace->WorkspaceMaxNumberDof;	
	int WorkspaceMaxNodes			=elemWorkspace->WorkspaceMaxNodes;
	int WorkspaceMaxNumberOfStrains	=elemWorkspace->WorkspaceMaxNumberOfStrains;

	Assert(elemWorkspace->quadPointStresses    = new double * [WorkspaceMaxQuadPoints]);
	Assert(elemWorkspace->rotatedStresses      = new double * [WorkspaceMaxQuadPoints]); 
	Assert(elemWorkspace->extrapolatedStresses = new double * [WorkspaceMaxQuadPoints]);      

	Assert(elemWorkspace->quadPointStrains     = new double * [WorkspaceMaxQuadPoints]);
	Assert(elemWorkspace->rotatedStrains       = new double * [WorkspaceMaxQuadPoints]);
	Assert(elemWorkspace->extrapolatedStrains  = new double * [WorkspaceMaxQuadPoints]);

	Assert(elemWorkspace->quadPointFluxes		= new double * [WorkspaceMaxQuadPoints]);
	Assert(elemWorkspace->rotatedFluxes			= new double * [WorkspaceMaxQuadPoints]);
	Assert(elemWorkspace->extrapolatedFluxes   = new double * [WorkspaceMaxQuadPoints]);

	Assert(elemWorkspace->quadPointGradients	= new double * [WorkspaceMaxQuadPoints]);
	Assert(elemWorkspace->rotatedGradients		= new double * [WorkspaceMaxQuadPoints]);
	Assert(elemWorkspace->extrapolatedGradients	= new double * [WorkspaceMaxQuadPoints]);

	for(i=0; i<WorkspaceMaxQuadPoints; i++)
	{
	 Assert(elemWorkspace->quadPointStresses[i]    = new double [WorkspaceMaxStresses]);
	 Assert(elemWorkspace->rotatedStresses[i]      = new double [WorkspaceMaxStresses]); 
	 Assert(elemWorkspace->extrapolatedStresses[i] = new double [WorkspaceMaxStresses]); 

	 Assert(elemWorkspace->quadPointStrains[i]     = new double [WorkspaceMaxStresses]);
	 Assert(elemWorkspace->rotatedStrains[i]       = new double [WorkspaceMaxStresses]);
	 Assert(elemWorkspace->extrapolatedStrains[i]  = new double [WorkspaceMaxStresses]);

	 Assert(elemWorkspace->quadPointFluxes[i]		= new double [WorkspaceMaxFluxes]);
	 Assert(elemWorkspace->rotatedFluxes[i]			= new double [WorkspaceMaxFluxes]);
	 Assert(elemWorkspace->extrapolatedFluxes[i]    = new double [WorkspaceMaxFluxes]);

	 Assert(elemWorkspace->quadPointGradients[i]	= new double [WorkspaceMaxFluxes]);
	 Assert(elemWorkspace->rotatedGradients[i]		= new double [WorkspaceMaxFluxes]);
	 Assert(elemWorkspace->extrapolatedGradients[i]	= new double [WorkspaceMaxFluxes]);

	}
	//...........
	
	Assert(elemWorkspace->b      = new double *[WorkspaceMaxNumberOfStrains]);
	Assert(elemWorkspace->b1      = new double *[WorkspaceMaxNumberOfStrains]);
	Assert(elemWorkspace->b2      = new double *[WorkspaceMaxNumberOfStrains]);
	for(i=0; i<WorkspaceMaxNumberOfStrains; i++){
		Assert(elemWorkspace->b[i] = new double [WorkspaceMaxNumberDof]); 
		Assert(elemWorkspace->b1[i] = new double [WorkspaceMaxNumberDof]); 
		Assert(elemWorkspace->b2[i] = new double [WorkspaceMaxNumberDof]); 
	}

	//Assert(elemWorkspace->dofList = new int [WorkspaceMaxNumberDof]);

	//..............................ESA Assembly.......................................
	if(parallelAssembly_ESA) {
		elemWorkspace->Fe        = 0;
		elemWorkspace->initialFe = 0;
		elemWorkspace->Ke	     = 0;
		elemWorkspace->dofList   = 0;
		elemWorkspace->Ke_original = 0;
	}
	//.................................................................................
	else{
		Assert(elemWorkspace->Fe        = new double [WorkspaceMaxNumberDof]);
		Assert(elemWorkspace->initialFe = new double [WorkspaceMaxNumberDof]);
		Assert(elemWorkspace->Ke	    = new Matrix(WorkspaceMaxNumberDof,WorkspaceMaxNumberDof));
		Assert(elemWorkspace->dofList   = new int [WorkspaceMaxNumberDof]);
		Assert(elemWorkspace->Ke_original      = new Matrix(WorkspaceMaxNumberDof,WorkspaceMaxNumberDof));
	}

	Assert(elemWorkspace->F_original       = new double [WorkspaceMaxNumberDof]);
	
	if (!StoreAllMe){
		elemWorkspace->StoreAllMe=false;
		if(parallelAssembly_ESA){
			elemWorkspace->Me = 0;
		}
		else{
			Assert(elemWorkspace->Me               = new Matrix(WorkspaceMaxNumberDof,WorkspaceMaxNumberDof));
		}
	}else{//if going the 'parallel' way, this section needs to be looked at in more detail.
		//storing all the Me in memory
		if(parallelAssembly_ESA){
		BETA_OUT<<"Have not addressed storeAllMe option with parallel assembly... exiting"<<endl;
		exit(1);
		}

		elemWorkspace->StoreAllMe=true;
		Assert(elemWorkspace->MeList=new Array<Matrix>);
		elemWorkspace->MeList->reserve(mesh->numElements);
		Matrix *Melem_i=0;
		BasicElement *e=0;
		for( i=0; i<mesh->numElements; i++)
		{
			e=&mesh->element[i];
			e->setDofListPointer(elemWorkspace->dofList);
			e->calculateDofList();//needed to calculate numDof
			int elemRank = e->getNumDof();
			Melem_i = new Matrix(elemRank,elemRank);
			elemWorkspace->MeList->add(*Melem_i);
		}
	}
	
	elemWorkspace->allocated=true;
}
//========================================================================

void BasicModel::assignElementWorkspacePointers(ElementWorkspace * &elemWorkspace)
{
	//set element pointers
	BasicElement *e=0;
	for(int i=0; i<mesh->numElements; i++) //Check element static array sizes
	{
		e=&mesh->element[i];
		e->setPointers(elemWorkspace);

	}
}
//========================================================================
void BasicModel::allocateEquationsMatrix()
{
	//have to apply restraints before calling this function!!!
	BasicElement *e=0;
	for(int i=0; i<mesh->numElements; i++) //Check element static array sizes
	{
		e=&mesh->element[i];
		e->specifyNonZeroLocation();
	}
	equations.allocate();
	equations.zeroMatrix();
}
//========================================================================
double BasicModel::CalculateResidual(double *theSolution)
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

		e->update("F", theSolution);

		int numDof   = e->getNumDof();
		int *dofList = e->calculateDofList();
		elementForce = e->getForceVector();
		equations.assembleResultantVector(dofList,numDof,
			                              elementForce);
	}//end of loop on elements

	return equations.calculateResidualNorm();
	//BETA_OUT<<"Residual norm= "<<residualNorm<<endl;
}//end of CalculateResidual
//===============================================================
void BasicModel::getGlobalForces( double *internalForce,
							    double *theSolution,
								double *mappedResidualVector,
								double &maxResidual,
								int    &dofForMaxResidual)
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
		e->update("F", theSolution);

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
void BasicModel::getGlobalForces( double *internalForce, double *theSolution, ElementGroup *elist, int NonZeroStressComponent)
{
	BasicElement *e;
	int i;

	double * elementForce;
	for(i=0;i<totalNumDof;i++) internalForce[i]=0;

	equations.zeroResultantVector();	

	int numE=elist->getNumElements();

	for(i=0;i<numE;i++) {
		e=&(*elist)[i];
		e->initializeSummary(); // clear worksapce ?

		//((ElasticityElement3D*)e)->ForceSingleNonZeroStressComponent=NonZeroStressComponent;
		//currently have to force this typecast... see discussion in towtracer on how to handle ForceSingleNonZeroStressComponent
		e->update("F", theSolution);
		int numDof   = e->getNumDof();
		int *dofList = e->calculateDofList();
		elementForce = e->getForceVector();

		equations.assembleResultantVector(dofList,numDof,
										  elementForce);
		for(int ii=0; ii<numDof; ii++)// InternalForce[] is not mapped
		{ internalForce[dofList[ii]] += elementForce[ii];	}

	}
/* //uncomment this block to write Global Forces
OS<<"Global Forces (restraint forces): " << endl;
OS<<"nodenum dof force" << endl;
Node *n;
	for(i=0;i<num;i++) {
		ci=componentList[i];
		if (ci->component->isActiveModel) {
			OS << " Model: " << ci->component->modelName << endl;
			for(int ii=0; ii<ci->component->mesh->numNodes; ii++){
				n=&ci->component->mesh->node[ii];
				if(n->isActive){
					for(j=0; j<n->numDof; j++){
						int gdof=n->getFirstDof() +j;
						if ( fabs(internalForce[gdof]) >= 1e-4){
							OS << n->nodeNum  << '\t' << j << '\t' << internalForce[gdof] << endl;
						}
					}
				}
			}
		}
	}
*/

}//end of getGlobalForces
//===============================================================
#undef F1
#undef F2
//========================================================================
void BasicModel::DoAnalysis(istream * inStream)
{
	BETA_OUT<<"Start DoAnalysis"<<endl;
	int i;

	FormWriter fr(BETA_OUT);
	
	createAssembler();
	allocateForAnalysis();
	
	// Send loads to equation object
	equations.zeroLoadVector();
	
	for( i=0; i<(int)loads.size(); i++){
		loads[i]->applyTo(equations,&mesh->node,mesh->numNodes, this);
	}

	checkModel();
	SetInitialState();
	InitializeOptionalOutput();
		
	BETA_OUT<<"Assemble stiffness matrix"<<endl;
	displayTime("", BETA_OUT);
	
    assembler->assembleInitialF();
	displayTime("Time to assembleInitialF", BETA_OUT);
	assembler->assembleK(NULL);
	displayTime("Time to assembleK", BETA_OUT);

	BETA_OUT<<"Start solving equations (" << totalNumDof << ") :"<<endl;
	displayTime("", BETA_OUT);

	assembler->postAssemblyOperations(this);

	cout << "solving (" << totalNumDof << " eqns)"<<endl;
	solution = equations.solve();
	displayTime("Time to solve", BETA_OUT);

	for(i=0; i<totalNumDof; i++){currentSolution[i] += solution[i];}

//commenting out the following block.. use ReadOptionalOutput for any getting any output
//	ofstream *os=filemanager->OpenOutputStream("disp");
//	mesh->printNodalData(currentSolution, "displacements", os );
//	filemanager->CloseOutputStream(os);
	
	haveSolution = true;

	WriteOptionalOutput();
	FinalizeOptionalOutput();

	if( BasicElement::getAnalysisType() == 1) iterate();
} // end of DoAnalysis
//=========================================================================
void BasicModel::DoPostProcess(istream *inStream)
{
cout<<"Start PostProcessing"<<endl;
SetInitialState();
WriteOptionalOutput();
//FinalizeOptionalOutput();
}
//=========================================================================
bool BasicModel::readDisplacementsIntoCurrentSolution(char *filename)
{
	createAssembler();
    //Don't allocate the equations.
	allocateForAnalysis(true);

	//this is only to make sure that the file exists in the right location and that beta can open it
	ifstream *is=filemanager->OpenInputStream(filename);
	filemanager->CloseInputStream(is);

	setVectorFromNodalValueFile(filename,currentSolution, mesh->numNodes, mesh->node);
	return true;	
}
//=========================================================================
void BasicModel::iterate()     // Move this to nonlinear model class
{
 double *currentSolutionInEquations;
 int    dofForMaxResidual;
 double maxResidual;
 int i;

	banner("Start iteration", BETA_OUT);
	double * deltaSolution;
	double * mappedResidualVector = new double[totalNumDof];
	double * Ro                   = new double[totalNumDof];
	
	//Eventually set in input data
	double scaleFactor, So, S1, RoRo, RoR1, R1R1;
	int variableRelaxation = None;  // other options: Felippa    Crisfield;

scaleFactor = 1.;	
variableRelaxation = None;

	//printNodalData(appliedForces, "Applied Forces" );

	for(int iteration = 0; iteration<400; iteration ++) {
//	for(int iteration = 0; iteration<20; iteration ++) {
		//first iteration (i.e. iteration=0 refers to the initial linear solution (from doAnalaysis())
		//if the initial linear solution is good enough, it does not solve anymore! 

		//handles restart
		if(doRestart && !restartInitiated){
			setVectorFromNodalValueFile(restart_dispfilename, currentSolution, mesh->numNodes, mesh->node);
			currentSolutionInEquations = equations.getSolution();
			copyVector(currentSolution,currentSolutionInEquations,
			       totalNumDof);
			iteration=restart_iterationNumber;
			restartInitiated=true;
		}

		getGlobalForces(internalForce, currentSolution,mappedResidualVector,
		                maxResidual,dofForMaxResidual );
		
		BETA_OUT<<" Iteration:"<<iteration<< "   Max. abs(Residual) = "<< 
				maxResidual<<":   ";

		cerr <<" Iteration:"<<iteration<< "   Max. abs(Residual) = "<< 
				maxResidual<< endl;
		

		if(iteration==0||variableRelaxation==None) BETA_OUT<<endl;
		
		if(maxResidual< maxAllowableResidual){BETA_OUT<<endl; break;}

		//handles creating backup files for restart
		if(createFilesForRestart  && iteration%restart_backup_interval==0){
			char fndump[256];
			sprintf(fndump,"disp_NL.%d",iteration);
			cout << "backing up to:" << fndump << endl;
			ofstream *os=filemanager->OpenOutputStream(fndump);
			*os << "displacements" << endl;
			mesh->printNodalData(currentSolution, "  ", os );
			filemanager->CloseOutputStream(os);
		}

		assembler->assembleK(currentSolution);

//Check this...what should I send?
		deltaSolution=equations.incrementalSolutionUsingResidual(currentSolution);
// Eventually change this ... make local copy rather than
// perhaps accidentally changing data in equations

		//printNodalData(deltaSolution, "delta Nodal Displacements" );
//This is where variable deltaSolution would be scaled & then sent 
///*
		if(iteration>0){
		for(i=0; i< totalNumDof; i++) {
//			deltaSolution[i] = scaleFactor * deltaSolution[i]; 
			deltaSolution[i] = 1.0 * deltaSolution[i]; 
			}
		}//if 
//*/
		equations.addIncrementalDisplacements(deltaSolution);

//Get solution from equations and make copy 
       currentSolutionInEquations = equations.getSolution();
	   copyVector(currentSolutionInEquations,currentSolution,
			       totalNumDof);

		//printNodalData(currentSolution,"currentSolution..after update");

		
		getGlobalForces(internalForce, currentSolution,mappedResidualVector,
					    maxResidual,dofForMaxResidual);
		
		R1R1 = dotProduct( mappedResidualVector, mappedResidualVector, totalNumDof);
		
		S1   = dotProduct( mappedResidualVector, currentSolution, totalNumDof);
		
		scaleFactor = 1.0;
		if(iteration>0 && variableRelaxation != None) {
			if(variableRelaxation==Crisfield) {
				scaleFactor = -So/(S1 - So);
				BETA_OUT<<"ScaleFactor = "<<scaleFactor<<endl;
			}
			if(variableRelaxation==Felippa) {
				RoR1 = dotProduct( Ro, mappedResidualVector, totalNumDof);
				scaleFactor = (RoRo-RoR1)/(RoRo-2. * RoR1 + R1R1);
				BETA_OUT<<"ScaleFactor = "<<scaleFactor<<endl;
			}
			
			if(scaleFactor<.2 || scaleFactor> 2.0) { scaleFactor = .5; }
			scaleFactor = scaleFactor - 1.;
			//for(i=0; i< totalNumDof; i++) {
            // Fix this...see "standard relaxation

			//	currentSolution[i] += scaleFactor * deltaSolution[i]; 
			//}
		}
		
		So = S1;
		RoRo = R1R1;
		//Save residual vector
		copyVector( mappedResidualVector, Ro, totalNumDof);
		
		//printNodalData(currentSolution, "Nodal Displacements" );
	}//end of iterate loop
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	printNodalData(currentSolution, "Nodal Displacements" );
	ofstream *os=filemanager->OpenOutputStream("disp_NL");
	mesh->printNodalData(currentSolution, "displacements", os );
	filemanager->CloseOutputStream(os);

	banner("Element stresses",BETA_OUT);
	getGlobalForces(internalForce, currentSolution,mappedResidualVector,
					maxResidual,dofForMaxResidual);
	
	mesh->printNodalData(internalForce, "Global force vector" );

    delete [] mappedResidualVector; //xt062901	
	delete [] Ro;//xt062901	
}
//=========================================================================
void BasicModel::ReadMesh(istream *inputStream)
{
	mesh->setFileManager(filemanager);
	mesh->ReadMesh(inputStream);
	haveReadMesh=true;
}
//=========================================================================
void BasicModel::ReadElementOrientations(istream* inputStream)
{
    for(int i=0;i<mesh->numElements;++i){
       (&mesh->element[i])->readElementOrientation(inputStream);
    }
}
//=========================================================================
void BasicModel::CreateNodeGroup(istream *inputStream)
{
    string command;

	while(true) {
        //Command should either be "exitCreateNodeGroup" or 
        //the name of a NodeSet to create or modify
        (*inputStream) >> command;

        //Check for exit string...
        if(COMPARE(command.c_str(),"exitCreateNodeGroup") == 0){
            BETA_OUT << "Finished Reading Node Sets - summary of current sets:" << endl;
            for(map<const string,NodeSet>::const_iterator Pair = mesh->getNodeSets().begin(); Pair != mesh->getNodeSets().end();++Pair){
                BETA_OUT << endl << "Node Set: " << Pair->first;
                Pair->second.print(BETA_OUT);
            }
            return;
        }

        //Note, getAndCreateNodeSet returns a reference to the nodeset if it
        //already exists, or makes a new nodeset and returns a reference to it
        //if it does not already exist.  This allows you to create a complex
        //nodeset by adding multiple types of geometry.  
        //Names are not case sensitive.

        mesh->getNodeSet(command).read(inputStream,mesh);

	} //end of while
}
//=========================================================================
void BasicModel::CreateElementGroup(istream *inputStream)
{
	char *localTokenList[20];
	int  numberOfTokens;

	ElementGroup *elementgroup=0;

	while(2==2) {
		if( getLineAndTokenize(inputStream,"exitCreateElementGroup",localTokenList, numberOfTokens)==1) {
			return;
		}
		string groupname=localTokenList[0];

		if( COMPARE(localTokenList[1],"ElementNumber")==0 ){
			int elemnumber=-1;
			if(numberOfTokens != 3)FatalError("Error in specifying Elementgroup with 'ElementNumber'.\n");
			elemnumber=atoi(localTokenList[2]);
			if(elemnumber < 0 || elemnumber >= mesh->numElements)FatalError("Invalide Element number specified for Elementgroup.\n");
			elementgroup= new ElementGroup;
			elementgroup->reserve(1);
			BasicElement *e;
			e=&(mesh->element)[elemnumber];
			elementgroup->add(*e);
			if(elementgroup == 0) FatalError("Error in creating ElementGroup.\n");
		}
		else if( COMPARE(localTokenList[1],"Elematfile")==0 ){
			if(numberOfTokens != 3)FatalError("Error in specifying Elementgroup with 'Elematfile'.\n");
			string elematfilename=localTokenList[2];
			elementgroup=CreateElementGroupFromElematFile(elematfilename, mesh, filemanager);
			if(elementgroup == 0) FatalError("Error in creating ElementGroup.\n");
		}
		else
			EXIT_BETA("Cannot recognize the command for creating ElementGroup");


		BETA_OUT << "Created ElementGroup called *" << localTokenList[0] << "* containing " << elementgroup->getNumElements() << " element(s). Listing the elementnumbers:" << endl;
		for(int i=0;i<elementgroup->getNumElements();i++){
			if(i && i%10==0) 
				BETA_OUT << endl;
			BETA_OUT << (*elementgroup)[i].elementNumber << '\t';
		}
		BETA_OUT << endl;

		elementgroup->SetName(groupname.c_str());
		mesh->ElementGroupsList.add(*elementgroup);

	} //end of while
}
//=========================================================================
void BasicModel::setSolutionVectorFromMeshDisplacements(BasicMesh *mesh)
{
	solution=equations.getSolution();
	setVectorFromMeshDisplacements(mesh, solution);
} 
//========================================================================
bool BasicModel::RestartParameters(ifstream * defaultStream)
{
	char *localTokenList[20];
	int  numberOfTokens;
	int foundMatch;

	ifstream *currentStream;
    BETA_OUT<<"Current stream in 'BasicModel::RestartParameters' ="<<(*defaultStream)<<endl;
	currentStream = defaultStream;

	while(2==2) {
		foundMatch=0;		
		currentStream = defaultStream;		
		if( getLineAndTokenize(currentStream,"end",localTokenList, numberOfTokens)==1) {
			goto _ExitChoose;
		}
	
		foundMatch=processCommand(localTokenList,currentStream, numberOfTokens); // decendent class first

		if(foundMatch ==0) {	exit(1);}
		
	} //end of while

	_ExitChoose:
	doRestart=true;
	BETA_OUT << "Finished reading Restart Parameters"<<endl; 
	return true;
}
//=========================================================================
void BasicModel::initializeRestartFlags()
{
	isRestart = false;
	createFilesForRestart=false;
	doRestart	=false;
	restartInitiated=false;
	restart_backup_interval=5;
	restart_iterationNumber=0;
}
//--------------------------------------------------------------------
bool BasicModel::RestartBackupParameters(ifstream * defaultStream)
{
	char *localTokenList[20];
	int  numberOfTokens;
	int foundMatch;

	ifstream *currentStream;
    BETA_OUT<<"Current stream in 'BasicModel::RestartBackupParameters' ="<<(*defaultStream)<<endl;
	currentStream = defaultStream;

	while(2==2) {
		foundMatch=0;		
		currentStream = defaultStream;		
		if( getLineAndTokenize(currentStream,"end",localTokenList, numberOfTokens)==1) {
			goto _ExitChoose;
		}
	
		foundMatch=processCommand(localTokenList,currentStream, numberOfTokens); // decendent class first

		if(foundMatch ==0) {	exit(1);}
		
	} //end of while
	createFilesForRestart=true;

	_ExitChoose:
	BETA_OUT << "Finished reading Restart Backup Parameters"<<endl; 
	return true;

}
//========================================================================
void BasicModel::CalculateTotalGlobalForcesForNodeGroup(double * theSolution, const NodeSet &nlist)
{
	int    dofForMaxResidual;
	double maxResidual;
	double * mappedResidualVector = new double[totalNumDof];
	getGlobalForces(internalForce, theSolution, mappedResidualVector,
					maxResidual,dofForMaxResidual );

	int numN=nlist.getNumNodes();
	int numDof=nlist[0]->getNumDof();
	double *total=0;
	total=new double[numDof];
	int i,j;
	for(i=0;i<numDof; i++) total[i]=0.0;

	for(i=0;i<numN; i++){
		for(j=0;j<numDof; j++){
			total[j] += internalForce[nlist[i]->getFirstDof()+j];
		}
	}

	BETA_OUT << "Total Global Forces from nodegroup:" << endl;
	for(i=0;i<numDof; i++) {
		BETA_OUT << DOUBLE_FORMAT << total[i] << endl;
	}

	delete [] total;
	delete [] mappedResidualVector;

}

//==========================================================
bool BasicModel::ReadOptionalOutput(ifstream *defaultStream)
{
	char *localTokenList[20];
	int  numberOfTokens;
	int foundMatch;

	ifstream *currentStream;
    BETA_OUT<<"Current stream(readCommands)="<<(*defaultStream)<<endl;
	currentStream = defaultStream;

	while(2==2) {
		foundMatch=0;		
		currentStream = defaultStream;		
		if( getLineAndTokenize(currentStream,"exitReadOptionalOutput",localTokenList, numberOfTokens)==1) {
			goto _ExitChoose;
		}
	
		foundMatch=processOptionalOutput(localTokenList,currentStream, numberOfTokens); // decendent class first

		if(foundMatch ==0) {	exit(1);}
		
	} //end of while

	_ExitChoose:
	BETA_OUT << "Exiting readOptionalOutput"<<endl; 
	return true;
}
//===================================================================================
int BasicModel::processOptionalOutput(char** localTokenList,ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;

	char token[200];
	int  foundMatch=0;

strcpy(token,localTokenList[0]);
BETA_OUT<<"\nCommand processor = 'BasicModel::processOptionalOutput'          "<<token<<endl;

//========
	/* example:
	if(COMPARE(token,"concentration")==0){
		concentration_out=true;
		concentration_out_Step_Interval=1;
		if(numToken==2)
			concentration_out_Step_Interval=atoi(localTokenList[1]);
		OK 
	}
	*/

	if(foundMatch ==0){	
		BETA_OUT << "Command: " << token << endl;
		BETA_OUT << "Did not recognize command: " <<token<<endl;
		exit(1);
		//eventually make basicmodel the parent class for this function.
	}

	return foundMatch;
}
//======================================================================
void BasicModel::createAssembler()
{
	if (ParallelAssembly){

		if(parallelAssembly_ESA){
			assembler = new ESA_Assembler;
		}
		if(parallelAssembly_OFA){
			assembler = new OFA_Assembler;
		}
	}
	else{
		assembler = new SerialAssembler;
	}

}
//====================================================================

void BasicModel::setAssemblerPointers()
{
	assembler->setEquationsPointer(&equations);
	assembler->setFileManagerPointer(filemanager);
	assembler->setMeshPointer(mesh);
	assembler->setZeroVectorPointer(zeroVector);
	assembler->setMaterialListPointer(&materialList);
}
