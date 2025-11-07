#include "stdafx.h"

#include "ElasticityModel.hpp"
#include "utility/formWriter.hpp"
#include "math/gaussQuadrature/gaussPointList.hpp"
#include "mesh/MeshUtility.h"
#include "mesh/NodeSet.hpp"

#include "elements/3D/ElasticityElement3D.hpp"
//=========================================================================
#include "BCs/DisplacementLoad.hpp"
#include "BCs/NodeGroupDisplacementLoad.hpp"
#include "BCs/NodeGroupDisplacementLoadWithRotation.hpp"
#include "BCs/PlaneDisplacement.hpp"
//=========================================================================
BasicElement* getElementContainingNodes(bool linear, Node **n1);
void ExtractSubMesh(BasicMesh &mesh, istream  &is);
void ReadMPCFile(string filename, BasicMesh &mesh, Equations &equations);

extern double MPCTransformationTolerance;
//-------------------------------------------------------------------
CreateErrorHandler(ElasticityModel);
//--------------------------------------------------------------------
ElasticityModel::ElasticityModel(void)
{
//	OutputOptions=0;
//	globalstressout	= &cout;
//	globalstrainout	= &cout;
//	localstressout	= &cout;
//	localstrainout	= &cout;
//	quadstressout	= &cout;
//	quadstrainout	= &cout;
//	volumeAverageout = &cout;
//	GlobalForcesout = &cout;

	dispRequired=false;
	NodalLCS_stressRequired=false;
	NodalGCS_stressRequired=false;
	Quad_stressRequired=false;
	NodalLCS_strainRequired=false;
	NodalGCS_strainRequired=false;
	Quad_strainRequired=false;
	Volume_AveragesRequired=false;
	Global_ForcesRequired=false;
    Nodal_Forces_By_ElementRequired=false;
    strainEnergyRequired=false;
    QuadPointsRequired=false;
    QuadPointOrientaitonRequired=false;


}
//--------------------------------------------------------------------
ElasticityModel::~ElasticityModel(void)
{
}
//--------------------------------------------------------------------
int ElasticityModel::processCommand(char** localTokenList,ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;
	char newFileName[80];

	char token[200];
	int  foundMatch=0;

	strcpy(token,localTokenList[0]);
	BETA_OUT<<"\nCommand processor = 'ElasticityModel::processCommand'     "<<token<<endl;

	if( COMPARE(token,"openFile")==0 ) {  banner("openFile", BETA_OUT);  //This is for this pass only
		//Should be 3 tokens
		strcpy(newFileName,localTokenList[1]);
		currentStream = openFileForInput(newFileName);
		strcpy(token,localTokenList[2]);
		BETA_OUT<<"File for new stream = "<<newFileName<<endl;
		BETA_OUT<<"New stream  = "<<*currentStream<<endl;
		//OK //think about this
	}
	//----------- Insert new commands below this line ----------------------

	IS_COMMAND(ReadMultiPointConstraints)
	IS_COMMAND(CalculateStressAlongLine)

	IS_COMMAND(CalculateSubCellVolumeAverage)
		

	if( COMPARE(token,"ExtractSubMesh")==0 ){
		ExtractSubMesh(*mesh, *currentStream);
		OK 
	} 

	if( COMPARE(token,"ReadMPCFile")==0 ){
		string mpcfile(localTokenList[1]);
		ReadMPCFile(mpcfile, *mesh, equations);
		OK 
	} 

	if( COMPARE(token,"SetMPCTransformationTolerance")==0 ) {
		MPCTransformationTolerance=atof(localTokenList[1]);
		BETA_OUT	<< "NEW MPCTransformationTolerance = " << DOUBLE_FORMAT << MPCTransformationTolerance << endl;
		cout		<< "NEW MPCTransformationTolerance = " << DOUBLE_FORMAT << MPCTransformationTolerance << endl;
		OK
	}

	//----------- Insert new commands above this line ----------------------

	if(currentStream!=originalStream) {
		BETA_OUT<<"Close alternate input stream"<<endl;
		(*currentStream).close();
		delete currentStream;
		currentStream=originalStream;
		BETA_OUT<<"Revert to default stream "<< *originalStream
			<<" upon return to readCommands"<<endl;
	}

	if(foundMatch ==0){
		//BETA_OUT<<"\nNo match in 'ElasticityModel::processCommand' for command:  "<<token<<endl;
		//BETA_OUT<<"Use 'BasicModel::processCommand'\n\n";
        BETA_OUT<<"Process command by parent\n";
		foundMatch=BasicModel::processCommand(localTokenList,currentStream, numToken);
	}

	return foundMatch;
}
//======================================================================
void ElasticityModel::OptionalOutput(istream * inStream)
{
	EXIT_BETA("void ElasticityModel::OptionalOutput(istream * inStream) has been deprecated.\nUse 'ReadOptionalOutput' instead!");
/*
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
		else if(COMPARE(localTokenList[0],"localstress")==0)
			OutputOptions |= localstress;
		else if(COMPARE(localTokenList[0],"localstrain")==0)
			OutputOptions |= localstrain;
		else if(COMPARE(localTokenList[0],"globalstress")==0)
			OutputOptions |= globalstress;
		else if(COMPARE(localTokenList[0],"globalstrain")==0)
			OutputOptions |= globalstrain;
		else if(COMPARE(localTokenList[0],"quadstress")==0)
			OutputOptions |= quadstress;
		else if(COMPARE(localTokenList[0],"quadstrain")==0)
			OutputOptions |= quadstrain;
		else if(COMPARE(localTokenList[0],"volumeAverage")==0)
			OutputOptions |= volumeAverage;
		else if(COMPARE(localTokenList[0],"GlobalForces")==0)
			OutputOptions |= globalforces;
		else{
			FatalError("Incorrect Optional Output Specification!");
		}
	}//end of while

	openOutputFiles(filemanager, OutputOptions);

	//output to files
	outputToFiles(OutputOptions);

	//closing files
	closeOutputFiles(filemanager, OutputOptions);
*/
}
//==========================================================
/*
void ElasticityModel::openOutputFiles(FileManager *fm, int flags)
{
	if(flags == 0){
		BETA_OUT << "No Optional Output Parameters found!" << endl;
		return;
	}
	if(flags & localstress || flags & localstrain){  //BOTH stress and strain will be printed
		localstressout=fm->OpenOutputStream("stress");
		*localstressout << "stress" << endl;
		localstrainout=fm->OpenOutputStream("strain");
		*localstrainout << "stress" << endl;
	}
	if(flags & globalstress){
		globalstressout=fm->OpenOutputStream("globalstress");
		*globalstressout << "stress" << endl;
	}
	if(flags & globalstrain){
		globalstrainout=fm->OpenOutputStream("globalstrain");
		*globalstrainout << "stress" << endl;
	}
	if(flags & quadstress){
		quadstressout=fm->OpenOutputStream("quadstress");
	}
	if(flags & quadstrain){
		quadstrainout=fm->OpenOutputStream("quadstrain");
	}
	if(flags & volumeAverage){
		volumeAverageout=fm->OpenOutputStream("volumeAverage");
	}
	if(flags & globalforces){
		GlobalForcesout=fm->OpenOutputStream("GlobalForces");
	}
}
//=========================================================================
void ElasticityModel::closeOutputFiles(FileManager *fm, int flags)
{
	if(flags & localstress){
		fm->CloseOutputStream((ofstream*)localstressout);
	}
	if(flags & localstrain){
		fm->CloseOutputStream((ofstream*)localstrainout);
	}
	if(flags & globalstress){
		fm->CloseOutputStream((ofstream*)globalstressout);
	}
	if(flags & globalstrain){
		fm->CloseOutputStream((ofstream*)globalstrainout);
	}
	if(flags & quadstress){
		fm->CloseOutputStream((ofstream*)quadstressout);
	}
	if(flags & quadstrain){
		fm->CloseOutputStream((ofstream*)quadstrainout);
	}
	if(flags & volumeAverage){
		fm->CloseOutputStream((ofstream*)volumeAverageout);
	}
	if(flags & globalforces){
		fm->CloseOutputStream((ofstream*)GlobalForcesout);
	}

	if(flags == 0){
		BETA_OUT << "No Optional Output Parameters found!" << endl;
	}
}
//======================================================================
void ElasticityModel::outputToFiles(int flags)
{
	int i;
	double * theSolution=0;

	int totalNumIPs;
	GaussPointList  gaussPointList;

	theSolution=equations.getSolution();

	ElasticityElement3D *e=0;

	equations.zeroResultantVector();	
	for(i=0; i<mesh->numElements; i++)
	{
		e=(ElasticityElement3D*)(&mesh->element[i]);

		int group = e->getMaterialNumber();		
		e->update("F", theSolution);
		e->getQuadraturePoints(gaussPointList);
		totalNumIPs=gaussPointList.totalNumIPs;

		if(flags & localstress){
			e->extrapolateStresses(totalNumIPs, *localstressout, *localstrainout);
		}

		if(flags & globalstress){
			(*globalstressout) << e->elementNumber << " " << e->getMaterialNumber() << endl;
			e->extrapolateStresses(totalNumIPs, false, *globalstressout);
		}

		if(flags & globalstrain){
			(*globalstrainout) << e->elementNumber << " " << e->getMaterialNumber() << endl;
			e->extrapolateStrains(totalNumIPs, false, *globalstrainout);
		}
	
		if(flags & quadstress){
			e->printQuadStresses(totalNumIPs, *(quadstressout) );
		}

		if(flags & quadstrain){
			e->printQuadStrains(totalNumIPs, *(quadstrainout) );
		}
	}//end of loop on elements

	if(flags & volumeAverage)
	{
		BETA_OUT << "Printing Volume Averages..." << endl;
		mesh->element[0].initializeSummary(); // since getglobal forces is not called, i have to initialize the summary here 
		outVolumeStrConstituents(theSolution,0, 0, 0, true);
		mesh->element[0].outputSummary(); //no sense in printing summary since it just got initialized
	}

	if(flags & globalforces)
	{
		BETA_OUT << "Printing GlobalForces..." << endl;
		int    dofForMaxResidual;
		double maxResidual;
		double * mappedResidualVector = new double[totalNumDof];
	
		getGlobalForces(internalForce, currentSolution,mappedResidualVector, maxResidual,dofForMaxResidual );
		mesh->printNodalData(internalForce, "displacements", GlobalForcesout);
		delete [] mappedResidualVector;
	}

}
*/
//=========================================================================
bool ElasticityModel::createLoads(char* name, istream *inStream)
{
	bool exitLoop;
	Load *newLoad=0;

	IS_LOAD(DisplacementLoad)	
	IS_LOAD(NodeGroupDisplacementLoad)
	IS_LOAD(NodeGroupDisplacementLoadWithRotation)	
	IS_LOAD(PlaneDisplacement)	
	else
		return BasicModel::createLoads(name, inStream);
}
//======================================================================
//fix this function - the output is different from the alpha version FIXXXXXXXX
void ElasticityModel::outVolumeStrConstituents(double *disp,int matVol, int numCpt, double ** cpt, ostream *volumeAverageout, bool byConstituent=false)
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
	ElasticityElement3D *e=0;

	for(i=0;i<100;i++) {matID[i]=0;}

	if (byConstituent) {
		for(i=0;i<numElements;i++) {
			e=(ElasticityElement3D *)&mesh->element[i];

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

	int eListCount=0;
	if(!byConstituent){
		*(volumeAverageout) << "List of elements in the subcells"<< endl; 
	}

    for(i=0;i<numElements;i++) 
	{
		e=(ElasticityElement3D *)&mesh->element[i];
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
		if(!byConstituent && inRegion){
			*(volumeAverageout) << i <<"\t";
			eListCount++;
			if(eListCount % 10 == 0) *(volumeAverageout) << endl;
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

	if(!byConstituent){
		*(volumeAverageout) << "\nEnd of Element List\n" << endl;
	}

	ElasticMaterial * matPoint;
	e=(ElasticityElement3D*)(&mesh->element[0]); // Fix this... what should be done ??

	for (i=0; i<=numMats; i++) 
	{
		*(volumeAverageout) <<matID[i]<<"\t"<<Totalvolume[i]<<"\t";  
		matPoint = (ElasticMaterial * ) &(materialList[matID[i]]);				
		//if (matPoint->filemanager) //checking if it is an acutal material
		//{
		//	*(volumeAverageout) << matPoint->xtNumAngles<<"\t"<<matPoint->xtAngles[0]<<"\t"<<matPoint->xtAxis[0];
		//}
		//else 
            *(volumeAverageout) << -1<<"\t"<<0<<"\t"<<0;
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

	//JV081007
	//updating the element's vol avg data variable. 
	//so that element[0].outputSummary() prints the volume averaging info in 
	//the output file as well. this is just a convenince fix.
	//should decide a more formal procedure to do this later.
	ElementWorkspace* bag=e->getElementWorkspace();
	bag->totalVolume = Totalvolume[0];
	for(i=0; i<6; i++) {
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
#undef maxStresses
#undef FORMAT
}
//======================================================================
void ElasticityModel::CalculateEnergy(double *force, double *disp)
{
#define maxStresses 6
#define FORMAT setw(15)<<setprecision(6)
    double xEnergy=0., yEnergy=0., zEnergy=0., Energy =0.;
    double *volAvgStress, *volAvgStrain, volume, *SED, SEDall=0.;
    double Totalvolume=0.;
    double *TotalvolAvgStress, *TotalvolAvgStrain;
    int i, is, numDofAtNode;

    volAvgStress = new double [maxStresses];
    volAvgStrain = new double [maxStresses];
    SED = new double [3];
    TotalvolAvgStress = new double [maxStresses];
    TotalvolAvgStrain = new double [maxStresses];
    for(i=0;i<maxStresses;i++){
        TotalvolAvgStress[i]=0.0;
        TotalvolAvgStrain[i]=0.0;
    }
    // Calculate Volume
	ElasticityElement3D *e=0;
    for(i=0;i<mesh->numElements;i++) {
		e=(ElasticityElement3D*)(&mesh->element[i]);
        
        e->calculateDofList();
        volume=0.;
        e->getVolAvgValues(disp,volAvgStress,volAvgStrain, SED, &volume);      
        Totalvolume       += volume;
        SEDall += SED[0];
        for(is=0; is<maxStresses; is++){
           TotalvolAvgStress[is] += volAvgStress[is];
           TotalvolAvgStrain[is] += volAvgStrain[is];
        }
    }
    BETA_OUT<<endl;    
    BETA_OUT<<"-----------------------------------------------------------"<<endl;
    double value0;
    #define FF1 setw(7)<<setprecision(4)<<setiosflags(ios::fixed)

    value0=SEDall/Totalvolume;
     BETA_OUT<<"Strain Energy Density of the model  = "<<value0<<endl;

    
    BETA_OUT<<"Volume Averaged Strain and Stress (11,22,33,12,23,13 Order)"<<endl;
    for(is=0;is<maxStresses;is++)
    {   BETA_OUT<<FORMAT<<TotalvolAvgStrain[is]/Totalvolume;  } BETA_OUT<<endl;
    for(is=0;is<maxStresses;is++)
    {   BETA_OUT<<FORMAT<<TotalvolAvgStress[is]/Totalvolume;  } BETA_OUT<<endl;
    BETA_OUT<<"Volume of Model = "<<Totalvolume<<endl;

	

    cout<<"Calculate Energy=1/2*Force*disp"<<endl;  
    BETA_OUT<<"Calculate Energy=1/2*Force*disp"<<endl;   
    for(i=0;i<mesh->numNodes;i++)
    {           
        BETA_OUT<<setiosflags(ios::scientific);
        numDofAtNode =mesh->node[i].getNumDof(); 
        xEnergy += force[mesh->node[i].getFirstDof()  ]*disp[mesh->node[i].getFirstDof()  ];
        yEnergy += force[mesh->node[i].getFirstDof()+1]*disp[mesh->node[i].getFirstDof()+1];
        if(numDofAtNode==3){
        zEnergy += force[mesh->node[i].getFirstDof()+2]*disp[mesh->node[i].getFirstDof()+2];
        }
        else {zEnergy=0.;}  
    }
    
	Energy = 0.5 * (xEnergy + yEnergy + zEnergy);
    cout<<"EEngergy = "<<FORMAT<<Energy<<endl;
    BETA_OUT<<"EEngergy = "<<FORMAT<<Energy<<endl;
    BETA_OUT<<"StrainEnergyDensity/Volume = "<<FORMAT<<Energy/Totalvolume<<endl;
    BETA_OUT<<"From the equation: (SE/volume) = 1/2 * modulus * strain^2 "<<endl;
    BETA_OUT<<"When strain is 1%,  Modulus = "<<FORMAT<<Energy/Totalvolume*2./0.01/0.01<<endl;
    BETA_OUT<<"-----------------------------------------------------------"<<endl;

    delete [] volAvgStress; //xt062901
    delete [] volAvgStrain; //xt062901
    delete [] SED; //xt062901
    delete [] TotalvolAvgStress; //xt062901
    delete [] TotalvolAvgStrain; //xt062901

}
#undef FORMAT
//======================================================================
void ElasticityModel::ProcessSetElementProperty(string command, istream * inStream)
{

	int intCommand=0 ;
	enum {
		SET_ELEMENT_MOISTURE, //xtang 122000
	};

	ChangeToUpper(command)
	if(command=="SETELEMENTMOISTURE")
		{intCommand=SET_ELEMENT_MOISTURE;}
	else {
		return BasicModel::ProcessSetElementProperty(command, inStream);
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
	    exit(1);}				


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
//=========================================================================================
void ElasticityModel::CalculateStressAlongLine(istream * inStream)
{
	ifstream is;
	string outfilename;
	char StressFilename[256];
	Node Point1,Point2;

//*********************************************
//******Reading Commands From Script***********
//*********************************************
	char *localTokenList[20];
	int  numberOfTokens;
	while(true){
		if( getLineAndTokenize(inStream,"ExitCalculateStressAlongLine",localTokenList, numberOfTokens)==1) {
			break;
		}

		//read the points defining the line
		//Point1:
		if( COMPARE(localTokenList[0],"Point1")==0 ){
			if(numberOfTokens == 1){
				BETA_OUT << "Point1 definition incorrect - need to specify coordinates" << endl;
				exit(1);
			}
			if(numberOfTokens > 1) 
				Point1.x=atof(localTokenList[1]);
			if(numberOfTokens > 2) 
				Point1.y=atof(localTokenList[2]);
			if(numberOfTokens > 3) 
				Point1.z=atof(localTokenList[3]);
			continue;
		}
		//Point2:
		if( COMPARE(localTokenList[0],"Point2")==0 ){
			if(numberOfTokens == 1){
				BETA_OUT << "Point2 definition incorrect - need to specify coordinates" << endl;
				exit(1);
			}
			if(numberOfTokens > 1) 
				Point2.x=atof(localTokenList[1]);
			if(numberOfTokens > 2) 
				Point2.y=atof(localTokenList[2]);
			if(numberOfTokens > 3) 
				Point2.z=atof(localTokenList[3]);
			continue;
		}

		if( COMPARE(localTokenList[0],"StressFile")==0 ){
			strcpy(StressFilename, localTokenList[1]);
			continue;
		}
		if( COMPARE(localTokenList[0],"OutputFile")==0 ){
			outfilename=localTokenList[1];
			//strcpy(outfilename, localTokenList[1]);
			continue;
		}
		//cannot find command:
		EXIT_BETA("Incorrect input for CalculateStressAlongLine!");
	}
//*************************************************
//*******End Reading Commands From Script**********
//*************************************************

//add a check to see if you have got all the needed information from the script file

    //get the list of nodes on the line defined by point1 and point2
    NodeSet nlist(mesh->getNodeSearchTolerance());
    double P1[] = {Point1.x, Point1.y, Point1.z};
    double P2[] = {Point2.x, Point2.y, Point2.z};
    nlist.addOrderedNodesOnLineSegment(mesh,P1,P2);
    //printing the sorted nodes for testing
    //for (int i=0;i<nlist.getNumNodes();i++){BETA_OUT << nlist[i]->nodeNum << endl;}

	//this determines which elements are connected to a node
	mesh->setupNodeDoc();

	int numMaterials=materialList.getNumElements();
	Array<ofstream> fileStreams;
	fileStreams.reserve(numMaterials);

	char suffix[256];
	sprintf(suffix,".txt");
	string fn=outfilename + suffix;
	fileStreams.add(*(filemanager->OpenOutputStream(fn)));
	fileStreams[0] << "x \t y \t z \t stress " << endl;

	for(int i=1; i<numMaterials; i++){
		sprintf(suffix,"_%d.txt",i);
		string fn=outfilename + suffix;
		fileStreams.add(*(filemanager->OpenOutputStream(fn)));
		fileStreams[i] << "x \t y \t z \t stress " << endl;
	}

	Node locatedNode;
	Node *n;

	int numStresses=6;
	int numLines=0;
	//loop to calc how many lines in the output file
	for(int i=0; i<nlist.getNumNodes(); i++){
		n=nlist[i];
		for(int j=0;j<n->numElementsAtNode;j++)numLines++;
	}
	//allocating memory for storing stress values for each line
	double **data;
	data= new double *[numLines];
	for(int i=0;i<numLines;i++){
		data[i]=new double[numStresses];
	}

	//reading stress file
	mesh->contourDataSet.doAveraging=true;
	mesh->contourDataSet.doNormalize=false;
	mesh->readContourData(StressFilename);

	int linectr=0;
	int stressIndex=0;
	for(stressIndex=0; stressIndex<numStresses; stressIndex++){//loop over stress components
		if(!mesh->contourDataSet.switchColumn(stressIndex+1)){//column numbering starts from 1
			cout << "error in switching columns" << endl;
			exit(1);
		}
		linectr=0;
		for(int i=0; i<nlist.getNumNodes(); i++){
			n=nlist[i];
			for(int j=0;j<n->numElementsAtNode;j++){
				BasicElement *e=n->eList[j];
				double value = e->nodalValues[getIndex(n, e->node)];
				data[linectr][stressIndex]=value;
				linectr++;
			}
		}
	}

	//outputting to file
	linectr=0;
	fileStreams[0].setf(ios::scientific);
	for(int i=0; i<nlist.getNumNodes(); i++){
		n=nlist[i];
		for(int j=0;j<n->numElementsAtNode;j++){
			BasicElement *e=n->eList[j];
			fileStreams[0] << DOUBLE_FORMAT << n->x << "  " << DOUBLE_FORMAT<< n->y << "  " << DOUBLE_FORMAT<< n->z << "  ";
			fileStreams[e->getMaterialNumber()] << DOUBLE_FORMAT << n->x << "  " << DOUBLE_FORMAT<< n->y << "  " << DOUBLE_FORMAT<< n->z << "  ";
			for(int k=0;k<numStresses; k++){
				fileStreams[0] << DOUBLE_FORMAT<< data[linectr][k] << "  ";
				fileStreams[e->getMaterialNumber()] << DOUBLE_FORMAT<< data[linectr][k] << "  ";
			}
			fileStreams[0] << setw(8)<< e->getMaterialNumber() << endl;	
			fileStreams[e->getMaterialNumber()] << setw(8)<< e->getMaterialNumber() << endl;	
			linectr++;
		}
	}

//_Exit_CalculateStressAlongLine:
	for(int i=0; i<numMaterials; i++){
		filemanager->CloseOutputStream(&fileStreams[i]);
	}
	fileStreams.clear();

	for(int i=0;i<numLines;i++) delete [] data[i];
	delete [] data;

} // end of CalculateStressAlongLine
//===========================================================
void ElasticityModel::CalculateSubCellVolumeAverage(istream * inStream)
{
	ifstream is;
	string outfilename;
	Node Point1,Point2;
	int numSubCells=0;
	int iSubCell=0;

	double ** cpt=0;
//************************************************************************
//***********************Reading Commands From Script*********************
//************************************************************************
	char *localTokenList[20];
	int  numberOfTokens;
	while(true){
		if( getLineAndTokenize(inStream,"ExitCalculateSubCellVolumeAverage",localTokenList, numberOfTokens)==1) {
			break;
		}

		if( COMPARE(localTokenList[0],"NumberOfSubCells")==0 ){
			numSubCells=atoi(localTokenList[1]);
			cpt=new double*[numSubCells];
			for(int i=0; i<numSubCells; i++){
				cpt[i]=new double [6];
			}
			iSubCell=0;
			continue;
		}
		//read the points defining the subcells
		if( COMPARE(localTokenList[0],"SubCell")==0 ){
			if(numberOfTokens != 7){
				BETA_OUT << "Subcell definition incorrect - need to specify coordinates for two points" << endl;
				exit(1);
			}
			cpt[iSubCell][0]=atof(localTokenList[1]);
			cpt[iSubCell][1]=atof(localTokenList[2]);
			cpt[iSubCell][2]=atof(localTokenList[3]);
			cpt[iSubCell][3]=atof(localTokenList[4]);
			cpt[iSubCell][4]=atof(localTokenList[5]);
			cpt[iSubCell][5]=atof(localTokenList[6]);
			iSubCell++;
			continue;
		}
		if( COMPARE(localTokenList[0],"OutputFile")==0 ){
			outfilename=localTokenList[1];
			//strcpy(outfilename, localTokenList[1]);
			continue;
		}
		//cannot find command:
		BETA_OUT << "Incorrect input for CalculateSubCellVolumeAverage!" << endl;
		exit(1);
	}
//************************************************************************
//***********************End Reading Commands From Script*****************
//************************************************************************
	ofstream *volumeAverageout=filemanager->OpenOutputStream(outfilename);
	outVolumeStrConstituents(currentSolution,0, 1, cpt,volumeAverageout, false);
	filemanager->CloseOutputStream(volumeAverageout);

}
//===================================================================================
int ElasticityModel::processOptionalOutput(char** localTokenList,ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;

	char token[200];
	int  foundMatch=0;

strcpy(token,localTokenList[0]);
BETA_OUT<<"\nCommand processor = 'ElasticityModel::processOptionalOutput'          "<<token<<endl;

//========
	if( COMPARE(token,"Displacements")==0 )
	{
		dispRequired=true;
		OK 
	}

	if( COMPARE(token,"StrainEnergy")==0 )
	{
		strainEnergyRequired=true;
		OK
	}

	if( COMPARE(token,"LCS_stress")==0 ||
		COMPARE(token,"stress")==0 )
	{
		NodalLCS_stressRequired=true;
		OK 
	}

	if( COMPARE(token,"GCS_stress")==0 )
	{
		NodalGCS_stressRequired=true;
		OK 
	}

	if( COMPARE(token,"Quad_stress")==0 )
	{
		Quad_stressRequired=true;
		OK 
	}

	if( COMPARE(token,"LCS_strain")==0 ||
		COMPARE(token,"strain")==0 )
	{
		NodalLCS_strainRequired=true;
		OK 
	}

	if( COMPARE(token,"GCS_strain")==0 )
	{
		NodalGCS_strainRequired=true;
		OK 
	}

	if( COMPARE(token,"Quad_strain")==0 )
	{
		Quad_strainRequired=true;
		OK 
	}

	if( COMPARE(token,"volumeAverage")==0 )
	{
		Volume_AveragesRequired=true;
		OK 
	}

	if( COMPARE(token,"GlobalForces")==0 )
	{
		Global_ForcesRequired=true;
		OK 
	}

    if( COMPARE(token,"NodalForcesByElement")==0 )
	{
		Nodal_Forces_By_ElementRequired=true;
		OK 
	}

    if( COMPARE(token,"QuadPoints")==0 )
	{
		QuadPointsRequired=true;
		OK 
	}

    if( COMPARE(token,"QuadPointOrientation")==0 )
	{
		QuadPointOrientaitonRequired=true;
		OK 
	}

	if(foundMatch ==0){	
		foundMatch=BasicModel::processOptionalOutput(localTokenList,currentStream, numToken);
	}

	return foundMatch;
}
//==========================================================
void ElasticityModel::WriteOptionalOutput()
{
//	SET_SCIENTIFIC(BETA_OUT);

	ofstream *ost;
	int i;
	int numElements=mesh->numElements;
	ElasticityElement3D *e=0;

//////////////////////////////////////////////////////////////////
	if(dispRequired)	{
		string filename= "disp" ;
		ost = filemanager->OpenOutputStream(filename);
		mesh->printNodalData(currentSolution, "displacements", ost );
		filemanager->CloseOutputStream(ost);
	}


//checking if elemental output is needed
	bool doElementLoop=false;
	bool printNodalLCS_stress=false;
	bool printNodalGCS_stress=false;
	bool printQuad_stress=false;
	bool printNodalLCS_strain=false;
	bool printNodalGCS_strain=false;
	bool printQuad_strain=false;
    bool printQuadPoints=false;
    bool printQuadPointOrientation=false;
    bool printNodalForcesByElement=false;

	ofstream *NodalLCS_stress_ost;
	ofstream *NodalGCS_stress_ost;
	ofstream *Quad_stress_ost;
	ofstream *NodalLCS_strain_ost;
	ofstream *NodalGCS_strain_ost;
	ofstream *Quad_strain_ost;
    ofstream *QuadPoints_ost;
    ofstream *QuadPointOrientation_ost;
    ofstream *NodalForcesByElement_ost;

	if(NodalLCS_stressRequired){
		doElementLoop=true;
		printNodalLCS_stress=true;
		string filename = "stress" ;
		NodalLCS_stress_ost = filemanager->OpenOutputStream(filename);
		(*NodalLCS_stress_ost)<<"stress"<<endl;
	}
	if(NodalGCS_stressRequired){
		doElementLoop=true;
		printNodalGCS_stress=true;
		string filename = "globalstress" ;
		NodalGCS_stress_ost = filemanager->OpenOutputStream(filename);
		(*NodalGCS_stress_ost)<<"stress"<<endl;
	}
	if(Quad_stressRequired){
		doElementLoop=true;
		printQuad_stress=true;
		string filename = "Quad.stress" ;
		Quad_stress_ost = filemanager->OpenOutputStream(filename);
		SET_SCIENTIFIC(*Quad_stress_ost);
	}
	if(NodalLCS_strainRequired){
		doElementLoop=true;
		printNodalLCS_strain=true;
		string filename = "strain" ;
		NodalLCS_strain_ost = filemanager->OpenOutputStream(filename);
		(*NodalLCS_strain_ost)<<"stress"<<endl;
	}
	if(NodalGCS_strainRequired){
		doElementLoop=true;
		printNodalGCS_strain=true;
		string filename = "globalstrain" ;
		NodalGCS_strain_ost = filemanager->OpenOutputStream(filename);
		(*NodalGCS_strain_ost)<<"stress"<<endl;
	}
	if(Quad_strainRequired){
		doElementLoop=true;
		printQuad_strain=true;
		string filename = "Quad.strain" ;
		Quad_strain_ost = filemanager->OpenOutputStream(filename);
		SET_SCIENTIFIC(*Quad_strain_ost);
	}
    if(QuadPointsRequired){
        doElementLoop=true;
        printQuadPoints=true;
        string filename = "QuadPoints" ;
        QuadPoints_ost = filemanager->OpenOutputStream(filename);
    }
    if(QuadPointOrientaitonRequired){
        doElementLoop=true;
        printQuadPointOrientation=true;
        string filename = "QuadPointOrientation" ;
        QuadPointOrientation_ost = filemanager->OpenOutputStream(filename);
    }
    if(Nodal_Forces_By_ElementRequired){
        doElementLoop=true;
        printNodalForcesByElement=true;
        string filename = "NodalForcesByElement" ;
        NodalForcesByElement_ost = filemanager->OpenOutputStream(filename);
    }
	if (doElementLoop==false)
		goto SKIP_UPDATE_F_LOOP;


//writing elemental output
	GaussPointList  gaussPointList;
	int totalNumIPs;
	equations.zeroResultantVector();	
	for(i=0; i<numElements; i++)
	{
		e=(ElasticityElement3D*)(&mesh->element[i]);
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
		if(printNodalGCS_strain){
			e->printNodalStrainsInGCS(*NodalGCS_strain_ost);
		}
		if(printNodalLCS_strain){
			e->printNodalStrainsInLCS(*NodalLCS_strain_ost);
		}
		if(printQuad_strain){
			e->printQuadStrains(totalNumIPs, *Quad_strain_ost);
		}
        if(printQuadPoints){
            (*QuadPoints_ost) << e->elementNumber << "\t" << group << "\t" << e->getTotalNumberOfIPs() << endl;
            e->printQuadPoints(totalNumIPs, *QuadPoints_ost);
        }
        if(printQuadPointOrientation){
            e->printQuadOrientations(*QuadPointOrientation_ost);
        }
        if(printNodalForcesByElement){
            e->printNodalForcesByElement(*NodalForcesByElement_ost);
        }
	}
	if(printNodalGCS_stress)	filemanager->CloseOutputStream(NodalGCS_stress_ost);
	if(printNodalLCS_stress)	filemanager->CloseOutputStream(NodalLCS_stress_ost);
	if(printQuad_stress)		filemanager->CloseOutputStream(Quad_stress_ost);
	if(printNodalGCS_strain)	filemanager->CloseOutputStream(NodalGCS_strain_ost);
	if(printNodalLCS_strain)	filemanager->CloseOutputStream(NodalLCS_strain_ost);
	if(printQuad_strain)		filemanager->CloseOutputStream(Quad_strain_ost);
    if(printQuadPoints)         filemanager->CloseOutputStream(QuadPoints_ost);
    if(printQuadPointOrientation) filemanager->CloseOutputStream(QuadPointOrientation_ost);
    if(printNodalForcesByElement) filemanager->CloseOutputStream(NodalForcesByElement_ost);

SKIP_UPDATE_F_LOOP:
	if(Volume_AveragesRequired){
		string filename = "volumeAverage" ;
		ost = filemanager->OpenOutputStream(filename);
		BETA_OUT << "Printing Volume Averages..." << endl;
		mesh->element[0].initializeSummary(); // since getglobal forces is not called, i have to initialize the summary here 
		outVolumeStrConstituents(currentSolution,0, 0, 0, ost, true);
		mesh->element[0].outputSummary(); //no sense in printing summary since it just got initialized
		filemanager->CloseOutputStream(ost);
	}

	if(strainEnergyRequired){
		CalculateEnergy(internalForce,currentSolution);
	}

	if(Global_ForcesRequired){
		string filename = "GlobalForces" ;
		ost = filemanager->OpenOutputStream(filename);
		BETA_OUT << "Printing GlobalForces..." << endl;
		int    dofForMaxResidual;
		double maxResidual;
		double * mappedResidualVector = new double[totalNumDof];
	
		getGlobalForces(internalForce, currentSolution,mappedResidualVector, maxResidual,dofForMaxResidual );
		mesh->printNodalData(internalForce, "displacements", ost);
		delete [] mappedResidualVector;
		filemanager->CloseOutputStream(ost);
	}

}
//==========================================================
