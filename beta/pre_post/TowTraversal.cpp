#include "stdafx.h"

#include "models/BasicModel.hpp"
#include "TowTraversalUtility.hpp"
#include "mesh/MeshUtility.h"

#define TOWPATH_TOLERANCE 1e-4

#include <cmath>
#include <string>
#include <sstream>

bool readMaterialGroupsFromFile(const char *name, BasicMesh* mesh);
Node calculateCentroid(NodeGroup * nodeList);

//=============================================================
#include "models/ElasticityModel.hpp"
//====================================================================
//========================================================
#include "elements/3D/ElasticityElement3D.hpp"
//fix this ! why does basic model need elasticity element ???
//========================================================

//==========================================================================
void TowTraversalUtility::setFileManager(FileManager *fm)
{
	ElasticityModel::setFileManager(fm);
	towpath.setFileManager(fm);
}
//==========================================================================
int TowTraversalUtility::processCommand(charlist localTokenList, ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;
	
	char newFileName[80];

	char token[200];
	int  foundMatch=0;

strcpy(token,localTokenList[0]);
BETA_OUT<<"\nCommand processor = 'TowTraversalUtility::processCommand'          "<<token<<endl;
	int i;

	if( COMPARE(token,"openFile")==0 ) {   //This is for this pass only
		//Should be 3 tokens
		strcpy(newFileName,localTokenList[1]);
		currentStream = openFileForInput(newFileName);
		strcpy(token,localTokenList[2]);
		BETA_OUT<<"Input stream for command ***"<<token<<"*** = "<<newFileName<<endl;
		BETA_OUT<<"Stream(openFile) = "<<currentStream<<endl;
		BETA_OUT<<"Stream(openFile) = "<<(*currentStream)<<endl;  
		//OK //think about this
	}


	//========
	if( COMPARE(token,"ReadElematFile")==0 ) {
		strcpy(ElematFile,localTokenList[1]);
		readMaterialGroupsFromFile(ElematFile, mesh);
		OK
	}

	if( COMPARE(token,"ElementsOnTowSection")==0 ) {
		int elementNumber;
		*currentStream >> numElementsOnSection;
		for(i=0;i<numElementsOnSection;i++){
			*currentStream >> elementNumber;
			elementList.add(mesh->element[elementNumber]);
		}
		OK
	}

	if( COMPARE(token,"numberOfElementsOnTowSection")==0 ) {
		numElementsOnSection=atoi(localTokenList[1]);
		OK
	}

	if( COMPARE(token,"readMasterNode")==0 ) {
		MasterNode=atoi(localTokenList[1]);
		OK
	}

	if( COMPARE(token,"SectionOrientation")==0 ) {
		if( COMPARE(localTokenList[1],"Vertical")==0 ) {
			SectionOrientation=1;
		}
		OK
	}

	if( COMPARE(token,"TextileType")==0 ) {
		if( COMPARE(localTokenList[1],"PlainWeave")==0 ) {
			TextileType=1;
		}
		OK
	}

	if( COMPARE(token,"DefineTowPath")==0 ) {
		towpath.readCommands(currentStream);
		OK
	}


	if( COMPARE(token,"FindAllCrossSections")==0 ) {
		FindAllCrossSections();
		OK
	}

	if( COMPARE(token,"getGlobalForcesAtCrossSections")==0 ) {
		getGlobalForcesAtCrossSections(currentStream);
		OK
	}

	if(currentStream!=originalStream) {
		BETA_OUT<<"Close alternate input stream"<<endl;
		(*currentStream).close();
		delete currentStream;
		currentStream=originalStream;
		BETA_OUT<<"Automatically reverts to default stream "<<currentStream
			<<" upon return to readCommands"<<endl;
	}

	if(foundMatch ==0){	
		foundMatch=ElasticityModel::processCommand(localTokenList,currentStream,numToken);
	}

	return foundMatch;
}

void TowTraversalUtility::FindAllCrossSections()
{
	cout << "started FindAllCrossSections()" << endl;
	int i,j,k;
	CrossSectionInfo *csi;
	BasicElement *e;

	int numE;
	
	mesh->calculateElementCentroid();

	//finding all elements in the cross-sections 
	Array<BasicElement> *eList1;
	eList1=towpath.findElementsUsingMaterialGroup(mesh);

	int direction=0; //hardcoded for now !!!
	numE=1;
	Node *ec;
	double centroidCoord,tempCoord;
	int nnpe;
	int numNodesInNodeList;
	e=&(*eList1)[0];
	nnpe=e->numNodesPerElement;
	ec=&(e->centroid);
	centroidCoord=(*ec)(direction);

	for (i=1;i<eList1->getNumElements();i++){
		e=&(*eList1)[i];
		ec=&(e->centroid);
		tempCoord=(*ec)(direction);
		if( fabs(tempCoord-centroidCoord) <= TOWPATH_TOLERANCE)numE++;
		else break;
	}

	numElementsOnSection=numE;
	cout << "\nNumber of elements in single layer: " << numElementsOnSection << endl;

	int numCS=eList1->getNumElements()/numE;
	cout << "Number of layers: " << numCS << endl;
	int counter=0;

	//storing the cross section information
	CrossSectionList.setCapacity(numCS);
	for (i=0;i<numCS;i++){
		csi=new CrossSectionInfo;
		csi->elementList.setCapacity(numE);
		csi->nodeList.setCapacity(numE*nnpe);
		for (j=0;j<numE;j++){
			e=&(*eList1)[counter];
			csi->elementList.add(*e);
			for(k=0;k<e->numNodesPerElement;k++){
				if(csi->nodeList.isInByAddress(&(e->node[k])) == -1)
					csi->nodeList.add(e->node[k]);
			}
			counter++;
		}
		CrossSectionList.add(*csi);
	}

	//sort the node list in each csi
	int numElementsInList;
	Node *n1,*n2;
	for (k=0;k<numCS;k++){
		csi=&CrossSectionList[k];
		numElementsInList=csi->nodeList.getNumElements();
		for(i=0;i<numElementsInList-1;i++){
			for(j=i+1;j<numElementsInList;j++){
				n1=&csi->nodeList[i];
				n2=&csi->nodeList[j];
				if ( (*n1)(direction) > (*n2)(direction) ) 
					csi->nodeList.swap(i,j);

			}
		}
	}

	//for getting an estimate of the capacity needed for the individual layers
	numNodesInNodeList=CrossSectionList[0].nodeList.getNumElements();

	//split the nodes into the 3 layers in the csi
	int ctr=0;
	Node *prevn;
	for (k=0;k<numCS;k++){
		ctr=0;
		csi=&CrossSectionList[k];
		csi->nList[0].setCapacity(numNodesInNodeList);
		csi->nList[1].setCapacity(numNodesInNodeList);
		csi->nList[2].setCapacity(numNodesInNodeList);
		numElementsInList=csi->nodeList.getNumElements();
		for(i=0;i<numElementsInList;i++){
			n1=&csi->nodeList[i];
			if ( i!=0 && ((*n1)(direction) > (*prevn)(direction))  ) {
				ctr++;
				prevn=n1;
			}
			if(i==0) prevn=n1;
			(csi->nList[ctr]).add(*n1);
		}
	}
	
	//calculate the centroid for each cross section layer
	Node centroid;
	for (k=0;k<numCS;k++){
		csi=&CrossSectionList[k];
		centroid=calculateCentroid(&(csi->nList[0]));
		csi->centroid[0]=centroid;
		centroid=calculateCentroid(&(csi->nList[1]));
		csi->centroid[1]=centroid;
		centroid=calculateCentroid(&(csi->nList[2]));
		csi->centroid[2]=centroid;
	}

	//printing Cross Section information
	for (i=0;i<CrossSectionList.getNumElements();i++){
		csi=&CrossSectionList[i];
		BETA_OUT << "\nCrossSection Number:" << i << endl;
		BETA_OUT << "\nElements:" << endl;
		for (j=0;j<csi->elementList.getNumElements();j++){
			BETA_OUT << csi->elementList[j].elementNumber << endl;
		}
		BETA_OUT << "\nNodes:" << endl;
		for (j=0;j<csi->nodeList.getNumElements();j++){
			BETA_OUT << csi->nodeList[j].nodeNum << " ";
			if((j+1)%10 == 0) BETA_OUT << endl;
		}
	}
	BETA_OUT << endl;

	eList1->clear();
	delete eList1;
	cout << "finished FindAllCrossSections()" << endl;

}
//=========================================================
void TowTraversalUtility::getGlobalForcesAtCrossSections(istream *currentStream)
{
	char token[256];

	(*currentStream) >> token;

	if( COMPARE(token,"normal")==0 ) {   
		getGlobalForcesAtCrossSectionsNormal();			
	}
	else if( COMPARE(token,"splitup")==0 ) {   
		getGlobalForcesAtCrossSectionsSplitup();			
	}

}
//=========================================================
void TowTraversalUtility::getGlobalForcesAtCrossSectionsSplitup()
{
	cout << "started getGlobalForcesAtCrossSectionsSplitup()" << endl;

//	allocateForAnalysis(false);//not allocating equations //JV030409: did not allocate eqns in hs4rae, but in beta currently we have to
	allocateForAnalysis(true);
	setSolutionVectorFromMeshDisplacements(mesh);
	//setMaterial();
	

	Node *n;
	int i,j,k,p;
	int numNodes;

	CrossSectionInfo *csi, *csi2;
	double *globalForces=internalForce;


	double fx,fy,fz;
	double mxZ, mxY, myX, myZ, mzX, mzY;
	int firstDof;
	int direction=0; // hardcoded for now

//	ofstream os("ForceMoments_SplitUp.xls");
	stringstream number;
	number << towpath.MaterialGroup;
	string filename="ForceMoments_SplitUp_material_" + number.str() + ".xls";
	ofstream *ost = filemanager->OpenOutputStream(filename.c_str());
	ofstream &os=*ost;


//	os << "x\tfx\tfy\tfz\tmxZ\tmxY\tmyX(=0)\tmyZ(twist)\tmzX(=0)\tmzY(twist)\tMtwist=MyZ+MzY\n" << endl;
//	os << "x\tfx\tfz\tmxZ\n" << endl;
	//writing header
	os << "x\t";
	for (p=0;p<6;p++){
		os << "fx(from stress-" << p << ")\t";
	}
	os << "fx-total\t";
	for (p=0;p<6;p++){
		os << "fz(from stress-" << p << ")\t";
	}
	os << "fz-total\t";
	for (p=0;p<6;p++){
		os << "mxZ(from stress-" << p << ")\t";
	}
	os << "mxZ-total\t";
	os << "ForceResultant-offset\t";
	os << "\n";
	
for (p=0;p<6;p++){
	
	// this prints only the values from the "left end"
	for (i=0;i<CrossSectionList.getNumElements();i++){
		csi=&CrossSectionList[i];
		getGlobalForces(globalForces, solution, &(csi->elementList), p);
		cout << "."; // just to see progress...

		for (k=0;k<3;k++){
			if(k==1) continue;
			if(k==2) continue; // skipping the right end
			fx=fy=fz=0.0;
			mxZ=mxY=myX=myZ=mzX=mzY=0.0;
			numNodes=csi->nList[k].getNumElements();
			for(j=0;j<numNodes;j++){
				n=&csi->nList[k][j];
				firstDof=n->getFirstDof();
				fx += globalForces[firstDof];
				fy += globalForces[firstDof+1];
				fz += globalForces[firstDof+2];

				mxZ += globalForces[firstDof]   * (n->z - csi->centroid[k].z);
				mxY += globalForces[firstDof]   * (n->y - csi->centroid[k].y);

				myX += globalForces[firstDof+1] * (n->x - csi->centroid[k].x);
				myZ += globalForces[firstDof+1] * (n->z - csi->centroid[k].z);

				mzX += globalForces[firstDof+2] * (n->x - csi->centroid[k].x);
				mzY += globalForces[firstDof+2] * (n->y - csi->centroid[k].y);
			}
			//os << csi->centroid[k](direction) << "\t";
			/*
			os << fx << "\t" << fy << "\t" << fz << "\t";
			os << mxZ << "\t" << mxY << "\t";
			os << myX << "\t" << myZ << "\t";
			os << mzX << "\t" << mzY << "\t" << myZ+mzY << endl;
			*/
			csi->data[6*p + 0]= -fx;
			csi->data[6*p + 1]= -fz;
			csi->data[6*p + 2]= -mxZ;
			//os << -fx << "\t" << -fz << "\t";
			//os << -mxZ << endl;
		}
		//BETA_OUT << "=============" << endl;
	}
	
	// this prints only the values from the "right end"
	//os << endl;
	for (i=0;i<CrossSectionList.getNumElements();i++){
		csi=&CrossSectionList[i];
		getGlobalForces(globalForces, solution, &(csi->elementList), p);
		for (k=0;k<3;k++){
			if(k==1) continue;
			if(k==0) continue; // skipping the left end
			fx=fy=fz=0.0;
			mxZ=mxY=myX=myZ=mzX=mzY=0.0;
			numNodes=csi->nList[k].getNumElements();
			for(j=0;j<numNodes;j++){
				n=&csi->nList[k][j];
				firstDof=n->getFirstDof();
				fx += globalForces[firstDof];
				fy += globalForces[firstDof+1];
				fz += globalForces[firstDof+2];

				mxZ += globalForces[firstDof]   * (n->z - csi->centroid[k].z);
				mxY += globalForces[firstDof]   * (n->y - csi->centroid[k].y);

				myX += globalForces[firstDof+1] * (n->x - csi->centroid[k].x);
				myZ += globalForces[firstDof+1] * (n->z - csi->centroid[k].z);

				mzX += globalForces[firstDof+2] * (n->x - csi->centroid[k].x);
				mzY += globalForces[firstDof+2] * (n->y - csi->centroid[k].y);
			}
			//os << csi->centroid[k](direction) << "\t";
			/*
			os << fx << "\t" << fy << "\t" << fz << "\t";
			os << mxZ << "\t" << mxY << "\t";
			os << myX << "\t" << myZ << "\t";
			os << mzX << "\t" << mzY << "\t" << myZ+mzY << endl;
			*/
			csi->data[6*p + 3]= fx;
			csi->data[6*p + 4]= fz;
			csi->data[6*p + 5]= mxZ;

			//os << fx << "\t" << fz << "\t";
			//os << mxZ << endl;
		}
		//BETA_OUT << "=============" << endl;
	}
}

	double totalfx,totalfz,totalmxz;
	int numCSI=CrossSectionList.getNumElements();
	for (i=0;i<numCSI;i++){
		csi=&CrossSectionList[i];
		os << csi->centroid[0](direction) << "\t";
		totalfx=0;
		for(j=0;j<6;j++){//for fx
			if(i==0){
				os << csi->data[6*j + 0] << "\t";
				totalfx += csi->data[6*j + 0];
			}else{
				csi2=&CrossSectionList[i-1];
				os << (csi->data[6*j + 0] + csi2->data[6*j + 3])/2 << "\t";
				totalfx += (csi->data[6*j + 0] + csi2->data[6*j + 3])/2;
			}				
		}
		os << totalfx << "\t";
		totalfz=0;
		for(j=0;j<6;j++){//for fz
			if(i==0){
				os << csi->data[6*j + 1] << "\t";
				totalfz += csi->data[6*j + 1];
			}else{
				csi2=&CrossSectionList[i-1];
				os << (csi->data[6*j + 1] + csi2->data[6*j + 4])/2 << "\t";
				totalfz += (csi->data[6*j + 1] + csi2->data[6*j + 4])/2;
			}				
		}
		os << totalfz << "\t";
		totalmxz=0;
		for(j=0;j<6;j++){//for mxZ
			if(i==0){
				os << csi->data[6*j + 2] << "\t";
				totalmxz += csi->data[6*j + 2];
			}else{
				csi2=&CrossSectionList[i-1];
				os << (csi->data[6*j + 2] + csi2->data[6*j + 5])/2 << "\t";
				totalmxz += (csi->data[6*j + 2] + csi2->data[6*j + 5])/2;
			}				
		}
		os << totalmxz << "\t";
		os << totalmxz/totalfx << "\t";
		os << endl;

		//for last cross section layer:
		if(i==numCSI-1){
			os << csi->centroid[2](direction) << "\t";
			totalfx=0;
			for(j=0;j<6;j++){//for fx
				os << csi->data[6*j + 3] << "\t";
				totalfx += csi->data[6*j + 3];
			}
			os << totalfx << "\t";
			totalfz=0;
			for(j=0;j<6;j++){//for fz
				os << csi->data[6*j + 4] << "\t";
				totalfz += csi->data[6*j + 4];
			}
			os << totalfz << "\t";
			totalmxz=0;
			for(j=0;j<6;j++){//for mxZ
				os << csi->data[6*j + 5] << "\t";
				totalmxz += csi->data[6*j + 5];
			}
			os << totalmxz << "\t";
			os << totalmxz/totalfx << "\t";
			os << endl;
		}
	}

	cout << endl;

//	os.close();
	filemanager->CloseOutputStream(ost);
}
//=====================================================================
//this used to be the older version - without split up based on stress component
void TowTraversalUtility::getGlobalForcesAtCrossSectionsNormal()
{
	cout << "started getGlobalForcesAtCrossSectionsNormal()" << endl;

//	allocateForAnalysis(false);//not allocating equations //JV030409: did not allocate eqns in hs4rae, but in beta currently we have to
	allocateForAnalysis(true);
	setSolutionVectorFromMeshDisplacements(mesh);
	//setMaterial();
	

	Node *n;
	int i,j,k;
	int numNodes;

	CrossSectionInfo *csi;
	double *globalForces=internalForce;


	double fx,fy,fz;
	double mxZ, mxY, myX, myZ, mzX, mzY;
	int firstDof;
	int direction=0; // hardcoded for now

//	ofstream os("ForceMoments_Normal.xls");
	stringstream number;
	number << towpath.MaterialGroup;
	string filename="ForceMoments_Normal_material_" + number.str() + ".xls";
	ofstream *ost = filemanager->OpenOutputStream(filename.c_str());
	ofstream &os=*ost;


	os << "x\tfx\tfy\tfz\tmxZ\tmxY\tmyX(=0)\tmyZ(twist)\tmzX(=0)\tmzY(twist)\tMtwist=MyZ+MzY\n" << endl;
	
	// this prints only the values from the "left end"
	os << "Left face of layers" << endl;
	for (i=0;i<CrossSectionList.getNumElements();i++){
		csi=&CrossSectionList[i];
		getGlobalForces(globalForces, solution, &(csi->elementList), -1);
		cout << "."; // just to see progress...

		for (k=0;k<3;k++){
			if(k==1) continue;
			if(k==2) continue; // skipping the right end
			fx=fy=fz=0.0;
			mxZ=mxY=myX=myZ=mzX=mzY=0.0;
			numNodes=csi->nList[k].getNumElements();
			for(j=0;j<numNodes;j++){
				n=&csi->nList[k][j];
				firstDof=n->getFirstDof();
				fx += globalForces[firstDof];
				fy += globalForces[firstDof+1];
				fz += globalForces[firstDof+2];

				mxZ += globalForces[firstDof]   * (n->z - csi->centroid[k].z);
				mxY += globalForces[firstDof]   * (n->y - csi->centroid[k].y);

				myX += globalForces[firstDof+1] * (n->x - csi->centroid[k].x);
				myZ += globalForces[firstDof+1] * (n->z - csi->centroid[k].z);

				mzX += globalForces[firstDof+2] * (n->x - csi->centroid[k].x);
				mzY += globalForces[firstDof+2] * (n->y - csi->centroid[k].y);
			}
			os << csi->centroid[k](direction) << "\t";
			os << fx << "\t" << fy << "\t" << fz << "\t";
			os << mxZ << "\t" << mxY << "\t";
			os << myX << "\t" << myZ << "\t";
			os << mzX << "\t" << mzY << "\t" << myZ+mzY << endl;
		}
		//BETA_OUT << "=============" << endl;
	}
	os << endl;
	// this prints only the values from the "right end"
	os << "Right face of layers" << endl;
	for (i=0;i<CrossSectionList.getNumElements();i++){
		csi=&CrossSectionList[i];
		getGlobalForces(globalForces, solution, &(csi->elementList), -1);
		for (k=0;k<3;k++){
			if(k==1) continue;
			if(k==0) continue; // skipping the left end
			fx=fy=fz=0.0;
			mxZ=mxY=myX=myZ=mzX=mzY=0.0;
			numNodes=csi->nList[k].getNumElements();
			for(j=0;j<numNodes;j++){
				n=&csi->nList[k][j];
				firstDof=n->getFirstDof();
				fx += globalForces[firstDof];
				fy += globalForces[firstDof+1];
				fz += globalForces[firstDof+2];

				mxZ += globalForces[firstDof]   * (n->z - csi->centroid[k].z);
				mxY += globalForces[firstDof]   * (n->y - csi->centroid[k].y);

				myX += globalForces[firstDof+1] * (n->x - csi->centroid[k].x);
				myZ += globalForces[firstDof+1] * (n->z - csi->centroid[k].z);

				mzX += globalForces[firstDof+2] * (n->x - csi->centroid[k].x);
				mzY += globalForces[firstDof+2] * (n->y - csi->centroid[k].y);
			}
			os << csi->centroid[k](direction) << "\t";
			os << fx << "\t" << fy << "\t" << fz << "\t";
			os << mxZ << "\t" << mxY << "\t";
			os << myX << "\t" << myZ << "\t";
			os << mzX << "\t" << mzY << "\t" << myZ+mzY << endl;
		}
		//BETA_OUT << "=============" << endl;
	}

	cout << endl;
	filemanager->CloseOutputStream(ost);
}

//=============================================================
/*
void TowTraversalUtility::getGlobalForcesAtCrossSections()
{
	ModelHandler_FE *mh;
	mh=(ModelHandler_FE *)modelHandler;
	mh->allocateForStaticAnalysis(mh->elementWorkspace);
	mh->setSolutionVectorFromMeshDisplacements(mesh);
	mh->setMaterial();
	

	BasicElement *e;
	Node *n,*oldn;
	int i,j,k;
	int numNodes;

	CrossSectionInfo *csi;
	double *globalForces=mh->internalForce;


	double fx,fy,fz;
	double mx,my,mz;
	double prevValue=9999999999.999999999999;
	int firstDof;
	int direction=0; // hardcoded for now

	BETA_OUT << "\nForces and moments along the tow path\n" << endl;
	BETA_OUT << "\n\tx\tfx\tfy\tfz\tmx\tmy\tmz\n" << endl;

	for (i=0;i<CrossSectionList.getNumElements();i++){
		fx=fy=fz=0.0;
		mx=my=mz=0.0;
		csi=CrossSectionList[i];
		mh->getGlobalForces(globalForces,&(csi->elementList));
		numNodes=csi->nodeList.getNumElements();
		for(j=0;j<numNodes;j++){
			n=csi->nodeList[j];
			if(j==0) prevValue = (*n)(direction);
			if(prevValue != (*n)(direction) && (j!=0) ){
				BETA_OUT << (*oldn)(direction) << "\t" << fx << "\t" << fy << "\t" << fz << "\t";
				BETA_OUT << mx << "\t" << my << "\t" << mz << endl;
				fx=fy=fz=0.0;
				mx=my=mz=0.0;
				prevValue = (*n)(direction);
			}
			firstDof=n->getFirstDof();
			fx += globalForces[firstDof];
			fy += globalForces[firstDof+1];
			fz += globalForces[firstDof+2];

			mx += globalForces[firstDof] * (n->x - csi->centroid.x);
			my += globalForces[firstDof+1] * (n->y - csi->centroid.y);
			mz += globalForces[firstDof+2] * (n->z - csi->centroid.z);


			if( j==numNodes-1 ){
				BETA_OUT << (*n)(direction) << "\t" << fx << "\t" << fy << "\t" << fz <<"\t";
				BETA_OUT << mx << "\t" << my << "\t" << mz << endl;
			}
			oldn=n;
		}
		BETA_OUT << "=============" << endl;
	}
}
*/
Node calculateCentroid(NodeGroup * nodeList)
{
	//this include all the nodes in the layer of elements at the Cross section.
	//so there will be 3 centroids.

	Node centroid;
	int i,numNodes;
	double x,y,z;
	x=y=z=0.0;
	numNodes=nodeList->getNumElements();
	for(i=0;i<numNodes;i++){
		x += (*nodeList)[i].x;
		y += (*nodeList)[i].y;
		z += (*nodeList)[i].z;
	}
	x /= numNodes;
	y /= numNodes;
	z /= numNodes;

	centroid.x=x;
	centroid.y=y;
	centroid.z=z;

	return centroid;

}

#undef TOWPATH_TOLERANCE

//===============================================================
void TowTraversalUtility::getGlobalForces(double *internalForce, double *theSolution, ElementGroup *elist, int NonZeroStressComponent)
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

		((ElasticityElement3D*)e)->ForceSingleNonZeroStressComponent=NonZeroStressComponent;
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
