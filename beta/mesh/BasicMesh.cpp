#include "stdafx.h"

#include "BasicMesh.hpp"
#include "utility/utility.h"
#include "MeshUtility.h"
#include "NodeSet.hpp"

int getNodeIndex(NodeGroup *nlist, double x, double y, double z);
bool readDisplacements(istream &fp, int num_nodes, int dims, Node **displacements);
bool readMaterialGroupsFromFile(const char *filename, BasicMesh* mesh);
void GetMinMax(NodeGroup *nlist, Node *Min, Node *Max, bool initializeMinMax);
void GetMinMax(ElementGroup *elist, Node *Min, Node *Max, bool initializeMinMax);
BoundaryNodeList * getMeshBoundaryNodes(ElementGroup *ele, int numElements, int numNodes);
//-------------------------------------------------------------------
CreateErrorHandler(BasicMesh);
//--------------------------------------------------------------------
BasicMesh::BasicMesh()
{
initializeBasicMesh();
//initializeStatusFlags();

}


void BasicMesh::initializeBasicMesh() 
{
	filemanager			=0;
	numElements        = 0;
	numNodes = numDims = 0;
	element.SetName("Mesh_ElementList");
	node.SetName("Mesh_NodeList");

	ElementGroupsList.SetName("Mesh_ElementGroupsList");

	haveSetElementType=haveReadMesh=false;

	nodeOffset=elementOffset=dofOffset=0;
	xOffset=yOffset=zOffset=0;
	displacements		=0;
    NodeSearchTolerance = 1.0e-8;
}


BasicMesh::~BasicMesh(void)
{
	clearMesh(); //deletes NodeGroupLists and ElementGroupLists

	// initializing data
	numElements        = 0;
	numNodes = numDims = 0;
	haveSetElementType=haveReadMesh=false;

	if(displacements)		free(displacements);	displacements=0;
}
//=========================================================
void BasicMesh::clearMesh()
{
	int i;
	for(i=0;i<ElementGroupsList.getNumElements();i++){
		ElementGroupsList[i].clear();//so that it doesn't try to delete the elements, 
	}
}
//=========================================================
void BasicMesh::ReadMesh(char *filename)
{
	ifstream is;
	is.open(filename);
	ReadMesh(&is);
	is.close();
}
//=========================================================
BasicElement* BasicMesh::createElement_adhoc(int numDims, int nnpe)
{
	BasicElement *e=0;
	e = new BasicElement;
	return e;
}
//=========================================================
void BasicMesh::ReadMesh(istream *inputStream)
{
	CoordinateReadOption readOption;	
	int i, number;

	// initializing data 
	numNodes = numDims = 0;


// read in total numbers
	(*inputStream)>>numNodes;
	(*inputStream)>>number;
	(*inputStream)>>numDims; 

	if(numDims==1){readOption = X;}
	if(numDims==2){readOption = XY;}
	if(numDims==3){readOption = XYZ;}

//Check for special commands

	string wholeLine;
	string token;
	wholeLine=peekAtNextLine(inputStream);
	if( wholeLine== "processCommands"){
		BETA_OUT<<"processCommands\n";
		(*inputStream)>>token; // advance file position
		while(2==2){
			(*inputStream)>>token;  BETA_OUT<<token<<endl; 
			if(token =="endCommands")  break;
			if(token== "setCoordinateReadOption") 
			(*inputStream)>>token;
			BETA_OUT<<"readOption ="<<token<<endl;
			if(token == "xyz"){readOption = XYZ;}
			if(token == "xy" ){readOption = XY;}
			if(token == "yz" ){readOption = YZ;}
			if(token == "x"  ){readOption = X;}
		};//while
	}//if processCommands
//=======================================

	Node			*n;
	BasicElement	*e=0;
	bool useAdhocFactory=false;

	if(number != numElements){
		if (numElements==0) {
			cout << "elements not created prior to ReadMesh(). Creating elements now..." << endl;
			numElements=number;
			element.reserve(numElements);
			useAdhocFactory=true;
		}else{
			cout << "The number of elements in the mesh(";
			cout << number << ") are not the same as the number of elements(";
			cout << numElements << ") created." << endl;
			exit(1);
		}
	}
	  	
	BETA_OUT<<"Number of Nodes    : "<<numNodes     <<endl;
	BETA_OUT<<"Number of Elements : "<<numElements  <<endl;
	BETA_OUT<<"Number of numDims  : "<<numDims      <<endl;

	if(abs(numDims) < 1 || abs(numDims) > 3) EXIT_BETA("numDims is invalid")

	// allocate memory for nodes
	node.reserve(numNodes);

// Read Nodes
	for(i=0; i< numNodes; i++)
	{
		(*inputStream)>>number;

		if ((*inputStream).eof()) {BETA_OUT << "Bad BasicMesh File !!! Fix !!!" << endl; exit(1);} //JV042204
		n = new Node();
		n->read(*inputStream, readOption);
		n->nodeNum=number;
		n->newNodeNum=number;
		node.add(*n);
	}
	if(node[0].nodeNum==1) nodeOffset=1; else nodeOffset=0;
	//nodeOffset is also set in void BasicElement::read(istream &istrm,Node *nodeList)

// Read Elements
    for(i=0; i< numElements; i++)
	{
		int nnpe;
		(*inputStream) >> number >> nnpe;
		if ((*inputStream).eof()) {EXIT_BETA("Bad BasicMesh File !!! Fix !!!");} //JV042204
		if(useAdhocFactory){
			e=createElement_adhoc(numDims, nnpe);
			element.putAt(number,e);
		}
		else
			e = &element[number];
		e->setFileManager(filemanager); // set filemanager
		e->setElementNumber(number);
		e->setNewElementNumber(number);
		e->read(*inputStream, &node, nnpe);
		e->setElementType(numDims);  // 042902
		e->setMidNodes(numDims);
		e->initialize();
	}
	if(element[0].elementNumber==1) elementOffset=1; else elementOffset=0;
	//this offset is of no use when it comes to plotting meshes.
	//in the current version, is there is an offset, it does not plot the mesh properly
	//this problem started when the data structures were modified from the original version.

    // get Min, Max of the mesh - this goes through ALL the nodes in the mesh
	GetMinMax(&node, &min, &max, true);
	// get Min, Max of the mesh - this ONLY goes through the nodes that are connected to elements in the mesh. (ignores dummy nodes)
//	GetMinMax(&element, &min, &max, true);
}
//=========================================================
bool BasicMesh::ReadMeshBinary(string filename)
{   
	WhoMethod("ReadMeshBinary(char *)");
	int i;
	const int SIZEOF_INT=sizeof(int);
	FILE *meshFile;
 
	if((meshFile = fopen(filename.c_str(),"rb"))==NULL)
		EXIT_BETA("Cannot open binary file:" << filename);

	int number;

// initializing data 
	numNodes = numDims = 0;

// read in total numbers
	fread(&numNodes	,SIZEOF_INT,1,meshFile);
	fread(&number	,SIZEOF_INT,1,meshFile);
	fread(&numDims	,SIZEOF_INT,1,meshFile);

	Node			*n;
	BasicElement	*e=0;
	bool useAdhocFactory=false;

	if(number != numElements){
		if (numElements==0) {
			cout << "elements not created prior to ReadMesh(). Creating elements now..." << endl;
			numElements=number;
			element.reserve(numElements);
			useAdhocFactory=true;
		}else{
			cout << "The number of elements in the mesh(";
			cout << number << ") are not the same as the number of elements(";
			cout << numElements << ") created." << endl;
			exit(1);
		}
	}
	  	
	BETA_OUT<<"Number of Nodes    : "<<numNodes     <<endl;
	BETA_OUT<<"Number of Elements : "<<numElements  <<endl;
	BETA_OUT<<"Number of numDims  : "<<numDims      <<endl;

	if(abs(numDims) < 1 || abs(numDims) > 3) EXIT_BETA("numDims is invalid");

	node.reserve(numNodes);
	for(i=0;i<numNodes;i++){
		n = new Node();
		n->readBinary(meshFile);
		n->newNodeNum=n->nodeNum;
		node.add(*n);
	}
	if(node[0].nodeNum==1) nodeOffset=1; else nodeOffset=0;

	for(i=0;i<numElements;i++){
		int nnpe;
		fread(&number	,SIZEOF_INT,1,meshFile);
		fread(&nnpe		,SIZEOF_INT,1,meshFile);

		if(useAdhocFactory){
			e=createElement_adhoc(numDims, nnpe);
			element.putAt(number,e);
		}
		else
			e = &element[number];
		e->setFileManager(filemanager); // set filemanager
		e->setElementNumber(number);
		e->setNewElementNumber(number);
		e->readBinary(meshFile, &node, nnpe);
		e->setElementType(numDims);  // 042902
		e->setMidNodes(numDims);
		e->initialize();
	}

	if(element[0].elementNumber==1) elementOffset=1; else elementOffset=0;
	GetMinMax(&node, &min, &max, true);
	return true;
}
//=========================================================
void BasicMesh::PrintMesh(string modelname, bool activeOnly)
{
	ofstream os;
	string filename(modelname + ".plt");
	os.open(filename.c_str());
	PrintMesh(&os, activeOnly);
	os.close();
}
//=========================================================
void BasicMesh::PrintMesh(ostream * ostr, bool activeOnly)
{
	int i,j=0, k;
	BasicElement	*e=0;
	Node		*n=0;

	if (activeOnly) {
		*ostr <<getNumActiveNodes()<<" " <<getNumActiveElements()<<" "<<numDims<<endl;
		// nodes
		for(i=0; i< numNodes; i++){
			n=&node[i];
			if (n->isActive) {
				*ostr << n->newNodeNum<<"\t"<< DOUBLE_FORMAT << n->x<<"\t"<<n->y;
				if (numDims==3 || numDims==-2 /*xtang 081802*/) {*ostr<<"\t"<<n->z<<endl;} else {*ostr<<endl;}
/*
				if (numDims==3) {*ostr<<"\t"<<node[i].z<<"\t";} else {*ostr<<"\t";}
				if (node[i].numElementsAtNode > 0) {
					*ostr <<node[i].numElementsAtNode<<"\t";
					for (k=0; k<node[i].numElementsAtNode;k++) {
						*ostr << node[i].eList[k]->elementNumber<<"\t";
					}
					*ostr << endl;
				}
*/
			}
		}
		// elements
		j=0;
		for(i=0; i< numElements; i++) {
			e=&element[i];
			if (e->activeElementFlag) {
				*ostr <<j<<"\t"<<e->numNodesPerElement<<"\t";
				for(k=0; k< e->numNodesPerElement; k++) 
 				        *ostr << e->node[k].newNodeNum<<" ";
				*ostr<<endl;
				j++;
			}
		}
	}
	else {
		*ostr <<numNodes<<" " <<numElements<<" "<<numDims<<endl;
		// nodes
		for(i=0; i< numNodes; i++) {
			n=&node[i];
			*ostr << n->nodeNum<<"\t"<<n->x<<"\t"<<n->y;
			if (numDims==3 ||numDims==-2 /*xtang 081802*/) {*ostr<<"\t"<<n->z<<endl;} else {*ostr<<endl;}
/*
				if (numDims==3) {*ostr<<"\t"<<node[i].z<<"\t";} else {*ostr<<"\t";}
				if (node[i].numElementsAtNode > 0) {
					*ostr <<node[i].numElementsAtNode<<"\t";
					for (k=0; k<node[i].numElementsAtNode;k++) {
						*ostr << node[i].eList[k]->elementNumber<<"\t";
					}
					*ostr << endl;
				}
*/
		}
		// elements
		for(i=0; i< numElements; i++) {
			e=&element[i];
			*ostr <<e->elementNumber<<"\t"<<e->numNodesPerElement<<"\t";
			for(k=0; k< e->numNodesPerElement; k++) 
 				  *ostr << e->node[k].nodeNum<<" ";
			*ostr<<endl;
		}
	}
}

int BasicMesh::getNumActiveNodes()
{
	int i,j=0;
    for(i=0; i< numNodes; i++) if (node[i].isActive) j++;
	return j;
}
//=========================================================
int BasicMesh::getNumActiveElements()
{
	int i,j=0;
    for(i=0; i< numElements; i++) if (element[i].activeElementFlag) j++;
	return j;
}
//=========================================================
NodeSet& BasicMesh::getNodeSet(const string &NodeSetName)
{
    //return reference to the nodeset with the input name
    //note that the behavior of the STL map is to create a new
    //object if there is no existing object with that name
    // First, set the tolerance to the current tolerance for the mesh...
    NodeSets[NodeSetName].setNodeSearchTolerance(NodeSearchTolerance);
    return NodeSets[NodeSetName];
}
//=========================================================
const NodeSet& BasicMesh::getExistingNodeSet(const string &NodeSetName) const
{
    //Returns a const reference to a NodeSet.  Gives error if the nodeset
    //doesn't already exist.
    map<const string, NodeSet>::const_iterator NodeSetPairIterator = NodeSets.find(NodeSetName);
    if(NodeSetPairIterator == NodeSets.end()){
        cout << "Cannot find nodeset by name \"" << NodeSetName << "\", exiting" << endl;
        exit(1);
    }
    return NodeSetPairIterator->second;
}
//=========================================================
const NodeSet* BasicMesh::getExistingNodeSetPtr(const string &NodeSetName) const
{
    //Returns a const pointer to a NodeSet.  Gives error if the nodeset
    //doesn't already exist.
    map<const string, NodeSet>::const_iterator NodeSetPairIterator = NodeSets.find(NodeSetName);
    if(NodeSetPairIterator == NodeSets.end()){
        cout << "Cannot find nodeset by name \"" << NodeSetName << "\", exiting" << endl;
        exit(1);
    }
    return &NodeSetPairIterator->second;
}
//=========================================================
NodeSet BasicMesh::getNodesOnPlane(int coordnum, double coordval)
{
    NodeSet result(NodeSearchTolerance);
    result.addNodesOnPlane(this,coordnum,coordval);
	return result;
}
//=========================================================
NodeSet BasicMesh::getNodesUsingCoordinates(double x, double y, double z)
{
    NodeSet result(NodeSearchTolerance);
    double coord[] = {x,y,z};
    result.addNodesOnPoint(this,coord);
	return result;
}
//==========================================================================
NodeSet BasicMesh::getNodesUsingCoordinates(Node *n)
{
    NodeSet result(NodeSearchTolerance);
    double coord[] = {n->x,n->y,n->z};
    result.addNodesOnPoint(this,coord);
	return result;
}
//==========================================================================
int BasicMesh::getNodeIndex(double x, double y, double z)
{
	return ::getNodeIndex(&node, x, y, z);
}
//==========================================================================
bool BasicMesh::readDisplacements(const char *filename)
{
	ifstream infile(filename);
	if(!infile) return false;

	//if(!(::readDisplacements(infile,numNodes,numDims, &displacements))) 
	//{infile.close();return false;} //JV052004 - have to stop if false !!!

//============================================start of copy
		int i,junk;
	double junkdisp;

	if(displacements) free(displacements);
	displacements = (Node *)malloc(numNodes * sizeof(Node));

	if (!infile) {
		cerr << "Error: <readDisplacement> - Bad istream.\n";
		return false;
	}


char * tokenList[20];
int numberOfColumns;
	getLineAndTokenize(&infile, "Dummy",tokenList,numberOfColumns); // flush first line
	getLineAndTokenize(&infile, "Dummy",tokenList,numberOfColumns);
if(numberOfColumns==2)
{
// (*displacements)[0].x = 0.f;
// (*displacements)[0].y = (float)atof(tokenList[1]);
// (*displacements)[0].z = 0.f;
//JV020905 temporarily changed to plot contours on 2d plots properly ( with 1 dofpn)
 displacements[0].x = (float)atof(tokenList[1]);
 displacements[0].y = 0.f;
 displacements[0].z = 0.f;
}

if(numberOfColumns==3)
{
 displacements[0].x = (float)atof(tokenList[1]);
 displacements[0].y = (float)atof(tokenList[2]);
 displacements[0].z = 0.f;
}
//if(numberOfColumns==4)
if(numberOfColumns >= 4) //JV0333005 temporary fix: what if there are more than 3 dofs/displacements per node??
{
 displacements[0].x = (float)atof(tokenList[1]);
 displacements[0].y = (float)atof(tokenList[2]);
 displacements[0].z = (float)atof(tokenList[3]);
}

//Start at 1 instead of zero becasue first line of data alread read.
	switch (numberOfColumns) {
		case 2: 
			for (i = 1; i < numNodes; i++) {
//				infile >> junk >> displacements[i].y ; // see similar comment above
				infile >> junk >> displacements[i].x ;
				if (!infile) {
					cerr << "Error: <readDisplacement> - Bad istream.\n";
					return false;
				}
				displacements[i].y = 0.0;
				displacements[i].z = 0.0;
			}
			break;
		case 3: 
			for (i = 1; i < numNodes; i++) {
				//infile >> junk >> displacements[i].y >> displacements[i].z; //JV072103
				infile >> junk >> displacements[i].x >> displacements[i].y;  //JV072103
				if (!infile) {
					cerr << "Error: <readDisplacement> - Bad istream.\n";
					return false;
				}
				//displacements[i].x = 0.0;  //JV072103
				displacements[i].z = 0.0;  //JV072103
			}
			break;
// If there are more four columns, then there must be extra dof ... ignore them
		case 4:  //this is temporary fix.   
		case 5:  //this is temporary fix. need a more robust way to plot contours if there are more than 3 'displacements'
		case 6: 
			for (i = 1; i < numNodes; i++) {
				if(numberOfColumns==4)
					infile >> junk >> displacements[i].x >> displacements[i].y >> displacements[i].z;
				else if(numberOfColumns==5)
					infile >> junk >> displacements[i].x >> displacements[i].y >> displacements[i].z >> junkdisp;

				if (!infile) {cerr << "Error: <readDisplacement> - Bad istream.\n";
					return false;}
			}//i
			break;
		default: 
			cerr << "Error: <readDisplacements> - numberOfColumns in readDisplacements = " << numberOfColumns << endl;
			return false;
	}
//============================================end of copy


	infile.close();
	return true;	
}
//==========================================================================
bool BasicMesh::readMaterialGroupsFromFile(const char *filename)
{
	return ::readMaterialGroupsFromFile(filename,this);
}
//==========================================================================
void BasicMesh::assignNewNodeNum()
{
	int i,j=0;
    for(i=0; i< numNodes; i++) {
		if (node[i].isActive) {
			node[i].newNodeNum=j;
			j++;
		}
		else {
			node[i].newNodeNum=-1;
		}
	}
}
//==========================================================================
void BasicMesh::assignNewElementNum()
{
	int i,j=0;
    for(i=0; i< numElements; i++) {
		if (element[i].activeElementFlag) {
			element[i].newElementNumber=j;
			j++;
		}
		else {
			element[i].newElementNumber=-1;
		}
	}
}
//==========================================================================
void BasicMesh::SetActiveState(bool flag)
{
	int i;
	for (i=0;i<numElements;i++)
		element[i].activeElementFlag = flag;
	for (i=0;i<numNodes;i++)
		node[i].isActive=flag;
}
//==================================================================
void BasicMesh::InitializeOutputSummary()
{
	for(int i=0;i<numElements;i++){
		element[i].InitializeOutputSummary();
	}
}
//==================================================================
void BasicMesh::OutputSummary()
{
	for(int i=0;i<numElements;i++){
		element[i].OutputSummary();
	}
}
//==================================================================
void BasicMesh::FinalizeSummary(string modelname)
{
	for(int i=0;i<numElements;i++){
		element[i].FinalizeSummary();
	}

	string filename;
	ofstream ofile;
	//this function has to know about the different types of files to be created.
	//elemat file :
	for(int i=0;i<numElements;i++) element[i].setConsolidate(false);
	filename = modelname + ".elemat";
	ofile.open(filename.c_str());
	ofile << "readmaterialgroup()" << endl;
	for(int i=0;i<numElements;i++){
		element[i].ConsolidateSummary("ELEMAT", ofile);
	}
	ofile << "-1 0 0" << endl;
	ofile.close();

	//mangles file :
	for(int i=0;i<numElements;i++) element[i].setConsolidate(false);
	filename = modelname + ".mangles";
	ofile.open(filename.c_str());
	ofile << "setElementMaterialAngles" << endl;
	for(int i=0;i<numElements;i++){
		element[i].ConsolidateSummary("MANGLES", ofile);
	}
	ofile.close();


}
//==================================================================
bool BasicMesh::readContourData(istream &infile, string filename)
{
	if(!infile) return false;

   contourDataSet.setPointers(this);
   contourDataSet.setFileManager(filemanager);
   contourDataSet.pickColumn = 1;
   contourDataSet.setFilename(filename);
   if (!contourDataSet.readElementNodalValues(infile))
	   return false;

//   contourDataSet.getMinMaxValues(minVal,maxVal);
   BETA_OUT<<"Set contourFilename : "<<filename<<endl;

//	statusFlags.haveData = true;
//	minContourAuto = minVal;
//	maxContourAuto = maxVal;
	//BETA_OUT<<"minVal="<<minVal<<endl;
    //BETA_OUT<<"maxVal="<<maxVal<<endl;
	return true;	
}
//==========================================================================
bool BasicMesh::readContourData(string filename)
{
	/*
	ifstream is;
	is.open(filename);
	if(!is) return false;
	bool result=readContourData(is, filename);
	is.close();
	return result;
	*/

	ifstream *is;
	is=filemanager->OpenInputStream(filename);
	if(is==0) return false;
	bool result=readContourData(*is, filename);
	filemanager->CloseInputStream(is);
	return result;
}
//==========================================================================
bool BasicMesh::allocateNodalValuesElements()
{
	if(element[0].nodalValues) return true;

	for(int i=0; i<numElements; i++){
		element[i].nodalValues=new double[element[i].numNodesPerElement];
	}
	return true;
}
//==========================================================================
void BasicMesh::setupNodeDoc(bool isActiveOnly)
{
	int i,j;
	Node *aN; 

	for(i=0; i< numNodes; i++)
	{
		if (node[i].numElementsAtNode > 0) {
			delete [] node[i].eList;
		    node[i].numElementsAtNode =0;
			node[i].eList=0;
		}
	}

	Node *n;
    for(i=0; i< numElements; i++)
	{
		n=0;
		 if ((isActiveOnly &&  element[i].activeElementFlag) || !isActiveOnly) 
		 for (j=0; j<element[i].numNodesPerElement; j++) {
			 if(n==&element[i].node[j])
				 continue;
			 else
				 element[i].node[j].numElementsAtNode ++;
			 n=&element[i].node[j];
		 }
	}

	for(i=0; i< numNodes; i++)
	{
		if (node[i].numElementsAtNode>0) {
			node[i].eList = new BasicElement * [node[i].numElementsAtNode];
			node[i].numElementsAtNode=0;
		}
	}
    for(i=0; i< numElements; i++)
	{
		n=0;
		 if ((isActiveOnly &&  element[i].activeElementFlag) || !isActiveOnly) 
		 for (j=0; j<element[i].numNodesPerElement; j++) {
			 if(n==&element[i].node[j])
				 continue;
			 else{
				 aN=&element[i].node[j];
				 aN->eList[aN->numElementsAtNode]=&element[i];
				 aN->numElementsAtNode ++;
			 }
			 n=&element[i].node[j];
		 }
	}
}
//===================================================================
void BasicMesh::calculateElementCentroid()
{
	::calculateElementCentroid(&element);
}
//===================================================================
BoundaryNodeList * BasicMesh::findBoundaryNodes()
{
	setupNodeDoc(true);
	BoundaryNodeList *m1=getMeshBoundaryNodes(&element, numElements, numNodes);
	if (m1) {
		string aName;
		aName="mesh_Nodes_" + meshFileName;
		m1->printBoundaryNodeList(aName.c_str());
	}
	return m1;
}
