#include "stdafx.h"

#include "BasicMesh.hpp"
#include "utility/utility.h"
#include "MeshUtility.h"

//======================================================================
ElementGroup*  getElementsInsideVolume(ElementGroup &elist, Node &n1, Node &n2);
void ExtractSubMeshUsingMaterialGroups(BasicMesh &mesh, Array<int> &materialgroups, string  filenamePrefix);
void ExtractSubMeshUsingVolume(BasicMesh &mesh, Node &n1, Node &n2, string  filenamePrefix);
void GetMinMax(NodeGroup *nlist, Node *Min, Node *Max, bool initializeMinMax);
void GetMinMax(ElementGroup *elist, Node *Min, Node *Max, bool initializeMinMax);
//======================================================================
void ExtractSubMesh(BasicMesh &mesh, istream  &is)
{
	int i;
	
	string modelname;
	is >> modelname;

	string command;
	is >> command;
	ChangeToUpper(command)
	if(command=="READMATERIALGROUPS"){
		//obtain list of material numbers 
		Array<int> materialgroups;
		int numMaterialGroups;
		is >> numMaterialGroups;
		materialgroups.reserve(numMaterialGroups);
		int *newMaterialNumber;
		int MaterialNumber;
		for(i=0;i<numMaterialGroups; i++){
			is >> MaterialNumber;
			newMaterialNumber=new int(MaterialNumber);
			materialgroups.add(*newMaterialNumber);
		}
		ExtractSubMeshUsingMaterialGroups(mesh, materialgroups, modelname);
	}
	else if(command=="READVOLUMEBOUNDS"){
		Node n1,n2;
		is >> n1.x >> n1.y >> n1.z;
		is >> n2.x >> n2.y >> n2.z;
		ExtractSubMeshUsingVolume(mesh, n1, n2, modelname);

	}else{
		cerr << "Unrecognized command in ExtractSubMesh()!" << endl;
		exit(1);
	}
}
//======================================================================
//decide where to put this function!
int getIndex(int numberToFind, int * intArray, int length) 
{// Determine index of integer in intArray that matches specified integer
	int i=0;
	while (i<length) {
		if (numberToFind == intArray[i]) return i;
		i++;
	}
	return -1;
}
//======================================================================
int getIndex(Node *n, NodeGroup &nodeArray) 
{// Determine index of node in nodeArray that matches specified node
	for(int i=0; i<nodeArray.getNumElements(); i++){
		if ( n == &(nodeArray[i]) ) return i;
	}
	return -1;
}
//======================================================================
// Change to use seek and tellg to figure out number of columns
bool readDisplacements(istream &fp, int num_nodes, int dims,
		Node **displacements)
{// This is currently mainly (only?) used by the plotter
	int i,junk;
	double junkdisp;

	if(*displacements) free(*displacements);
	*displacements = (Node *)malloc(num_nodes * sizeof(Node));

	if (!fp) {
		cerr << "Error: <readDisplacement> - Bad istream.\n";
		return false;
	}

// Could use notation (*displacements) instead of displacements[0] ??


// Get first line of data to determine how many columns there are
// Alternative:
// Use the info from alpha input about ndof per node & let nodes read the displacement data
char * tokenList[20];
int numberOfColumns;
	getLineAndTokenize(&fp, "Dummy",tokenList,numberOfColumns); // flush first line
	getLineAndTokenize(&fp, "Dummy",tokenList,numberOfColumns);
if(numberOfColumns==2)
{
// (*displacements)[0].x = 0.f;
// (*displacements)[0].y = (float)atof(tokenList[1]);
// (*displacements)[0].z = 0.f;
//JV020905 temporarily changed to plot contours on 2d plots properly ( with 1 dofpn)
 (*displacements)[0].x = (float)atof(tokenList[1]);
 (*displacements)[0].y = 0.f;
 (*displacements)[0].z = 0.f;
}

if(numberOfColumns==3)
{
 //(*displacements)[0].x = 0.f;
 //(*displacements)[0].y = (float)atof(tokenList[1]);
 //(*displacements)[0].z = (float)atof(tokenList[2]);
//JV072103 - the above 3 lines are wrong !!! 
	//when plotting, the coords are deformed BEFORE rotation,(rotation makes the y-z plane as drawing plane) 
	// so the coords should be read normally ( i.e. in this order: x,y,z)
 (*displacements)[0].x = (float)atof(tokenList[1]);
 (*displacements)[0].y = (float)atof(tokenList[2]);
 (*displacements)[0].z = 0.f;
}
//if(numberOfColumns==4)
if(numberOfColumns >= 4) //JV0333005 temporary fix: what if there are more than 3 dofs/displacements per node??
{
 (*displacements)[0].x = (float)atof(tokenList[1]);
 (*displacements)[0].y = (float)atof(tokenList[2]);
 (*displacements)[0].z = (float)atof(tokenList[3]);
}

//Start at 1 instead of zero becasue first line of data alread read.
	switch (numberOfColumns) {
		case 2: 
			for (i = 1; i < num_nodes; i++) {
//				fp >> junk >> (*displacements)[i].y ; // see similar comment above
				fp >> junk >> (*displacements)[i].x ;
				if (!fp) {
					cerr << "Error: <readDisplacement> - Bad istream.\n";
					return false;
				}
				(*displacements)[i].y = 0.0;
				(*displacements)[i].z = 0.0;
			}
			break;
		case 3: 
			for (i = 1; i < num_nodes; i++) {
				//fp >> junk >> (*displacements)[i].y >> (*displacements)[i].z; //JV072103
				fp >> junk >> (*displacements)[i].x >> (*displacements)[i].y;  //JV072103
				if (!fp) {
					cerr << "Error: <readDisplacement> - Bad istream.\n";
					return false;
				}
				//(*displacements)[i].x = 0.0;  //JV072103
				(*displacements)[i].z = 0.0;  //JV072103
			}
			break;
// If there are more four columns, then there must be extra dof ... ignore them
		case 4:  //this is temporary fix.   
		case 5:  //this is temporary fix. need a more robust way to plot contours if there are more than 3 'displacements'
		case 6: 
			for (i = 1; i < num_nodes; i++) {
				if(numberOfColumns==4)
					fp >> junk >> (*displacements)[i].x >> (*displacements)[i].y >> (*displacements)[i].z;
				else if(numberOfColumns==5)
					//this if block should handle quasi-3d as well
					fp >> junk >> (*displacements)[i].x >> (*displacements)[i].y >> (*displacements)[i].z >> junkdisp;

				if (!fp) {
					cerr << "Error: <readDisplacement> - Bad istream.\n";
					return false;
				}
			}//i
			break;
		default: 
			cerr << "Error: <readDisplacements> - dims in readDisplacements, dims = " << dims << endl;
			return false;
	}
	return true;
}
//=========================================================================
void calculateElementCentroid(ElementGroup *elements)
{
	double           tempx,tempy,tempz;
  	int             j,nodesPerElem,i;
	BasicElement	*ElemJ;
	Node			*n;
	int nodeNum;
	int numElems=elements->getNumElements();

	for (j = 0; j < numElems; j++) {
		ElemJ = &(*elements)[j];
    	nodesPerElem = ElemJ->numNodesPerElement;
		/* we add the coord for each element */
    	tempx = 0.0;
    	tempy = 0.0;
    	tempz = 0.0;
    	for (i = 0; i < nodesPerElem; i++) {
			n=&ElemJ->node[i];
			nodeNum = n->nodeNum;
      		tempx += n->x;
      		tempy += n->y;
      		tempz += n->z;
    	}
    	ElemJ->centroid.x = tempx / nodesPerElem;
    	ElemJ->centroid.y = tempy / nodesPerElem;
    	ElemJ->centroid.z = tempz / nodesPerElem;
  	}
}
//=========================================================================
void setVectorFromMeshDisplacements(BasicMesh *mesh, double *vector)
{
/* Extract solution from node objects & store in single vector

This function cannot handle setting the dofs for dummy nodes. 
or in general, any node that has more than 3 dofs per node.
if the goal is to read in a displacement/forces file, 
use setVectorFromNodalValueFile(char* filename, double *vector) instead
*/
	

	int i,firstDof;
	int nn=mesh->numNodes;
	int numdofpn;
	Node *n,*d;

	for (i=0;i<nn;i++){
		d=&(mesh->displacements[i]);
		n=&(mesh->node[i]);
		firstDof=n->getFirstDof();
		numdofpn=n->getNumDof();
		vector[firstDof]	= d->x;
		if(numdofpn>1)
			vector[firstDof+1]	= d->y;
		if(numdofpn>2)
			vector[firstDof+2]	= d->z;
	}
}
//====================================================================
void setVectorFromNodalValueFile(char* filename, double *vector, int nn, NodeGroup &node)
{
/* Read solution file (e.g. u v w) and store in a single vector
It uses the NodeGroup to determine the number of dof/node.
A slower version could count the tokens on each line in the 
file to determine the number of dof/node
*/

//opening files like this could be unsafe when used within BETA because of the current working directory might not match etc.
//eventually, maybe the fucntion should be changed so that filemanager is used to open the files.
	ifstream is;
	is.open(filename);
	if(is.good()==false){
		cout << "error opening file in setVectorFromNodalValueFile()";
		exit(1);
	}

	char dummy[256];
	int ctr,i,j,dummyi;
	ctr=0;
	is >> dummy;
	int numdofpn;
	
	for(i=0;i<nn;i++){
		numdofpn=node[i].getNumDof();
		is >> dummyi;
		for(j=0;j<numdofpn;j++){
			is >> vector[ctr++];
		}
	}
}
//====================================================================
bool readMaterialGroups(istream  &fp, BasicMesh *mesh)
{
	int i,numElem = mesh->numElements,first,last,increment,group;
	ElementGroup *element = &mesh->element;
	BasicElement *e=0;
	
	char option[255];
	bool flag;
	while(2==2) {
		fp.getline(option, 255);
		if( COMPARE(option,"selectElementMaterial")==0 ||
			COMPARE(option,"setElementMaterial")==0 ||
			COMPARE(option,"readmaterialgroup()")==0 ) continue; //JV091003 - making different versions of the elemat file compatible.
		sscanf(option," %d %d %d %d ",&first,&last,&increment,&group);
		if(first==-1) {
			return true;
		}
		flag = (group < 0)?true:false;
		first-=mesh->elementOffset;
		last-=mesh->elementOffset;
		if(group > 12) {
			cerr << "Warning: Color Number should not exceed 12.\n";
			group = 12;
		}
		for (i=first;i<=last && i<numElem;i+=increment) {
			e=&(*element)[i];
			e->setMaterialGroup(group);
		}
	}
}
//====================================================================
bool readMaterialGroupsFromFile(const char *name, BasicMesh* mesh)
{
	ifstream fp(name);
	bool result=readMaterialGroups(fp, mesh);
	fp.close();
	return result;
}
//====================================================================
void ExtractSubMeshUsingMaterialGroups(BasicMesh &mesh, Array<int> &materialgroups, string  filenamePrefix)
{
// Create new mesh from just those materials in the array of material numbers.
// There are multiple files created. "filenamePrefix " is the 
// initial part of each filename.
	int numE=mesh.getNumActiveElements();
	int numN=mesh.getNumActiveNodes();
	int numG=materialgroups.getNumElements();
	int i,j;
	//deactivate the whole mesh
	mesh.SetActiveState(false);

	for (i=0;i<numE;i++){
		for (j=0;j<numG;j++){
			if (mesh.element[i].getMaterialNumber() == materialgroups[j]){
				mesh.element[i].SetActiveFlag(true);
			}

		}
	}
	mesh.assignNewNodeNum();
	mesh.PrintMesh(filenamePrefix, true);

	mesh.InitializeOutputSummary();
	mesh.OutputSummary();
	mesh.FinalizeSummary(filenamePrefix);

}
//====================================================================
void ExtractSubMeshUsingVolume(BasicMesh &mesh, Node &n1, Node &n2, string  filenamePrefix)
{// Nodes n1 and n2 are the "extreme" vertices of a rectangular parallelipiped

	//deactivate the whole mesh
	mesh.SetActiveState(false);

	//obtain list of elements in the volume
	ElementGroup* elist=getElementsInsideVolume(mesh.element, n1, n2);

	int i;
	//activate all the elements (and its nodes) in the list 
	for (i=0; i<elist->getNumElements(); i++) (*elist)[i].SetActiveFlag(true);
	
	mesh.assignNewNodeNum();
	mesh.assignNewElementNum();

	mesh.PrintMesh(filenamePrefix, true);

	mesh.InitializeOutputSummary();
	mesh.OutputSummary();
	mesh.FinalizeSummary(filenamePrefix);
	
	if(elist)	elist->clear();
	if(elist)	delete elist;
}
//====================================================================
bool writeElematFile(BasicMesh &mesh, string outputfilename)
{
	ofstream os;
	os.open(outputfilename.c_str());

	int numE=mesh.numElements;
	int i;

	int first=0,last=0, inc=1, matgroup=-1,oldmatgroup=-1;

	os << "readmaterialgroup()" << endl;
	BasicElement *e=0;

	int ei=0;
	for(i=0;i<numE;i++){
		e=&mesh.element[i];
		if(e->activeElementFlag==false)
			continue;
		matgroup = e->getMaterialNumber();
		if(ei==0)
			oldmatgroup=matgroup;
		if(matgroup != oldmatgroup){
			os << first << " " << last << " " << inc << " " << oldmatgroup << endl;
			first = ei;
			oldmatgroup = matgroup;
		}else{
			last=ei;
		}
		ei++;
	}
	os << first << " " << last << " " << inc << " " << matgroup << endl;
	os << "-1 0 0" << endl;
	os.close();
	return true;
}
//====================================================================
BasicElement* getElementContainingNodes(bool ListHas2Nodes, Node **n)
/* 
Picks out the first element that has the 2 (for linear) or 3 (for quadratic) specified nodes.
Exploits a "new" node method that has previously identified all elements attached to each node.

Assumptions: works for 1d, 2d and 3d BUT returns the FIRST element that it finds which contains the specified nodes.
there could be more than one element that contains the specified nodes
*/
{
	BasicElement *e1=0,*e2=0,*e3=0;
	Node *n1=n[0];
	Node *n2=n[1];
	Node *n3=n[2];

	if(ListHas2Nodes){
		for(int i=0; i<n1->numElementsAtNode; i++){
			e1=n1->eList[i];
			for(int j=0; j<n2->numElementsAtNode; j++){
				e2=n2->eList[j];
				if(e1==e2)
					return e1;
			}
		}
		return 0;
	}
	else {//ListHas3Nodes
		for(int i=0; i<n1->numElementsAtNode; i++){
			e1=n1->eList[i];
			for(int j=0; j<n2->numElementsAtNode; j++){
				e2=n2->eList[j];
				if(e1==e2)
					for(int k=0; k<n3->numElementsAtNode; k++){
						e3=n3->eList[k];
						if(e1==e3)
							return e1;
					}
			}
		}
		return 0;
	}
}
//======================================================================
double CalculateLocationOfValueAlongLine(bool isLinear, double *xCoor, double *values, double specifiedValue)
{/*
 Given the nodal values and coordinates of the nodes along a straight line, 
 calulate the coordinate of the point with the specified value.
 */

	double xi=0.0;
	double x=0.0;
	if(isLinear==false) //quadratic
	{
		if( fabs(2.*(values[0]-2*values[1]+values[2])) >= 1e-8){
			//do a check so u don't get imaginary roots!
			double xi1 = (values[0]-values[2] + sqrt(values[0]*values[0] - 2.*values[0]*values[2] 
				   +values[2]*values[2] + 8.*values[0]*specifiedValue - 8.*values[0]*values[1]
				   -16.*values[1]*specifiedValue + 16.*values[1]*values[1] + 8.*values[2]*specifiedValue
				   - 8.*values[2]*values[1]))/(2.*(values[0]-2.*values[1]+values[2]));

			double xi2 = (values[0]-values[2] - sqrt(values[0]*values[0] - 2.*values[0]*values[2] 
				   +values[2]*values[2] + 8.*values[0]*specifiedValue - 8.*values[0]*values[1]
				   -16.*values[1]*specifiedValue + 16.*values[1]*values[1] + 8.*values[2]*specifiedValue
				   - 8.*values[2]*values[1]))/(2.*(values[0]-2*values[1]+values[2]));

			if ( fabs(xi1) > (1+1e-8) ){xi=xi2;}
			else if( fabs(xi2) > (1+1e-8) ){xi=xi1;}
			else {
				cout<<"Specified value not in this domain."<<endl;
				exit(1);
			}

			x = 0.5*(xi*(xCoor[2]-xCoor[0])+(xCoor[2]+xCoor[0]));
		}else{//the three points are colinear
			xi = (2*specifiedValue - values[0] - values[2])/(values[2]-values[0]);

			if ( fabs(xi) > (1+1e-8) ){
				cout<<"Specified value not in this domain."<<endl;
				exit(1);
			}

			x = 0.5*(xi*(xCoor[2]-xCoor[0])+(xCoor[2]+xCoor[0]));
		}
	}

	else //is linear
	{
		xi = (2*specifiedValue - values[0] - values[1])/(values[1]-values[0]);

		if ( fabs(xi) > (1+1e-8) ){
			cout<<"Specified value not in this domain."<<endl;
			exit(1);
		}

		x = 0.5*(xi*(xCoor[1]-xCoor[0])+(xCoor[1]+xCoor[0]));
	}

	return x;
}

bool findLocationOfValueAlongLine(bool isLinear, Node **n, double *values, int direction, double specifiedValue, Node &locatedNode)
{
	int numNodes;
	if(isLinear) numNodes=2;
	else numNodes=3;

	const double TOLERANCE=1e-5;

	if( (specifiedValue < values[0]) || (specifiedValue > values[numNodes-1]) ) return false;
	if( fabs(specifiedValue - values[0]) < TOLERANCE )			{
		locatedNode[direction-1]=(*n[0])(direction-1);
		return true;
	}
	if( fabs(specifiedValue - values[numNodes-1]) < TOLERANCE )	{
		locatedNode[direction-1]=(*n[numNodes-1])(direction-1);
		return true;
	}

	double xCoor[3];
	xCoor[0]=(*n[0])(direction-1);
	xCoor[1]=(*n[1])(direction-1);
	if(!isLinear)
		xCoor[2]=(*n[2])(direction-1);

	double location=CalculateLocationOfValueAlongLine(isLinear, xCoor, values, specifiedValue);
	locatedNode[direction-1]=location;
	return true;
}
//======================================================================
int getNumberOfElementsFromMesh(char* filename)
{
	ifstream is;
	is.open(filename);
	int numNodes,numElements,numDims;
	is >> numNodes >> numElements >> numDims;
	is.close();	
	return numElements;
}
//======================================================================
void GetMinMax(NodeGroup *nlist, Node *Min, Node *Max, bool initializeMinMax)
{
	NodeGroup &node=*nlist;
	Node &min=*Min;
	Node &max=*Max;

	int numNodes=node.getNumElements();
#define MINMAX(a,b,c) {if(a > b) b = a; else if(a < c) c = a;}
	{
		if(initializeMinMax){
			max.x = min.x = node[0].x;
			max.y = min.y = node[0].y;
			max.z = min.z = node[0].z;
		}
		for(int i=1;i<numNodes;i++) {
			MINMAX(node[i].x,max.x,min.x);
			MINMAX(node[i].y,max.y,min.y);
			MINMAX(node[i].z,max.z,min.z);
		}
	}
#undef MINMAX
}
//======================================================================
void GetMinMax(ElementGroup *elist, Node *Min, Node *Max, bool initializeMinMax)
{
	for(int i=0; i<elist->getNumElements(); i++){
		BasicElement* e=&(*elist)[i];
		GetMinMax(&(e->node), Min, Max, initializeMinMax);
		if(initializeMinMax) initializeMinMax=false;
	}
}
//======================================================================
