#include "stdafx.h"

#include "mesh/BasicMesh.hpp"
#include "utility/utility.h"
#include "mesh/NodeSet.hpp"

double findNodeTOLERANCE=1e-8;

//========================================================================
//function declarations

//Misc
int 		getNodeIndex(NodeGroup *nlist, double x, double y, double z);
bool		NodeInsideVolume(Node &n, Node &n1, Node &n2);
bool 		ElementInsideVolume(BasicElement &e, Node &n1, Node &n2);
void		sortElementGroupbyCoord(ElementGroup &eList,int direction);

//Node groups
// coordum = 1,2, or 3
NodeGroup*  getNodesOnPlane(BasicMesh *mesh,  int coordnum, double coordval);
NodeGroup*  getNodesOnPlane(NodeGroup *nlist, int coordnum, double coordval);
NodeGroup* createListOfCommonNodes(NodeGroup * refList,NodeGroup *compList);
			
//Element groups
ElementGroup*  createElementGroupFromNodeSet(const NodeSet &nodes);
ElementGroup*  createElementGroupUsingMaterialNumber(ElementGroup* eList, int matNum);
ElementGroup*  getElementsInsideVolume(ElementGroup &elist, Node &n1, Node &n2);


//=========================================================
ElementGroup* getElementGroup(char* name, Array<ElementGroup> *GroupList)
{
	int i;
	int numGroups=0;
	ElementGroup* group=0;
	numGroups=GroupList->getNumElements();
	for(i=0;i<numGroups;i++){
		group=&(*GroupList)[i];
		if(COMPARE(group->GetName(), name)==0) 
			return group;
	}
	return 0;
}
//=========================================================
ElementGroup* CreateElementGroupFromElematFile(string elematfilename, BasicMesh *mesh, FileManager *filemanager)
{
	char *localTokenList[20];
	int  numberOfTokens;
	ifstream *inStream=filemanager->OpenInputStream(elematfilename);
	if( getLineAndTokenize(inStream,"end",localTokenList, numberOfTokens)==1)
		EXIT_BETA("No information in elemat file");

	ElementGroup *elementgroup= new ElementGroup;
	elementgroup->reserve(mesh->numElements);

	int i,first,last,increment;

	while (2==2) {				
		if( getLineAndTokenize(inStream,"end",localTokenList, numberOfTokens)==1) {
			break;
		}
		if(numberOfTokens >= 2 && COMPARE(localTokenList[0],"all")==0 ){
			first=0;
			last=mesh->numElements-1;
			increment=1;
		}else if(numberOfTokens >= 4){
			first=atoi(localTokenList[0]);
			last=atoi(localTokenList[1]);
			increment=atoi(localTokenList[2]);
		}else{
			first=atoi(localTokenList[0]);
		}
		if(first<0) break;

		//Error checking...
		if(first < 0 || first    > mesh->numElements || first > last
					 || last     > mesh->numElements || increment<1 ){
			EXIT_BETA("Fatal error...input is bad\nFirst, last, increment="<<first<<" "<<last<<" "<<increment);
		}

		BasicElement *e=0;
		BETA_OUT << "\tAdding elements "<< first <<" to "<<last << " by " << increment << '\n';				
		for(i=first;i<=last;i+=increment) {
			e=&mesh->element[i];
			elementgroup->add(*e);
		}
	}
	if(elementgroup->getNumElements()==0) {
		delete elementgroup;
		elementgroup=0;
	}
	return elementgroup;
}
//=========================================================
NodeGroup*  getNodesOnPlane(BasicMesh *mesh, int coordnum, double coordval)
{
	return getNodesOnPlane(&mesh->node, coordnum, coordval);
}
//=========================================================
NodeGroup*  getNodesOnPlane(NodeGroup *nlist, int coordnum, double coordval)
{//coordnum starts from 1 -> 3
	int i;
	int numNodes=(int)nlist->getNumElements();
	NodeGroup* list=0;

	list= new NodeGroup;
	list->reserve(numNodes);

	Node* node=0;
	for(i=0;i<numNodes;i++)
	{
		node=&(*nlist)[i];
		if ( fabs((*node)(coordnum-1) - coordval) < findNodeTOLERANCE )
			list->add(*node);
	}
	if (list->empty() ){
		delete list;
		list=0;
		cout <<"could not find any nodes with coordnum =" << coordnum << " and coordval= " << DOUBLE_FORMAT << coordval <<  endl;
	}
	return list;
}
//=========================================================
void sortElementGroupbyCoord(ElementGroup &eList,int direction)
{
// If this is ever used on large lists, use lib routine qsort
	int numElements = eList.getNumElements();
	for(int i=0;i<numElements-1;i++){
		for(int j=i+1;j<numElements;j++){
			double vi=eList[i].centroid(direction-1);
			double vj=eList[j].centroid(direction-1);
			if(vj<vi)
				eList.swap(i,j);
		}
	}
}
//=========================================================
int getNodeIndex(NodeGroup *nlist, double x, double y, double z)
{
	int i;
	int numNodes=nlist->getNumElements();
	int index=-1;
	Node *n;
	for (i=0; i<numNodes; i++) {
		n=&(*nlist)[i];

		if ( (fabs(x-n->x) < findNodeTOLERANCE) && 
			 (fabs(y-n->y) < findNodeTOLERANCE) &&
			 (fabs(z-n->z) < findNodeTOLERANCE)  ) 
		{index=i;	break; 	}
	}
	if (index == -1){
		cout <<"could not find node with coordinates(" << DOUBLE_FORMAT << x << ',' << y << ',' << z << ')' <<  endl;
		//exit(1);
	}
	return index;
}
//=========================================================
bool NodeInsideVolume(Node &n, Node &n1, Node &n2)
{// Nodes n1 and n2 are "extreme" vertices of rectangular parallelipiped
	if(n2.x > n1.x){
		if ( n.x < (n1.x - findNodeTOLERANCE) ) return false;
		if ( n.x > (n2.x + findNodeTOLERANCE) ) return false;
	}else{
		if ( n.x < (n2.x - findNodeTOLERANCE) ) return false;
		if ( n.x > (n1.x + findNodeTOLERANCE) ) return false;
	}

	if(n2.y > n1.y){
		if ( n.y < (n1.y - findNodeTOLERANCE) ) return false;
		if ( n.y > (n2.y + findNodeTOLERANCE) ) return false;
	}else{
		if ( n.y < (n2.y - findNodeTOLERANCE) ) return false;
		if ( n.y > (n1.y + findNodeTOLERANCE) ) return false;
	}

	if(n2.z > n1.z){
		if ( n.z < (n1.z - findNodeTOLERANCE) ) return false;
		if ( n.z > (n2.z + findNodeTOLERANCE) ) return false;
	}else{
		if ( n.z < (n2.z - findNodeTOLERANCE) ) return false;
		if ( n.z > (n1.z + findNodeTOLERANCE) ) return false;
	}

	return true;
}
//=========================================================
bool ElementInsideVolume(BasicElement &e, Node &n1, Node &n2)
{// Nodes n1 and n2 are "extreme" vertices of rectangular parallelipiped
	int i;
	for(i=0; i<e.numNodesPerElement; i++){
		if (NodeInsideVolume(e.node[i], n1, n2) == false)
			return false;
	}
	return true;
}
//=========================================================
ElementGroup*  getElementsInsideVolume(ElementGroup &elist, Node &n1, Node &n2)
{
	int i;
	int numElements=elist.getNumElements();
	ElementGroup* list=0;

	list= new ElementGroup;
	list->reserve(numElements);

	BasicElement* e=0;
	for (i=0;i<numElements;i++){
		if(ElementInsideVolume(elist[i], n1, n2) == true)
			list->add(elist[i]);
	}

	return list;
}
//=========================================================
ElementGroup* createElementGroupUsingMaterialNumber(ElementGroup* eList,int matNum)
{
	ElementGroup *eListMat=NULL;
	eListMat = new ElementGroup;
	
	int numEl = eList->getNumElements();
    eListMat->setCapacity(numEl);
	//always remember to clear the array if you don't want its elements to be deleted
	//allocating more space than needed... but doing this will avoid the need for reallocs.
	//eListMat.setCapacity(mesh->numElements);

	BasicElement *e;

	for(int i=0;i<numEl;i++){
		e=&(*eList)[i];
		if(e->getMaterialNumber() == matNum)
		  eListMat->add(*e);
	}

return eListMat;
}
//========================================================================
NodeGroup*  createNodeGroupFromElementGroup(int numNodes,ElementGroup *eList)
{
	//does not reduce duplicate nodes shared by elements
	int numElements = eList->getNumElements();

	NodeGroup* nList;
	nList = new NodeGroup;
	nList->reserve(numNodes);


	BasicElement* e=0;
	for (int i=0;i<numElements;i++){
		e=&(*eList)[i];
		for(int j=0; j<e->getNumNodesPerElement();j++)
		{
			nList->add(e->node[j]);
		}
	}

	return nList;
	
}
//========================================================================
ElementGroup*  createElementGroupFromNodeSet(const NodeSet &nodes)
{
	int numNodes = nodes.getNumNodes();
	int maxElements=0;
	for(int i=0; i< numNodes; i++)
		maxElements += nodes[i]->numElementsAtNode;


	ElementGroup* eList;
	eList = new ElementGroup;
	eList->reserve(maxElements);

	for(int i=0; i< numNodes; i++)
	{
		for(int j=0; j< nodes[i]->numElementsAtNode; j++){
			if( eList->isInByAddress(nodes[i]->eList[j]) == -1)
				eList->add(*(nodes[i]->eList[j]));
		}
	}

	return eList;
}
//========================================================================
NodeGroup* createListOfCommonNodes(NodeGroup * refList,NodeGroup *compList)
{
int refSize = refList->getNumElements();
int compSize = compList->getNumElements();

NodeGroup *nList;
nList = new NodeGroup;
nList->reserve(compSize);

Node *refN;
Node *compN;

for(int i=0; i<compSize;i++)
{
	compN=&(*compList)[i];
   for(int j=0;j<refSize;j++)
   {
	   refN=&(*refList)[j];
	   if(compN->getNodeNum() == refN->getNodeNum()){
			nList->add(*refN);
			break;
	   }
   }
}

return nList;
}
//=========================================================
