#include "stdafx.h"

#include "utility/utility.h"
#include "utility/FileManager.hpp"
#include "MultiFieldMacroElement.hpp"

#include "utility/excepts.hpp"
#include "math/matrix.hpp"
#include "math/equation/equation.hpp"
#include "factory/Factory.hpp"
#include "mesh/NodeSet.hpp"

extern int verboseFlag;
extern Factory	*factory;

void shape3D(int numberOfInterp, double xi,  double eta, double zeta, double * S);
bool readScript(Factory *factory, BasicModel* &modelObject, string filename, ifstream *& inStream);
//=====================================================================
MultiFieldMacroElement::MultiFieldMacroElement(void)
{
}
//=====================================================================
MultiFieldMacroElement::~MultiFieldMacroElement()
{
}
//=====================================================================
Matrix * calcTransformationMatrix(BasicMesh *mesh, BoundaryNodeList *bn)
{
	int numSlaves=bn->numNodesInList; //this includes the masters as well
	int numMasters;
	if (mesh->numDims==3) 
		numMasters=20; // handle other cases later
	else 
		EXIT_BETA("calcTransformationMatrix() currently handles only 3D meshes");
	Matrix *T=new Matrix(numSlaves, numMasters);

//	InterpolationData       interpData;
//	interpData.numberOfInterp = numMasters;

	InterpAndSolDerivatives *interpData=new InterpAndSolDerivatives(numMasters,numMasters); interpData->allocate();

	int nodenum;
	int i,j;
	for (i=0;i<numSlaves;i++){
		nodenum = bn->nList[i]->nodeNum;
		shape3D(interpData->numberOfInterp, mesh->node[nodenum].x, mesh->node[nodenum].y, mesh->node[nodenum].z, interpData->S);
		for (j=0;j<numMasters;j++){
			(*T)(i,j) = interpData->S[j];
		}
	}
	delete interpData;
	return T;
}
//=====================================================================
void convertToGlobalCoordinates(double xi, double eta, double zeta, InterpAndSolDerivatives &interpData, double &x, double &y, double &z)
{
	shape3D(interpData.numberOfInterp, xi, eta, zeta, interpData.S);

	x=y=z=0.0;

	for(int kk=0; kk < interpData.numberOfInterp; kk++){
		x	+= interpData.S[kk] * interpData.xCoor[kk];
		y	+= interpData.S[kk] * interpData.yCoor[kk];
		z	+= interpData.S[kk] * interpData.zCoor[kk];
	}
}
//=====================================================================
void transformMesh(BasicMesh *mesh, BasicElement *templateElement)
// the mesh has to be in normalized coordinates
{
	int i;
	double x,y,z;
	Node *n;
//	InterpolationData       interpData;
//	interpData.numberOfInterp = templateElement->numNodesPerElement;

	InterpAndSolDerivatives *interpData=new InterpAndSolDerivatives(templateElement->numNodesPerElement,templateElement->numNodesPerElement); interpData->allocate();

	templateElement->extractNodalCoordinates(interpData->xCoor,interpData->yCoor, interpData->zCoor);

	for (i=0; i<mesh->numNodes; i++){
		n=&(mesh->node[i]);
		convertToGlobalCoordinates(n->x, n->y, n->z, *interpData, x,y,z);
		n->x=x;
		n->y=y;
		n->z=z;
	}
	delete interpData;
}
//=====================================================================
/*
Group* getReducedNodesList(BasicMesh *mesh, BasicElement *macroElement)
{
	Group *g=0;
	g=new Group("ReducedNodesList",NodeGroupType, 1);

	GroupItem *aItem;
	Node *n;
	NodeGroup* list=0;

	int i;
	int numNodes=macroElement->numNodesPerElement;
	g->items.reserve(numNodes);
	for (i=0; i< numNodes; i++){
		list = mesh->getNodesUsingCoordinates(&macroElement->node[i]);
		if (list==0)
			EXIT_BETA("did not find node !!!");
		n=&(*list)[0];
		list->clear();
		delete list;
		aItem=new GroupItem(n->nodeNum);
		g->items.add(*aItem);
	}
	return g;
}
*/
//=====================================================================
NodeSet getReducedNodesList(BasicMesh *mesh, BasicElement *macroElement)
{
    //This appears to eliminate duplicate nodes.  There are faster ways to do this.
	NodeSet g(mesh->getNodeSearchTolerance());
	Node *n=0;
	int numNodes=macroElement->numNodesPerElement;
	for(int i=0; i<numNodes; i++){
		NodeSet list(mesh->getNodeSearchTolerance());
		list = mesh->getNodesUsingCoordinates(&macroElement->node[i]);
		if (list.getNumNodes()==0)
			EXIT_BETA("did not find node !!!");
		n=list[0];
		g.addNode(n);
	}
	return g;
}
//=====================================================================
void MultiFieldMacroElement::solveForStiffnessMatrixColumn(int currentDof, int currentNode, int masterNodeIndex, int numMasterNodes, int numReducedDofs, BoundaryNodeList *bn, int *reducedDofList, Matrix *T, LargeMatrix *gK)
//this version is for the multifield macro element
//the column is written to the internal force vector.
{
	bool haveToAssemble;
	if(currentDof == reducedDofList[0]) 
		haveToAssemble=true;
	else
		haveToAssemble=false;

	int i,j;

	double maxResidual;
	int    dofForMaxResidual;
	double * mappedResidualVector = new double[model->totalNumDof];

	if(haveToAssemble){
		model->equations.cleanRestraints();
		model->equations.zeroLoadVector();
		model->equations.zeroResultantVector();
		model->equations.zeroSolution();
		model->equations.zeroMatrix();
	}else{
		model->equations.zeroPseudoLoadVector();
		model->equations.zeroLoadVector();
		model->equations.zeroSolution();
	}

	BasicMesh* mesh;
	//we are assuming for now that there is no hierarchy of models 
	//for the subdomain mesh. but this is something that can be added on later.
	mesh=this->model->mesh;
	int numSlaves=bn->numNodesInList;

	int nodenum,numDofAtNode,rdof;
	Node *n;

	//restrain all the slaves (including masters)
	for (i=0;i<numSlaves;i++){
		nodenum = bn->nList[i]->nodeNum;
		n=&(mesh->node[nodenum]);
		numDofAtNode = n->getNumDof();
		for( j=0; j<numDofAtNode; j++){
			rdof=n->getFirstDof()+j;
			model->equations.setRestraint(rdof,true);
		}
	}

	n=&(mesh->node[currentNode]);
	rdof=n->getFirstDof();
	int idof = currentDof - rdof;
	//now,set the values of the known slaves
	double disp;
	for (i=0;i<numSlaves;i++){
		nodenum = bn->nList[i]->nodeNum;
		n=&(mesh->node[nodenum]);
		rdof=n->getFirstDof()+idof;
		disp=(*T)(i,masterNodeIndex);
		model->equations.setSolution(disp,rdof); 
		if(!haveToAssemble){
			model->equations.updatePseudoLoadVectorUsingGK(disp,rdof,gK); 
		}
//		if(disp==1.0)
//			cout << "Say so !! " << endl;
	}
	//...finished setting constraints.

	if(haveToAssemble){
//		setMaterial();//JV071709
		model->assembler->assembleK(NULL);
		displayTime("Time to assembleKandF",BETA_OUT);
	}

	//FormWriter fr(BETA_OUT);
	//equations.printMatrix(fr);
	//equations.printPseudoLoadVector(fr);

	
	//Now solving the set of equations..
	BETA_OUT<<"Start solving equations...("<<model->totalNumDof<< " Equations)." <<endl;
	displayTime("",BETA_OUT);
	model->solution = model->equations.solve(); 
//	solution = equations.getSolution(); //just for testing

	displayTime("Time to solve",BETA_OUT);
	BETA_OUT<<"Done solving equations."<<endl;

	//FLEX94 uses a quicker method: it calculates the forces from only the boundary elements
	//right now, i am using the general function available to calculate the global forces
	//we might want to switch over to the idea used in flex94
	model->getGlobalForces(model->internalForce, model->solution,
					mappedResidualVector,
					maxResidual,dofForMaxResidual);

	displayTime("Time to calculate and print stresses",BETA_OUT);

	delete [] mappedResidualVector;
}
//==============================================================
/*
void MultiFieldMacroElement::getReactionForces(double *reaction, double *internalForce, Matrix *T, Group *masterNodeList, BoundaryNodeList *bn, int numNodes)
{
	int i,j,k,ii,m;
	double sum=0.0;
	int numSlaves = bn->numNodesInList;
	int nodenum,numDofAtNode,rdof;
	Node *n;

	BasicMesh* mesh;
	//we are assuming for now that there is no hierarchy of models 
	//for the subdomain mesh. but this is something that can be added on later.
	mesh=this->model->mesh;

	ii=0;
	int masterNodeIndex=0;
	G_ITERATOR(k,j,masterNodeList->items)
		n = &(mesh->node[j]);
		numDofAtNode = n->getNumDof();
		for( i=0; i<numDofAtNode; i++){
			sum=0.0;
			for (m=0;m<numSlaves;m++){
				nodenum = bn->nList[m]->nodeNum;
				n=&(mesh->node[nodenum]);
				rdof=n->getFirstDof()+i;
				sum += internalForce[rdof] * (*T)(m,masterNodeIndex);
			}
			reaction[ii]=sum;
			ii++;
		}
		masterNodeIndex++;
	G_END
}
*/
//=================================================
void MultiFieldMacroElement::getReactionForces(double *reaction, double *internalForce, Matrix *T, NodeSet &masterNodeList, BoundaryNodeList *bn, int numNodes)
{
//	int i,j,k,ii,m;
	double sum=0.0;
	int numSlaves = bn->numNodesInList;
	int nodenum,numDofAtNode,rdof;
	Node *n;

	BasicMesh* mesh;
	//we are assuming for now that there is no hierarchy of models 
	//for the subdomain mesh. but this is something that can be added on later.
	mesh=this->model->mesh;

	int ii=0;
	int masterNodeIndex=0;
	for(int j=0;j<masterNodeList.getNumNodes(); j++){
		n = masterNodeList[j];
		numDofAtNode = n->getNumDof();
		for(int i=0; i<numDofAtNode; i++){
			sum=0.0;
			for (int m=0;m<numSlaves;m++){
				nodenum = bn->nList[m]->nodeNum;
				n=&(mesh->node[nodenum]);
				rdof=n->getFirstDof()+i;
				sum += internalForce[rdof] * (*T)(m,masterNodeIndex);
			}
			reaction[ii]=sum;
			ii++;
		}
		masterNodeIndex++;
	}
}
//=================================================
void MultiFieldMacroElement::CalculateK(double *uGlobal,string filename)
{
	BasicMesh *SubDomainMesh=model->mesh;
	BoundaryNodeList *bn=0;
	bn=SubDomainMesh->findBoundaryNodes();
	bn->printBoundaryNodeList("SubstructureNodeList.txt");
	int numSlaves = bn->numNodesInList;

	Matrix *T;
	T=calcTransformationMatrix(SubDomainMesh,bn);
	T->print("Transformation Matrix", bn->numNodesInList, 20, &BETA_OUT); 

	transformMesh(SubDomainMesh, modelOwner); //now the mesh is in global coordinates //JV071609: modelOwner is also 'this' in this particular case.
	string TMeshName;
	TMeshName = "Transformed_" + SubDomainMesh->getMeshFileName();
	ofstream os(TMeshName.c_str());
	SubDomainMesh->PrintMesh(&os,true);
  SubDomainMesh->PrintMesh(filename,false);
	os.close();


	NodeSet	reducedNodeList;
	reducedNodeList=getReducedNodesList(SubDomainMesh, modelOwner);

	Node *n;
	int numDofAtNode;
	int currentDof;
	int numReducedDofs=0;
	int numNodes=0;

	//calculating numReducedDofs 

	for(int j=0;j<reducedNodeList.getNumNodes(); j++){
		n = reducedNodeList[j];
		numReducedDofs += n->getNumDof();
		numNodes++;
	}


	int ii;
	int *reducedDofList;
	reducedDofList= new int[numReducedDofs];
	MasterKe=new Matrix(numReducedDofs,numReducedDofs);
	//the matrix is already initialized to zero.


	ii=0;
	for(int j=0;j<reducedNodeList.getNumNodes(); j++){
		n = reducedNodeList[j];
		numDofAtNode = n->getNumDof();
		for(int i=0; i<numDofAtNode; i++){
			currentDof=n->getFirstDof()+i;
			reducedDofList[ii]=currentDof;
			ii++;
		}
	}

	double *reaction;
	reaction= new double[numReducedDofs];

	//assembling a global K without any restraints. 
	//therefore we have the actual assembles global K
	//and then create a copy so that it can be used to solve consecutive matrices
	//without having to assmeble them from scratch again.
//    setMaterial();//JV071709 materials are already set, right?
	model->assembler->assembleK(NULL);

	LargeMatrix *gK;
	gK=model->equations.createCopyOfMatrix();
	

	//FormWriter fr(BETA_OUT);
	//equations.printMatrix(fr);
	//mat->printMatrix("Kg mat = ",&OS);
	//gK.printMatrix("Kg (copy) = ",&OS);

	ii=0;
	int masterNodeIndex=0;
	for(int j=0;j<reducedNodeList.getNumNodes(); j++){
		n = reducedNodeList[j];
		numDofAtNode = n->getNumDof();
		for(int i=0; i<numDofAtNode; i++){
			currentDof=n->getFirstDof()+i;
			solveForStiffnessMatrixColumn(currentDof, n->nodeNum, masterNodeIndex, numNodes, numReducedDofs, bn, reducedDofList, T, gK);
			cout << "solved for dof number: " << ii << endl;
			//Now, transforming forces :
			getReactionForces(reaction, model->internalForce, T, reducedNodeList, bn, numNodes);
			for (int m=0;m<numReducedDofs;m++){
				(*MasterKe)[m][ii] = reaction[m];
			}	
			ii++;
		}
		masterNodeIndex++;
	}

	MasterKe->print("MasterKe Matrix",numReducedDofs,numReducedDofs, &BETA_OUT);
	MasterKe->printInMaple("Kg MF Maple= ",numReducedDofs,numReducedDofs,&BETA_OUT);

	delete [] reducedDofList;
	delete [] reaction;
	delete bn;
	delete T;
	delete gK;


}
//==========================================================
void MultiFieldMacroElement::update(char * action, double *uGlobal)
{
	if(strstr(action,"K") !=NULL) {
		if(model)
			Ke=modelOwner->MasterKe;
		else
			EXIT_BETA("you cannot call MultiFieldMacroElement::update() before the macroelement has been defined.");
	}else if(strstr(action,"I") !=NULL) {
		for(int i=0;i<numDof;i++) initialFe[i] = 0.;
	}else{
		cout << "calling MultiFieldMacroElement::update() with '"<< action << "'"<< endl;
		//EXIT_BETA("MultiFieldMacroElement::update() can currently calculate only the stiffness matrix");
    cout<<"WARNING: Incorrect stress output for macroelements...see SingleFieldMacroElement.cpp"<<endl; //ksm6-29-10
	}

}
//==========================================================
void MultiFieldMacroElement::readSpecialCommand(istream &inStream, ElementGroup *element, int numElements, char * command)
{
BETA_OUT<<"MultiFieldMacroElement::readSpecialCommand"<<endl;

if (COMPARE(command, "DefineMacroElementGroups")==0 ) {
	//calculate masterKe using first element in the ElementGroup
	MultiFieldMacroElement* master=(MultiFieldMacroElement*)(&(*element)[0]);

	BasicModel *modelObject =0; 
	ifstream *subdomainstream=0;
	string filename;
	inStream >> filename;
 
	readScript(factory, modelObject, filename,subdomainstream);
	displayTime("Time to finish running Subdomain script", *modelObject->filemanager->outStream);
	if(subdomainstream) delete subdomainstream;

	master->model=modelObject;
	master->modelOwner=master;
	modelObject->createAssembler();
	modelObject->allocateForAnalysis(BETA_ALLOCATE_EQUATIONS);
	master->setPointers(modelObject->assembler->threadWorkspaceList[0].eWorkspace);//using the subdomain model's workspace to calculate the master Ke
	master->calculateDofList();
	master->CalculateK(NULL,filename);

	//linking all the other elements to the subdomain model as well
	for(int i=1;i<numElements;i++) { 
		MultiFieldMacroElement* e=(MultiFieldMacroElement*)(&(*element)[i]);
		e->model=modelObject;
		e->modelOwner=master;
	}
	return;
}
// If no match, try super class
BETA_OUT<<"No match... try super class"<<endl;
MacroElement::readSpecialCommand(inStream,element,numElements, command);
}//readSpecialCommand
//==========================================================
