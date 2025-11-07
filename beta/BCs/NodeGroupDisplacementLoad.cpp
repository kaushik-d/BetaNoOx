#include "stdafx.h"

#include <string>
#include "NodeGroupDisplacementLoad.hpp"  
#include "models/BasicModel.hpp"
#include "mesh/NodeSet.hpp"

//==========================================================
NodeGroupDisplacementLoad::NodeGroupDisplacementLoad()
{
	WhoClass("NodeGroupDisplacementLoad");
	setType("Displacement");
	direction=-1;
	disp=0;
}
//==========================================================
bool NodeGroupDisplacementLoad::read(istream * inStream)
{
	char *localTokenList[20];
	int  numberOfTokens;

	if( getLineAndTokenize(inStream,"exitNodeGroupDisplacementLoad",localTokenList, numberOfTokens)==1) {
		return false;
	}

	if(numberOfTokens < 3){
		FatalError("Error in NodeGroupDisplacementLoad::read(istream * inStream)");
	}

	strcpy(NodeGroupName,localTokenList[0]);

	//numTokens has to be 3 or larger
	if(numberOfTokens >= 3){
		disp = atof(localTokenList[1]);
		direction = atoi(localTokenList[2]);
		if(direction<1)
			FatalError("Error in DisplacementLoad::read - bad direction number");
	}else{
		FatalError("Error in NodeGroupDisplacementLoad::read(istream * inStream)");
	}

	//if 4th token is "incremental", then set flag
	if(numberOfTokens == 4){
		string s="incremental";
		if (s==localTokenList[3])
			setIncFlag(true);
	}

	//this load is *really* setting a 'constraint', 
	//so you need to set a restraint if you want to take advantage of a reduced bandwidth
	const NodeSet nlist = model->mesh->getNodeSet(localTokenList[0]);
	Node *n=0;
	int globalDof=-1;
	int numNodes=nlist.getNumNodes();
	for(int i=0;i<numNodes;i++){
		n=nlist[i];
		if(n->numDof <= 0){
			BETA_OUT << "Node has no dofs. cannot apply displacement to this node - node num:" << n->getNodeNum() << endl;
		}else{
			globalDof = (*n).getFirstDof()+direction-1;
			model->equations.setRestraint(globalDof,true);
		}
	}

	BETA_OUT<<"NodeSet : "<<NodeGroupName <<endl;
	BETA_OUT<<"disp : "<<disp<<endl;
	BETA_OUT<<"direction : "<<direction<<endl;
	print(BETA_OUT);
	return true;
}
//==========================================================
void NodeGroupDisplacementLoad::print(ostream &ostrm)
{
//	ostrm << "Displacement Load - " << nodeNum << '\t' << disp << '\t' << 
//			direction << endl;
}
//==========================================================
void NodeGroupDisplacementLoad::applyTo(Equations &equations, NodeGroup *node, int numNodes, BasicModel *model)
{ 
	WhoMethod("applyTo");
	if(model==0) {
		cout << "model not passed as parameter to NodeGroupDisplacementLoad::applyTo(). fix the source code" << endl;
		exit(1);
	}

	const NodeSet nlist=model->mesh->getNodeSet(NodeGroupName);

	Node *n=0;
	int globalDof=-1;
	numNodes=nlist.getNumNodes();
	for(int i=0;i<numNodes;i++){
		n=nlist[i];
		if(n->numDof <= 0){
			BETA_OUT << "Node has no dofs. cannot apply displacement to this node - node num:" << n->getNodeNum() << endl;
		}else{
			globalDof = (*n).getFirstDof()+direction-1;
//			equations.setRestraint(globalDof,true); //shouldn't really have to do this since it is already been restrained during the read
			equations.setSolution(rate*disp,globalDof);  //xt:rate=1.0 by default;
		}
	}

}
//==========================================================
