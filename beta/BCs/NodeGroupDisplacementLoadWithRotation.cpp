#include "stdafx.h"

#include <string>
#include "NodeGroupDisplacementLoadWithRotation.hpp"  
#include "models/BasicModel.hpp"
#include "mesh/NodeSet.hpp"
NodeGroup* getNodeGroup(char* name, Array<NodeGroup> *GroupList);

//==========================================================
NodeGroupDisplacementLoadWithRotation::NodeGroupDisplacementLoadWithRotation()
{
	WhoClass("NodeGroupDisplacementLoadWithRotation");
	setType("Displacement");
	variableCoefficient=0;
	variableDirection=-1;
	direction=-1;
	disp=0;
}
//==========================================================
bool NodeGroupDisplacementLoadWithRotation::read(istream * inStream)
{
	char *localTokenList[20];
	int  numberOfTokens;

	if( getLineAndTokenize(inStream,"exitNodeGroupDisplacementLoadWithRotation",localTokenList, numberOfTokens)==1) {
		return false;
	}

	if(numberOfTokens < 5){ //modified by Sunil
		FatalError("Error in NodeGroupDisplacementLoadWithRotation::read(istream * inStream)");
	}

	strcpy(NodeGroupName,localTokenList[0]);

	//numTokens has to be 3 or larger
	if(numberOfTokens >= 5){ //modified by Sunil
		variableCoefficient = atof(localTokenList[1]); //modified by Sunil
		variableDirection = atoi(localTokenList[2]); //modified by Sunil
		disp = atof(localTokenList[3]); //modified by Sunil
		direction = atoi(localTokenList[4]); //modified by Sunil

		if(direction<1)
			FatalError("Error in DisplacementLoad::read - bad direction number");
	}else{
		FatalError("Error in NodeGroupDisplacementLoadWithRotation::read(istream * inStream)");
	}

	//if 4th token is "incremental", then set flag
	if(numberOfTokens == 6){ //modified by Sunil
		string s="incremental";
		if (s==localTokenList[5]) //modified by Sunil
			setIncFlag(true);
	}

	//this load is *really* setting a 'constraint', 
	//so you need to set a restraint if you want to take advantage of a reduced bandwidth
	//NodeGroup *nlist=getNodeGroup(localTokenList[0], &model->mesh->NodeGroupsList);
	const NodeSet nlist=model->mesh->getNodeSet(NodeGroupName);
	//if(nlist==0) FatalError("NodeGroupDisplacementLoadWithRotation : NodeGroup not found!\n");
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

	BETA_OUT<<"NodeGroup : "<<NodeGroupName <<endl;
	BETA_OUT<<"disp : "<<disp<<endl;
	BETA_OUT<<"direction : "<<direction<<endl;
	print(BETA_OUT);
	return true;
}
//==========================================================
void NodeGroupDisplacementLoadWithRotation::print(ostream &ostrm)
{
//	ostrm << "Displacement Load - " << nodeNum << '\t' << disp << '\t' << 
//			direction << endl;
}
//==========================================================
void NodeGroupDisplacementLoadWithRotation::applyTo(Equations &equations, NodeGroup *node, int numNodes, BasicModel *model)
{ 
	WhoMethod("applyTo");
	if(model==0) {
		cout << "model not passed as parameter to NodeGroupDisplacementLoad::applyTo(). fix the source code" << endl;
		exit(1);
	}

	const NodeSet nlist = model->mesh->getNodeSet(NodeGroupName);
	//if(nlist==0) FatalError("NodeGroupDisplacementLoad : NodeGroup not found!\n");

	Node *n=0;
	int globalDof=-1;
	double coordinate=0;
	numNodes=nlist.getNumNodes();
	for(int i=0;i<numNodes;i++){
		n=nlist[i];
		//added by Sunil ---------------------------------------
		if(variableDirection == 1)
			coordinate = (*n).x1;
		else if(variableDirection == 2)
			coordinate = (*n).x2;
		else if(variableDirection == 3)
			coordinate = (*n).x3;
		else
			FatalError("NodeGroupDisplacementLoad : invalid coordinate!\n");
		//------------------------------------------------------
		if(n->numDof <= 0){
			BETA_OUT << "Node has no dofs. cannot apply displacement to this node - node num:" << n->getNodeNum() << endl;
		}else{
			globalDof = (*n).getFirstDof()+direction-1;
//			equations.setRestraint(globalDof,true); //shouldn't really have to do this since it is already been restrained during the read
//			equations.setSolution(rate*disp,globalDof);  //xt:rate=1.0 by default;
			equations.setSolution(rate*(variableCoefficient*coordinate + disp),globalDof);  //modified by Sunil //xt:rate=1.0 by default;

		}
	}

}
//==========================================================
