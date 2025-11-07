#include "stdafx.h"

#include "BasicModel.hpp"
#include "BCs/BC.hpp"
#include "mesh/NodeSet.hpp"
using namespace std;

//====================================================================
#define IF_CREATE_CONSTRAINT(B) BETA_OUT<<"\tReading Constraint...\n";      \
							   constraint->B(inStream, this);	  \
							   return true;                       \
//====================================================================
void BasicModel::ConstraintFactory(istream *inStream)
{
	Constraint *newConstraint=0;
	int i,j;

	char *localTokenList[20];
	int  numberOfTokens;
	Node *n=0;

//=============================
	while (true) {
		if( getLineAndTokenize(inStream,"exitreadConstraints",localTokenList, numberOfTokens)==1) {
			return;
		}
		BETA_OUT<<"Name of constraint group = "<< localTokenList[0] <<endl;

		newConstraint= new Constraint;
		int NumDofsConstrained=0,numNodes=0;
		NumDofsConstrained=numberOfTokens-1;

		BETA_OUT<<"The dof at the nodes are u1 u2 etc\n";
		BETA_OUT<<"The dof constrained at each node in the group are\n";
		for(j=0;j<NumDofsConstrained;j++){
			BETA_OUT<<"u"<<atoi(localTokenList[j+1])<<"  ";
		}//j
		BETA_OUT<<endl;

		const NodeSet* nlist = mesh->getExistingNodeSetPtr(localTokenList[0]);
		if(nlist==0) FatalError("ConstraintFactory : NodeGroup not found!\n");
		numNodes=nlist->getNumNodes();
		//fill up the dofs in the constraint...
		for(i=0;i<numNodes;i++){
			n=(*nlist)[i];
			if(n->numDof <= 0){
				BETA_OUT << "Node has no dofs. cannot apply constraint to this node - node num:" << n->getNodeNum() << endl;
			}else
			for(j=0;j<NumDofsConstrained;j++){
				newConstraint->doflist.push_back(n->getFirstDof() + atoi(localTokenList[j+1]) - 1) ;
				//BETA_OUT << "Node " << n->getNodeNum() << " has been constrained - dof number(starting from 1):" << atoi(localTokenList[j+1]) << endl;
			}
		}
		//Save constraints to list...
		constraints.push_back(newConstraint);
		//Apply constraints...
		ApplyConstraints();
	}//end of while
}
//======================================================================
void BasicModel::ApplyConstraints()
{
	for(int iConstraint=0; iConstraint<(int)constraints.size(); iConstraint++){
		Constraint *c=constraints[iConstraint];
		int numdofs=(int)c->doflist.size();
		for(int i=0; i < numdofs; i++){
			c->SetConstraint(&equations,c->doflist[i]);
		}
	}
}
//======================================================================
