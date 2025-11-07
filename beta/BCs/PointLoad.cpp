#include "stdafx.h"

#include "PointLoad.hpp"
#include <cmath>

int getNodeIndex(NodeGroup *nlist, double x, double y, double z);

//==========================================================
PointLoad::PointLoad()
{
	setType("PointLoad");
	nodeNum=direction=globalDof=-1;
	force=0;
}

//==========================================================
bool PointLoad::read(istream * inStream)
{
	char  *localTokenList[20];
	int    numberOfTokens;
	coordSpecified=false;

	if( getLineAndTokenize(inStream,"end",localTokenList, numberOfTokens)==1) {
		return false;
	}

	if(numberOfTokens==5){
		point.x=atof(localTokenList[0]);
		point.y=atof(localTokenList[1]);
		point.z=atof(localTokenList[2]);
		force=atof(localTokenList[3]);
		direction=atoi(localTokenList[4]);
		coordSpecified=true;
	}
	else if(numberOfTokens==3){
		nodeNum=atoi(localTokenList[0]);
		if(nodeNum<0) return false;
		force=atof(localTokenList[1]);
		direction=atoi(localTokenList[2]);
	}
	else if(numberOfTokens==1){
		nodeNum=atoi(localTokenList[0]);
		if(nodeNum<0) return false;
	}

	print(BETA_OUT);
	return true;
}
//==========================================================
void PointLoad::print(ostream &ostrm)
{
	if(!coordSpecified)
	ostrm << "Point Load: " << nodeNum << '\t' <<
						force << '\t' << direction	<< '\n';
	else
	ostrm << "Point Load: coords " << point.x << '\t' <<
											point.y << '\t' <<
											point.z << '\n';
	ostrm <<"load: " <<force << "\t in " << "x"<<direction	<<" direction (x1 x2 x3 system)\n";
}

//==========================================================
void PointLoad::applyTo(Equations &equations, NodeGroup *node, int numNodes, BasicModel *model)
{                         
	WhoMethod("applyTo");
	if(coordSpecified==true){
//		for (int i=0; i< numNodes; i++){
//			if ( (point)==(*node)[i] ){
//				nodeNum=(*node)[i].nodeNum;
//			}
//		}
		nodeNum=getNodeIndex(node, point.x, point.y, point.z);
	}
	if(nodeNum==-1) {
		cerr << "Point load error: node not found with given coordinates!" << endl;
		FatalError("Bad Loading.\n");
	}

	if(nodeNum >= numNodes) {
		cerr << "Node " << nodeNum << " not in range of " << numNodes << ".\n";
		FatalError("Bad Loading.\n");
	}
	globalDof = (*node)[nodeNum].getFirstDof()+direction-1;
		 // BETA_OUT<<"Adding force: "<< force <<" to eqn. "
		//		<<( node[nodeNum].getFirstDof() +direction-1 )<<endl;

		//		BETA_OUT<< globalDof <<" "<<  node[nodeNum].getFirstDof()<<" "
		//			 <<  direction <<endl;
    // xtang 990418: rate=1.0 by default !!!
	equations.addLoadToEqn(force,globalDof);

}
//==========================================================
