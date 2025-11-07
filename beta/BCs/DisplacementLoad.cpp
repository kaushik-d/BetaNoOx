#include "stdafx.h"

#include <string>
#include "DisplacementLoad.hpp"  

//==========================================================
//eventually use nodegroups to impose displacement loads on ALL NODES
//========================================================
DisplacementLoad::DisplacementLoad()
{
	WhoClass("DisplacementLoad");
	setType("Displacement");
	nodeNum=direction=globalDof=-1;
	disp=0;
	allNodes = false;
}
//==========================================================
bool DisplacementLoad::read(istream * inStream)
{
	char *localTokenList[20];
	int  numberOfTokens;

	if( getLineAndTokenize(inStream,"end",localTokenList, numberOfTokens)==1) {
		FatalError("Error in DisplacementLoad::read(istream * inStream)");
	}

	//check if first token = -1 
	if(numberOfTokens >= 1){
		nodeNum = atoi(localTokenList[0]);
		if(nodeNum == -1 ) return false;
	}

	if(numberOfTokens >= 1){
		nodeNum = atoi(localTokenList[0]);
		if(nodeNum == -999 ) allNodes = true;
	}

	//numTokens has to be 3 or larger
	if(numberOfTokens >= 3){
		disp = atof(localTokenList[1]);
		direction = atoi(localTokenList[2]);
		if(direction<1)
			FatalError("Error in DisplacementLoad::read - bad direction number");
	}else{
		FatalError("Error in DisplacementLoad::read(istream * inStream)");
	}

	//if 4th token is "incremental", then set flag
	if(numberOfTokens == 4){
		string s="incremental";
		if (s==localTokenList[3])
			setIncFlag(true);
	}

	BETA_OUT<<"Node : "<<nodeNum<<endl;
	BETA_OUT<<"disp : "<<disp<<endl;
	BETA_OUT<<"direction : "<<direction<<endl;
	print(BETA_OUT);
	return true;
}
//==========================================================
void DisplacementLoad::print(ostream &ostrm)
{
	ostrm << "Displacement Load - " << nodeNum << '\t' << disp << '\t' << 
			direction << endl;
}
//==========================================================
void DisplacementLoad::applyTo(Equations &equations, NodeGroup *node, int numNodes, BasicModel *model)
{ 
	WhoMethod("applyTo");

	if(allNodes==true){
		for(int i=0; i<numNodes; i++){
		globalDof = (*node)[i].getFirstDof()+direction-1;
		equations.setRestraint(globalDof,true);
		equations.setSolution(rate*disp,globalDof);
		}//i
	}else{//	if(allNodes==false)
		if(nodeNum >= numNodes) { 
			cerr << "Node " << nodeNum << " in range " << numNodes << " does not exist.\n";
			BETA_OUT<<"Node does not exist.\n"<<endl; 
		}
		globalDof = (*node)[nodeNum].getFirstDof()+direction-1;
		equations.setRestraint(globalDof,true);
		equations.setSolution(rate*disp,globalDof);  //xt:rate=1.0 by default;
	}
}
//==========================================================
