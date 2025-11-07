#include "stdafx.h"

#include <string>
#include "PlaneDisplacement.hpp"  

NodeGroup*  getNodesOnPlane(NodeGroup *nlist, int coordnum, double coordval);
//==========================================================
PlaneDisplacement::PlaneDisplacement()
{
	WhoClass("PlaneDisplacement");
	setType("Displacement");
	planeCoordNum=direction=0;
	coordValue=displacement=0;
}
//==========================================================
bool PlaneDisplacement::read(istream * inStream)
{
   (* inStream) >> planeCoordNum;   //  1->3
	if( planeCoordNum < 0 ) return false;
	
	if(planeCoordNum < 1 || planeCoordNum > 3)        
		FatalError("Bad plane coordinate number...");
	(* inStream) >> coordValue ;
	(* inStream) >> displacement ;
	(* inStream) >> direction ;
	if(direction < 1)
		FatalError("Bad direction number...");
	print(BETA_OUT);
	return true;
}
//==========================================================
void PlaneDisplacement::print(ostream &ostrm)
{
	ostrm <<
		"\tPlane          = " << planeCoordNum <<
		"\n\tCoordinate   = " << coordValue    <<
		"\n\tdisplacement = " << displacement  <<
		"\n\tdirection    = " << direction     << '\n';
}
//==========================================================
void PlaneDisplacement::applyTo(Equations &equations, NodeGroup *node, int numNodes, BasicModel *model)
{ 
	int i,globalDof;
	WhoMethod("applyTo");

	NodeGroup *nodelist=0;
	nodelist = getNodesOnPlane(node, planeCoordNum, coordValue);
	if(nodelist == 0) {
		cout << "PlaneDisplacement::applyTo() - No nodes found in specified plane!" << endl;
		exit(1);
	}

	int nn;
	nn=nodelist->getNumElements();

	for(i=0;i<nn;i++) {
		if( (*nodelist)[i].getNumDof() < 1) {
			BETA_OUT << "Node has no dofs. cannot apply load to this node - node num:" << (*nodelist)[i].getNodeNum() << endl;
			continue;
		}
		globalDof = (*nodelist)[i].getFirstDof()+direction-1;
		equations.setRestraint(globalDof,true);
		//xt:rate=1.0 by default;
		equations.setSolution(displacement*rate,globalDof);
		BETA_OUT << "Node " << (*nodelist)[i].getNodeNum() << " has been loaded - dof number(starting from 1):" << direction << endl;

	};  

	//release memory for nodelist
	nodelist->clear();
	delete nodelist;
}
//==========================================================
