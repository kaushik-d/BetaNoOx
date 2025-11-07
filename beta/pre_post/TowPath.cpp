#include "stdafx.h"

#include "models/BasicModel.hpp"
#include "TowPath.hpp"
#include <cmath>

#define TOWPATH_TOLERANCE 1e-4
ElementGroup*  createElementGroupUsingMaterialNumber(ElementGroup* eList, int matNum);

TowPath::TowPath()
{
	RotationAxis=0;
	filemanager=0;
}

bool TowPath::readCommands(ifstream *defaultStream)
{
	char  *localTokenList[20];
	int    numberOfTokens;
	int    foundMatch;

	ifstream *currentStream;
    BETA_OUT<<"Current stream in 'TowPath::readCommands' ="<<(*defaultStream)<<endl;
	currentStream = defaultStream;

	while(2==2) {
		foundMatch=0;		
		currentStream = defaultStream;		
		if( getLineAndTokenize(currentStream,"end",localTokenList, numberOfTokens)==1) {
			goto _ExitChoose;
		}
		if (COMPARE(localTokenList[0],"exitModel")==0 || COMPARE(localTokenList[0],"exitUseModel")==0) {
			return true;  // xtang: how about the stream consistency???????
		}

		banner(localTokenList[0], BETA_OUT); 
		foundMatch=processCommand(localTokenList,currentStream,numberOfTokens); // decendent class first

		if(foundMatch ==0) {	
			exit(1);
		}
		
	} //end of while

	_ExitChoose:
	BETA_OUT << "Exiting reader"<<endl; 
	//defaultStream->close();				// FIX THIS !!!! 
	//delete defaultStream; //JV111002

	return true;
}


int TowPath::processCommand(charlist localTokenList, ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;
	
	ifstream istr;

	char token[200];
	int  foundMatch=0;

strcpy(token,localTokenList[0]);
BETA_OUT<<"\nCommand processor = 'TowPath::processCommand'          "<<token<<endl;

	//========
	if( COMPARE(token,"TraversalDirection")==0 ) {
		TraversalDirection.x=atof(localTokenList[1]);
		TraversalDirection.y=atof(localTokenList[2]);
		TraversalDirection.z=atof(localTokenList[3]);
		BETA_OUT << "Traversal Direction: " << TraversalDirection << endl;
		OK
	}

	if( COMPARE(token,"Rotation")==0 ) {
		if(COMPARE(localTokenList[1],"x")==0) RotationAxis = 1;
		if(COMPARE(localTokenList[1],"y")==0) RotationAxis = 2;
		if(COMPARE(localTokenList[1],"z")==0) RotationAxis = 3;
		RotationAngle=atof(localTokenList[2]);
		switch(RotationAxis){
		case 1:	TraversalDirection.y =cos(RotationAngle *M_PI/180);
				TraversalDirection.z =sin(RotationAngle *M_PI/180);
				break;
		case 2:	TraversalDirection.z =cos(RotationAngle *M_PI/180);
				TraversalDirection.x =sin(RotationAngle *M_PI/180);
				break;
		case 3:	TraversalDirection.x =cos(RotationAngle *M_PI/180);
				TraversalDirection.y =sin(RotationAngle *M_PI/180);
				break;	
		}
		c=cos(RotationAngle *M_PI/180);
		s=sin(RotationAngle *M_PI/180);

		BETA_OUT << "Traversal Direction: " << TraversalDirection << endl;
		OK
	}

	if( COMPARE(token,"PhaseDifference")==0 ) {
		PhaseDiff=atof(localTokenList[1]);
		OK
	}

	if( COMPARE(token,"Amplitude")==0 ) {
		Amplitude=atof(localTokenList[1]);
		OK
	}

	if( COMPARE(token,"PeriodicLength")==0 ) {
		PeriodicLength=atof(localTokenList[1]);
		OK
	}

	if( COMPARE(token,"CoordinateShift")==0 ) {
		CoordinateShift=atof(localTokenList[1]);
		OK
	}
	
	if( COMPARE(token,"MaterialGroup")==0 ) {
		MaterialGroup=atoi(localTokenList[1]);
		OK
	}

	if(foundMatch ==0){	
		BETA_OUT << "Command: " << token << endl;
		BETA_OUT << "Did not recognize command: " <<token<<endl;
		exit(1);
	}
	
	return foundMatch;
}


double TowPath::findOffset(Node *c)
{
	double x,y,z,z1,offset;

	if		(TraversalDirection.x==1.0)	{x=c->x;y=c->y;z=c->z;} // fix this !
	else if (TraversalDirection.y==1.0) {x=c->y;y=c->z;z=c->x;}
	else if (TraversalDirection.z==1.0) {x=c->z;y=c->x;z=c->y;}

	z1=Amplitude * sin(2 * M_PI * (x - CoordinateShift) + M_PI/2 );
	offset=z-z1;

	return offset;
}

bool TowPath::satisfyEquation(Node *c, double offset)
{
	double x,y,z,z1;

	if		(TraversalDirection.x==1.0)	{x=c->x;y=c->y;z=c->z;}
	else if (TraversalDirection.y==1.0) {x=c->y;y=c->z;z=c->x;}
	else if (TraversalDirection.z==1.0) {x=c->z;y=c->x;z=c->y;}

	z1=Amplitude * sin(2 * M_PI * (x - CoordinateShift) + M_PI/2 );
	z1 += offset;

	if( fabs(z-z1) <= TOWPATH_TOLERANCE)
		return true;
	else 
		return false;
}

//old function - does not work ! 
//element centroids do not satisfy the tow path equation !!
ElementGroup * TowPath::findElementsUsingCentroid(BasicMesh *mesh, Node *c)
{
	ElementGroup *eList=NULL;
	eList = new ElementGroup;

	ElementGroup temp;

	BasicElement *e,*e2;
	Node *ec,*ec2;

	double x,y,z,offset;

	offset=findOffset(c);

	int i,j;
	for(i=0;i<mesh->numElements;i++){
		e=&mesh->element[i];
		if(e->getMaterialNumber() != MaterialGroup) continue;
		ec=&(mesh->element[i].centroid);
		if(satisfyEquation(ec,offset))
			temp.add(*e);
	}

	//now sort in the tow path direction
	if		(TraversalDirection.x==1.0)	{x=c->x;y=c->y;z=c->z;}
	else if (TraversalDirection.y==1.0) {x=c->y;y=c->z;z=c->x;}
	else if (TraversalDirection.z==1.0) {x=c->z;y=c->x;z=c->y;}

	int numElementsInList=temp.getNumElements();

	for(i=0;i<numElementsInList-1;i++){
		for(j=i+1;j<numElementsInList;j++){
			e=&temp[i];
			ec=&(e->centroid);
			e2=&temp[j];
			ec2=&(e2->centroid);
			if ( (*ec)(0) > (*ec2)(0) ) //hardcoded to sort in x right now
				temp.swap(i,j);

		}
	}
	BETA_OUT << endl;
	for(i=0;i<numElementsInList;i++){
		eList->add(temp[i]);
	}
	

	return eList;

}
/*
ElementGroup * TowPath::findElementsUsingMaterialGroup(BasicMesh *mesh)
{
	ElementGroup *eList=NULL;
	eList = new ElementGroup;

	ElementGroup temp;
	//always remember to clear the array if you don't want its elements to be deleted
	//allocating more space than needed... but doing this will avoid the need for reallocs.
	temp.setCapacity(mesh->numElements);

	BasicElement *e,*e2;
	Node *ec,*ec2;

	int i,j;
	for(i=0;i<mesh->numElements;i++){
		e=&mesh->element[i];
		if(e->getMaterialNumber() != MaterialGroup) continue;
		temp.add(*e);
	}

	//now sort in the tow path direction
	int direction;
	if		(TraversalDirection.x==1.0)	{direction=0;}
	else if (TraversalDirection.y==1.0) {direction=1;}
	else if (TraversalDirection.z==1.0) {direction=2;}

	int numElementsInList=temp.getNumElements();

	for(i=0;i<numElementsInList-1;i++){
		for(j=i+1;j<numElementsInList;j++){
			e=&temp[i];
			ec=&(e->centroid);
			e2=&temp[j];
			ec2=&(e2->centroid);
			if(RotationAxis==0){
				if ( (*ec)(direction) > (*ec2)(direction) ) 
					temp.swap(i,j);
			}
			else{
				Point p1=(*ec).rotate(RotationAxis,RotationAngle);
				Point p2=(*ec2).rotate(RotationAxis,RotationAngle);
				if ( (*ec).rotate(RotationAxis,RotationAngle)(direction) > (*ec2).rotate(RotationAxis,RotationAngle)(direction) ) 
					temp.swap(i,j);
			}
		}
	}
	BETA_OUT << endl;

	eList->setCapacity(numElementsInList);
	for(i=0;i<numElementsInList;i++){
		eList->add(temp[i]);
	}

	temp.clear();

	return eList;
}
*/

ElementGroup * TowPath::findElementsUsingMaterialGroup(BasicMesh *mesh)
{
	ElementGroup *eList=NULL;
	eList = createElementGroupUsingMaterialNumber(&mesh->element, MaterialGroup);
	if(eList == NULL) return NULL;
 
	BasicElement *e,*e2;
	Node *ec,*ec2;
	int i,j;

	//now sort in the tow path direction
	int direction;
	if		(TraversalDirection.x==1.0)	{direction=0;}
	else if (TraversalDirection.y==1.0) {direction=1;}
	else if (TraversalDirection.z==1.0) {direction=2;}

	int numElementsInList=eList->getNumElements();

	for(i=0;i<numElementsInList-1;i++){
		for(j=i+1;j<numElementsInList;j++){
			e=&(*eList)[i];
			ec=&(e->centroid);
			e2=&(*eList)[j];
			ec2=&(e2->centroid);
			if(RotationAxis==0){
				if ( (*ec)(direction) > (*ec2)(direction) ) 
					eList->swap(i,j);
			}
			else{
				Point p1=(*ec).rotate(RotationAxis,RotationAngle);
				Point p2=(*ec2).rotate(RotationAxis,RotationAngle);
				if ( (*ec).rotate(RotationAxis,RotationAngle)(direction) > (*ec2).rotate(RotationAxis,RotationAngle)(direction) ) 
					eList->swap(i,j);
			}
		}
	}
	BETA_OUT << endl;

	return eList;
}

#undef TOWPATH_TOLERANCE
