#include "stdafx.h"

#include "utility/utility.h"
#include "utility/FileManager.hpp"
#include "MacroElement.hpp"

#include "utility/excepts.hpp"
#include "math/matrix.hpp"
#include "math/equation/equation.hpp"
#include "factory/Factory.hpp"

extern int verboseFlag;
extern Factory	*factory;
bool readScript(Factory *factory, BasicModel* &modelObject, string filename, ifstream *& inStream);
//=====================================================================
MacroElement::MacroElement(void)
{
	model=0;
	modelOwner=0;
	MasterKe=0;
}
//=====================================================================
MacroElement::~MacroElement()
{
	if(modelOwner==this){delete model;model=0;}
	if(MasterKe)  {delete MasterKe;MasterKe=0;}
}
//=====================================================================
/*
void MacroElement::readSpecialCommand(istream &inStream, ElementGroup *element, int numElements, char * command)
{
BETA_OUT<<"MacroElement::readSpecialCommand"<<endl;

if (COMPARE(command, "DefineMacroElementGroups")==0 ) {
	//calculate masterKe using first element in the ElementGroup
	MacroElement* master=(MacroElement*)(&(*element)[0]);

	BasicModel *modelObject =0; 
	ifstream *subdomainstream=0;
	string filename;
	inStream >> filename;
	readScript(factory, modelObject, filename, subdomainstream);
	displayTime("Time to finish running Subdomain script", *modelObject->filemanager->outStream);
	if(subdomainstream) delete subdomainstream;

	master->model=modelObject;
	master->modelOwner=master;
	modelObject->allocateForAnalysis(BETA_DO_NOT_ALLOCATE_EQUATIONS);
	master->setPointers(&modelObject->elementWorkspaceList[0]);//using the subdomain model's workspace to calculate the master Ke
	master->calculateDofList();
	master->CalculateK(NULL);

	//linking all the other elements to the subdomain model as well
	for(int i=1;i<numElements;i++) { 
		MacroElement* e=(MacroElement*)(&(*element)[i]);
		e->model=modelObject;
		e->modelOwner=master;
	}
	return;
}
// If no match, try super class
BETA_OUT<<"No match... try super class"<<endl;
ElasticityElement3D::readSpecialCommand(inStream,element,numElements, command);
}//readSpecialCommand
*/
//==========================================================
