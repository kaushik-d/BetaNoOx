#include "stdafx.h"

#include "Factory.hpp"
#include "models/ElasticityModel.hpp"
#include "pre_post/TowTraversalUtility.hpp"

/*========================================================

To remove a element from this project, do the following
	* delete the include file
	* delete the corresponding "IF_CREATE( ....." statement below

========================================================*/
//====================================================================
void Factory::ModelFactory(istream *inStream, BasicModel *&model)
{
	char *localTokenList[20];
	int  numberOfTokens;

//=============================
	while (true) {
		if( getLineAndTokenize(inStream,"exitcreateModel",localTokenList, numberOfTokens)==1) {
			goto _EndSetModelType;
		}
		cout <<"ModelType = "<< localTokenList[0] <<endl;
		//Create Models...
		if (!createModels(model, localTokenList[0]))  {
			//No match found !
			cout <<"You specified type= "<< localTokenList[0] <<endl;
			cout << "This factory cannot create this model object" << endl;
			exit(1);
		}
		model->initialize();
		if(numberOfTokens == 2)
			model->setFileManagerType(localTokenList[1]);
	}//end of while

_EndSetModelType:
	return;
}  // end of CreateElements
//======================================================================
bool Factory::createModels(BasicModel *&model, char* name)
{
	IF_CREATE_MODEL("BasicModel",       BasicModel)		else
    IF_CREATE_MODEL("ElasticityModel",		ElasticityModel)		else
    IF_CREATE_MODEL("TowTraversalUtility",	TowTraversalUtility)	else
	
		return false;

}
//======================================================================
