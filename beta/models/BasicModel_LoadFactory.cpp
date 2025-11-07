#include "stdafx.h"

#include "BasicModel.hpp"
#include "BCs/BC.hpp"
//=========================================================================
#include "BCs/PointLoad.hpp"
//=========================================================================
void BasicModel::LoadFactory(istream *inStream)
{
	char name[80];
	while (true) {
		(*inStream)>>name;
		BETA_OUT<<"LoadType = "<< name<<endl;
		if(COMPARE(name,"exitReadLoads")==0)goto _EndReadLoadType;

		//Create Loads...
		if ( !createLoads(name, inStream) )  {
			//No match found !
			BETA_OUT<<"You specified type= "<<name<<endl;
			BETA_OUT << "This type is not defined " << endl;
			exit(1);
		}
	}//end of while

_EndReadLoadType:
	return;
}
//======================================================================
bool BasicModel::createLoads(char* name, istream *inStream)
{
	bool exitLoop;
	Load *newLoad=0;

	IS_LOAD(PointLoad)	
	else
		return false;
}
//======================================================================
