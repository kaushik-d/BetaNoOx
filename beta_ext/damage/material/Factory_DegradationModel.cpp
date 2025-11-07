#include "stdafx.h"

#include "utility/utility_defines.h"

#include "../material/BasicDegradationModel.hpp"
#if !defined NoOxDamage
#include "../material/OxDegradationModel.hpp"
#endif
//========================================================================
BasicDegradationModel* CreateDegradationModel(string DegradationModelName, FileManager *filemanager)
{
	if(COMPARE(DegradationModelName.c_str(),"BasicDegradationModel")==0){
		return new BasicDegradationModel();
    #if !defined NoOxDamage
	}else if(COMPARE(DegradationModelName.c_str(),"OxDegradationModel")==0){
		return new OxDegradationModel();
    #endif
	}else{
		EXIT_BETA("Don't know how to create specified DegradationModel");
	}
}
//========================================================================
