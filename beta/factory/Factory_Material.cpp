#include "stdafx.h"

#include "utility/utility.h"
#include "Factory.hpp"

#include "materials/Material.hpp"
#include "materials/1D/RodMaterial.hpp"
#include "materials/3D/elastic.hpp"
#include "materials/2D/pStress.hpp"
#include "materials/2D/pStrain.hpp"

//======================================================================
bool Factory::MaterialFactory(istream * currentStream, Array<Material> *materialList, FileManager* filemanager)
{
	Material *newMaterial=0;
	char  label[80], command[80];
	int group;

	char *localTokenList[20];
	int  numberOfTokens;

top:
	if( getLineAndTokenize(currentStream,"exitReadMaterials",localTokenList, numberOfTokens)==1) {
		BETA_OUT << "Finished reading materials" << endl;
		return true;
	}
	strcpy(command, localTokenList[0]);
	banner(command, BETA_OUT);
    if(COMPARE(command,"exitReadMaterials")==0) {return true;}	
	(*currentStream) >>group;
			  			  
	(*currentStream).getline( label, 80,'\n'); 
	BETA_OUT <<"Group "<<group<<" "<<label<<endl;

	if(!createMaterials(newMaterial, command)){
		BETA_OUT <<"Material Type: "<<command<<endl;;
		banner("Unrecognized material type", BETA_OUT);
		exit(1);
	}

	newMaterial->setFileManager(filemanager);
	newMaterial->setGroupName(label);
	newMaterial->setGroupNum(group);
	newMaterial->initialize();

	newMaterial->read(currentStream);
	materialList->putAt(group,newMaterial);
	goto top;
}  //...end of function
//======================================================================
bool Factory::createMaterials(Material * &newMaterial, char* command)
{
    Create_material_ifTokenIs(Material)		else
    Create_material_ifTokenIs(RodMaterial)	else
    Create_material_ifTokenIs(ElasticMaterial)	else
    Create_material_ifTokenIs(PlaneStress)		else
    Create_material_ifTokenIs(PlaneStrain)		else
		return false;

}
//======================================================================


