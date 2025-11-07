#include "stdafx.h"

#include <typeinfo>

#include "models/BasicModel.hpp"
#include "Factory.hpp"
#include "utility/utility.h"

void OpenTimeLogFile(FileManager *filemanager);
extern double findNodeTOLERANCE;
FileManager* createFileManager(int type);

bool readScript(Factory *factory, BasicModel* &modelObject, string filename, ifstream *& inStream)
{
//	banner(filename, *filemanager->outStream);
	if( doesFileExist(filename) == false)
		{cout<<"< " << filename <<" >  not found: I quit!\n"; exit(1);}

	inStream= new ifstream(filename.c_str());
	char command[256];
	*inStream >> command;
//	banner(command,  *filemanager->outStream);


	if( COMPARE(command,"createModel")==0 ){
		factory->ModelFactory(inStream,modelObject);
	}else{
		inStream->close(); delete inStream;
		inStream= new ifstream(filename.c_str());
		cout << "No createModel command found! \nCreating BasicModel by default..." <<endl;
		modelObject = new BasicModel;  
	}
	FileManager* filemanager=createFileManager(modelObject->FileManagerType);
	filemanager->initialize(filename);    // Change the name of this method... setDirectoryLocation?
	filemanager->setDefaultOutputStream();
	OpenTimeLogFile(filemanager);

	modelObject->setFileManager(filemanager); //set filemanager
	modelObject->equations.setOuputStream(filemanager->outStream); //set filemanager

	BETA_OUT	<< "Model object created using " << typeid(*modelObject).name() << endl;

	BETA_OUT	<< "DEFAULT findNodeTOLERANCE = " << DOUBLE_FORMAT << findNodeTOLERANCE << endl;
	cout		<< "DEFAULT findNodeTOLERANCE = " << DOUBLE_FORMAT << findNodeTOLERANCE << endl;

	modelObject->readCommands(inStream);

	return true;
}