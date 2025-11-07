#include "stdafx.h"

#include "BasicAssembler.hpp"
#include "models/BasicModel.hpp"
#include <omp.h>

extern int verboseFlag;

//======================================================================
void BasicAssembler::allocateAndInitialize(BasicModel *model)
{
	BETA_num_threads_for_assembly = model->BETA_num_threads;

	for(int i=0;i<BETA_num_threads_for_assembly;i++){
	
		threadWorkspace *tWorkspace=0;
		tWorkspace = new threadWorkspace;
		tWorkspace->initialize();
		threadWorkspaceList.add(*tWorkspace);
		ElementWorkspace *eWorkspace=0;
		model->CreateAndAllocateWorkspace(eWorkspace); // Maybe we can move this into Assembler
		threadWorkspaceList[i].eWorkspace = eWorkspace;  

	
		FileManager *filemanagerThread=0;
		filemanagerThread = new FileManager;
		filemanagerThread->initialize(filemanager->inputFilename);
		threadWorkspaceList[i].filemanager=filemanagerThread;
		threadWorkspaceList[i].setThreadFileName(i);
		threadWorkspaceList[i].setDefaultOutputStreamForFileManager();
		copyMaterialList(threadWorkspaceList[i].materialList);
	
		}

}

//=================================================================================================
void BasicAssembler::copyMaterialList(Array<Material> &newList)
{
	int numMats=materialList->getNumElements();
	newList.reserve(numMats);
	for(int i=0; i<numMats; i++){
		Material* newMat=(*materialList)[i].clone();
		newList.add(*newMat);
	}	
}