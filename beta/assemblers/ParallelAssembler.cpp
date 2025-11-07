#include "stdafx.h"

#include "ParallelAssembler.hpp"
#include "models/BasicModel.hpp"

extern int verboseFlag;

//========================================================================
void ParallelAssembler::splitEquationsAmongThreads(int totalNumDof)
{
	int equationsPerThread = totalNumDof/BETA_num_threads_for_assembly;

	for(int i=0;i<BETA_num_threads_for_assembly;i++)
	{
		threadWorkspaceList[i].minEquationNumber = i*equationsPerThread;
		if( (i+1)==BETA_num_threads_for_assembly)
			threadWorkspaceList[i].maxEquationNumber = totalNumDof - 1;
		else
			threadWorkspaceList[i].maxEquationNumber = (i+1)*(equationsPerThread) - 1;
	}

}
//=======================================================================
void ParallelAssembler::assignElementsToThreadsForAssembly()
{
	BasicElement *e=0;
	for(int i=0; i<mesh->numElements; i++) 
	{
		e=&mesh->element[i];
		e->AssignElementsToThreads(threadWorkspaceList);
	}


}
//=======================================================================
void ParallelAssembler::makeAssignments(int totalNumDof)
{
	splitEquationsAmongThreads(totalNumDof);
	assignElementsToThreadsForAssembly();
	
	for(int i=0;i<BETA_num_threads_for_assembly;i++){
		threadWorkspaceList[i].printThreadAssignments();
	}

}
//========================================================================
