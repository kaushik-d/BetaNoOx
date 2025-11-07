#include "stdafx.h"

#include "ESA_Assembler.hpp"
#include "models/BasicModel.hpp"
#include <omp.h>

extern int verboseFlag;

void ESA_Assembler::assembleK(double * theSolution)
{
// Assemble coefficient matrix
	equations->zeroMatrix();  //This zeroes K and pseudoLoadVector
	
	IntegrationDataType   idt;
	int numNodesPerElement = mesh->element[0].numNodesPerElement;
	idt.id=new InterpAndSolDerivatives(numNodesPerElement,numNodesPerElement); idt.id->allocate();
	idt.id->numberOfInterp = numNodesPerElement;
	mesh->element[0].getQuadraturePoints(idt.gaussPointList); 

	omp_set_num_threads(BETA_num_threads_for_assembly);   //moved BCO 07/24/09

	#pragma omp parallel                   //added bco 07/24/09
	{
	#pragma omp for		//added bco 07/24/09
	for(int i=0; i<BETA_num_threads_for_assembly; i++)
	{
		threadWorkspace *thread=&(threadWorkspaceList[i]);
		BasicElement *e=0;	
		for(int j=thread->minElementCalculationNumber;j<=thread->maxElementCalculationNumber;j++) {
			e=&mesh->element[j];
			e->update("K", theSolution);

		}
	}
	}//parallel

	//Do Assembly In Parallel

	omp_set_num_threads(BETA_num_threads_for_assembly);
	#pragma omp parallel
	{
		#pragma omp for
		for(int i=0; i<BETA_num_threads_for_assembly; i++)
		{
		threadWorkspace *thread = &(threadWorkspaceList[i]);
		int numElAssignedtoThread = thread->elementList.getNum();
		BasicElement *te=0;

			for(int j=0;j<numElAssignedtoThread;j++)
			{
				te = &mesh->element[thread->elementList[j]];
				
				
				Matrix *Ke = te->getStiffnessMatrix();
				int *dofList = te->getDofList();
				int numDof = te->getNumDof();
				equations->addMatrix_Parallel(*Ke,
				                 dofList,
				                 numDof,
				                 NULL,thread);

								 
			}
		
		}
	}

}
//========================================================================
void ESA_Assembler::assembleKandF(double * theSolution)
{
// Assemble coefficient matrix
	equations->zeroMatrix();  //This zeroes K and pseudoLoadVector
	
	IntegrationDataType   idt;
	int numNodesPerElement = mesh->element[0].numNodesPerElement;
	idt.id=new InterpAndSolDerivatives(numNodesPerElement,numNodesPerElement); idt.id->allocate();
	idt.id->numberOfInterp = numNodesPerElement;
	mesh->element[0].getQuadraturePoints(idt.gaussPointList); 

	omp_set_num_threads(BETA_num_threads_for_assembly);   //moved BCO 07/24/09

	#pragma omp parallel                   //added bco 07/24/09
	{
	#pragma omp for		//added bco 07/24/09
	for(int i=0; i<BETA_num_threads_for_assembly; i++)
	{
		threadWorkspace *thread=&(threadWorkspaceList[i]);
		BasicElement *e=0;	
		for(int j=thread->minElementCalculationNumber;j<=thread->maxElementCalculationNumber;j++) {
			e=&mesh->element[j];
			e->update("KI", theSolution);

		}
	}
	}//parallel


	//Do Assembly In Parallel

	omp_set_num_threads(BETA_num_threads_for_assembly);
	#pragma omp parallel
	{
		#pragma omp for
		for(int i=0; i<BETA_num_threads_for_assembly; i++)
		{
		threadWorkspace *thread = &(threadWorkspaceList[i]);
		int numElAssignedtoThread = thread->elementList.getNum();
		BasicElement *te=0;

			for(int j=0;j<numElAssignedtoThread;j++)
			{
				te = &mesh->element[thread->elementList[j]];
				
				
				Matrix *Ke = te->getStiffnessMatrix();
				double *Fe = te->getInitialForceVector();
				int *dofList = te->getDofList();
				int numDof = te->getNumDof();
				equations->addMatrix_Parallel(*Ke,
				                 dofList,
								 numDof,
				                 Fe,thread);
								 
			}
		
		}
	} //parallel

}
//========================================================================

void ESA_Assembler::assembleInitialF()
{
	// This is no good for incremental analysis !!!! ..... fix !!


// Current displacement state is required to properly calculate the
// nonlinear B matrix.... also, what about strain calculation?

// Assemble coefficient matrix

	equations->zeroMatrix();  //This zeroes K and pseudoLoadVector
	
	IntegrationDataType   idt;
	int numNodesPerElement = mesh->element[0].numNodesPerElement;
	idt.id=new InterpAndSolDerivatives(numNodesPerElement,numNodesPerElement); idt.id->allocate();
	idt.id->numberOfInterp = numNodesPerElement;
	mesh->element[0].getQuadraturePoints(idt.gaussPointList); 

	omp_set_num_threads(BETA_num_threads_for_assembly);   //moved BCO 07/24/09
	#pragma omp parallel               //added bco 07/24/09
	{
	#pragma omp for			
		for(int i=0; i<BETA_num_threads_for_assembly; i++)
		{
			BasicElement *e=0;
			threadWorkspace *thread = &(threadWorkspaceList[i]);
			for(int j=thread->minElementCalculationNumber;j<=thread->maxElementCalculationNumber;j++)
			{
			e=&mesh->element[j];
			e->update("I", NULL);

			}
		}
	}

	//Do Assembly In Parallel
	omp_set_num_threads(BETA_num_threads_for_assembly);
	#pragma omp parallel
	{
		#pragma omp for
		for(int i=0; i<BETA_num_threads_for_assembly; i++)
		{
		threadWorkspace *thread = &(threadWorkspaceList[i]);
		int numElAssignedtoThread = thread->elementList.getNum();
		BasicElement *te=0;

			for(int j=0;j<numElAssignedtoThread;j++)
			{
				te = &mesh->element[thread->elementList[j]];
				
				double *Fe = te->getInitialForceVector();
				int *dofList = te->getDofList();
				int numDof = te->getNumDof();

				equations->addToLoadVector_Parallel(dofList,
					              numDof,
					              Fe,thread);
								 
			}
		
		}
	}

}
//========================================================================
void ESA_Assembler::assignElementWorkspacePointers(ElementWorkspace * elemWorkspace, int minElement, int maxElement)
{
	//This function was added to allow certian ranges of elements to point to a common ElementWorkspace.
	//This was done so that this range of elements would not compete with another range of elements
	//when run in parallel. phh 7/17/08

	//set element pointers
	BasicElement *e=0;
	
	BETA_OUT<<"Workspace: "<<elemWorkspace<<endl;

	for(int i=minElement; i<=maxElement; i++) //Check element static array sizes
	{
		e=&mesh->element[i];
		e->setPointers(elemWorkspace);
	}
}
//========================================================================
void ESA_Assembler::assignElementFilemanagerPointers(FileManager *filemanagerptr, int minElement,int maxElement)
{
	//This function was added to allow certian ranges of elements to point to a common filemanger.
	//This was done so that this range of elements would not compete with another range of elements
	//when run in parallel. phh 7/17/08

	//set element pointers
	BasicElement *e=0;

	for(int i=minElement; i<=maxElement; i++) //Check element static array sizes
	{
		e=&mesh->element[i];
		e->setFileManager(filemanagerptr);
	}
}
//========================================================================
void ESA_Assembler::assignElementMaterialPointers(Array<Material> &mList,int minElement, int maxElement)
{
	BasicElement *e = 0;
	for(int i = minElement; i <= maxElement; i++)
	{
		e = &mesh->element[i];
		int matNum=e->getMaterialGroup();
		Material *matPoint = &(mList[matNum]);
		e->setMaterial(matPoint);
	}
}

//========================================================================
void ESA_Assembler::assignElementsToThreadsForCalculation()
{
	int elementsPerThread = (mesh->numElements)/BETA_num_threads_for_assembly;

	for(int i=0;i<BETA_num_threads_for_assembly;i++)
	{
		threadWorkspaceList[i].minElementCalculationNumber = i*elementsPerThread;
		if( (i+1)==BETA_num_threads_for_assembly)
			threadWorkspaceList[i].maxElementCalculationNumber = (mesh->numElements) - 1;
		else
			threadWorkspaceList[i].maxElementCalculationNumber = (i+1)*(elementsPerThread) - 1;
	}

}
//========================================================================
void ESA_Assembler::allocateAndInitialize(BasicModel *model)
{
	ParallelAssembler::allocateAndInitialize(model);

	int WorkspaceMaxNumberDof = threadWorkspaceList[0].eWorkspace->WorkspaceMaxNumberDof;

	for(int i=0;i<mesh->numElements;i++)
	{
		mesh->element[i].requiresParallelESAConsideration = true;
		mesh->element[i].allocateElementalDataMemory(WorkspaceMaxNumberDof);
	}

	
}

//========================================================================
void ESA_Assembler::makeAssignments(int totalNumDof)
{
	splitEquationsAmongThreads(totalNumDof);
	assignElementsToThreadsForAssembly();

	assignElementsToThreadsForCalculation();

	for(int i=0;i<BETA_num_threads_for_assembly;i++){
		threadWorkspaceList[i].printThreadAssignments();
		
		threadWorkspaceList[i].printThreadCalculationAssignments();
		assignElementWorkspacePointers(threadWorkspaceList[i].eWorkspace,
									   threadWorkspaceList[i].minElementCalculationNumber,
									   threadWorkspaceList[i].maxElementCalculationNumber);
		assignElementFilemanagerPointers(threadWorkspaceList[i].filemanager,
										 threadWorkspaceList[i].minElementCalculationNumber,
										 threadWorkspaceList[i].maxElementCalculationNumber);
    	assignElementMaterialPointers(threadWorkspaceList[i].materialList,
										 threadWorkspaceList[i].minElementCalculationNumber,
										 threadWorkspaceList[i].maxElementCalculationNumber);
		}
}
//========================================================================
void ESA_Assembler::postAssemblyOperations(BasicModel *model)
{
//This functions resets the first threadWorkspace to be that of a 
//typical serial analysis. This is required to do serial post-processing.

	for (int i=0;i<mesh->numElements;i++){
		mesh->element[i].releaseElementalDataMemory();
		mesh->element[i].requiresParallelESAConsideration=false;
	}

	model->parallelAssembly_ESA = false;
	model->CreateAndAllocateWorkspace(threadWorkspaceList[0].eWorkspace);
	model->assignElementWorkspacePointers(threadWorkspaceList[0].eWorkspace);
}

//========================================================================