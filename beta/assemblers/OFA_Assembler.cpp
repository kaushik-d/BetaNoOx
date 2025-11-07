#include "stdafx.h"

#include "OFA_Assembler.hpp"
#include <omp.h>

extern int verboseFlag;
//========================================================================
void OFA_Assembler::assembleInitialF()
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
	
	omp_set_num_threads(BETA_num_threads_for_assembly);
	#pragma omp parallel               
	{
	#pragma omp for			
		for(int i=0; i<BETA_num_threads_for_assembly; i++)
		{
			//Every thread executes the following code concurrently ********************************************
			bool iSetIsBusyTrue = false;
			bool threadFinished = false;

			BasicElement *e=0;
			threadWorkspace *thread = &(threadWorkspaceList[i]);

			while(threadFinished == false)
			{
			for(int j=0;j<thread->secondaryElementList.getNum();j++)
			{
				if(thread->secondaryElementList[j]!=-1)
				{
				iSetIsBusyTrue = false;
				e=&mesh->element[thread->secondaryElementList[j]];
			
				#pragma omp critical (lock1)
				{
				//only one thread can enter this section at any time
					if (e->isBusy == false)
					{
					   e->isBusy = true;
					   iSetIsBusyTrue = true;
					}
				}// end critical section 

				if (iSetIsBusyTrue == true)
				{

					//set element's pointers to the thread's workspace, filemanager, and material list
					assignElementWorkspacePointers(thread->eWorkspace, e->elementNumber);
					assignElementFilemanagerPointers(thread->filemanager, e->elementNumber);
					assignElementMaterialPointers(thread->materialList, e->elementNumber);

					//calculate Fe and Assemble
					e->calculateDofList();
					e->update("I", NULL);
					double *Fe = e->getInitialForceVector();

					equations->addToLoadVector_Parallel(e->calculateDofList(),
												 e->getNumDof(),
												 Fe,thread);
					

					thread->secondaryElementList[j] = -1;
					
					//release element so other processors may access it
					e->isBusy = false;
				}
				}
			//***************************************************************************************************				
			} //end of element loop
			threadFinished = thread->isThreadFinishedF();
			} //end of while loop (is thread finished check)
		} //end of omp for
	}// end omp parallel

}
//========================================================================

void OFA_Assembler::assembleK(double *theSolution)
{
	// Assemble coefficient matrix
	equations->zeroMatrix();  //This zeroes K and pseudoLoadVector
	
	//threadWorkspace *thread=0;
	

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
			//Every thread executes the following code concurrently ********************************************
			bool iSetIsBusyTrue = false;
			bool threadFinished = false;

			BasicElement *e=0;
			threadWorkspace *thread = &(threadWorkspaceList[i]);

			while(threadFinished == false)
			{
			for(int j=0;j<thread->elementList.getNum();j++)
			{
				
				
				if(thread->elementList[j]!=-1)
				{
				iSetIsBusyTrue = false;
				e=&mesh->element[thread->elementList[j]];
			
				#pragma omp critical (lock2)
				{
				//only one thread can enter this section at any time
					if (e->isBusy == false)
					{
					   e->isBusy = true;
					   iSetIsBusyTrue = true;
					}
				}// end critical section 

				if (iSetIsBusyTrue == true)
				{

					//set element's pointers to the thread's workspace, filemanager, and material list
					assignElementWorkspacePointers(thread->eWorkspace, e->elementNumber);
					assignElementFilemanagerPointers(thread->filemanager, e->elementNumber);
					assignElementMaterialPointers(thread->materialList, e->elementNumber);

					int *doflist=e->calculateDofList();
					e->update("K", theSolution);
					Matrix *Ke = e->getStiffnessMatrix();


					equations->addMatrix_Parallel((*Ke), 
				                 e->calculateDofList(),
				                 e->getNumDof(),NULL,thread); 
					

					thread->elementList[j] = -1;
					
					//release element so other processors may access it
					e->isBusy = false;
				}
				}
			//***************************************************************************************************				
			} //end of element loop
			threadFinished = thread->isThreadFinishedK();
			} //end of while loop (is thread finished check)
		} //end of omp for
	}// end omp parallel
}
//========================================================================

//========================================================================

void OFA_Assembler::assembleKandF(double *theSolution)
{
	// Assemble coefficient matrix
	equations->zeroMatrix();  //This zeroes K and pseudoLoadVector
	
	//threadWorkspace *thread=0;
	

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
			//Every thread executes the following code concurrently ********************************************
			bool iSetIsBusyTrue = false;
			bool threadFinished = false;

			BasicElement *e=0;
			threadWorkspace *thread = &(threadWorkspaceList[i]);

			while(threadFinished == false)
			{
			for(int j=0;j<thread->elementList.getNum();j++)
			{
				
				
				if(thread->elementList[j]!=-1)
				{
				iSetIsBusyTrue = false;
				e=&mesh->element[thread->elementList[j]];
			
				#pragma omp critical (lock3)
				{
				//only one thread can enter this section at any time
					if (e->isBusy == false)
					{
					   e->isBusy = true;
					   iSetIsBusyTrue = true;
					}
				}// end critical section 

				if (iSetIsBusyTrue == true)
				{

					//set element's pointers to the thread's workspace, filemanager, and material list
					assignElementWorkspacePointers(thread->eWorkspace, e->elementNumber);
					assignElementFilemanagerPointers(thread->filemanager, e->elementNumber);
					assignElementMaterialPointers(thread->materialList, e->elementNumber);

					int *doflist=e->calculateDofList();
					e->update("KI", theSolution);
					Matrix *Ke = e->getStiffnessMatrix();
					double *Fe = e->getInitialForceVector();


					equations->addMatrix_Parallel((*Ke), 
				                 e->calculateDofList(),
				                 e->getNumDof(),Fe,thread); 
					

					thread->elementList[j] = -1;
					
					//release element so other processors may access it
					e->isBusy = false;
				}
				}
			//***************************************************************************************************				
			} //end of element loop
			threadFinished = thread->isThreadFinishedK();
			} //end of while loop (is thread finished check)
		} //end of omp for
	}// end omp parallel

	resetThreadElementLists();
}
//========================================================================

void OFA_Assembler::assignElementWorkspacePointers(ElementWorkspace * elemWorkspace,int elNum)
{
	//This function was added to allow certian ranges of elements to point to a common ElementWorkspace.
	//This was done so that this range of elements would not compete with another range of elements
	//when run in parallel. phh 7/17/08

	//set element pointers
	BasicElement *e=0;
	e=&mesh->element[elNum];
	e->setPointers(elemWorkspace);

}

//========================================================================
void OFA_Assembler::assignElementFilemanagerPointers(FileManager *filemanagerptr, int elNum)
{
	//This function was added to allow certian ranges of elements to point to a common filemanger.
	//This was done so that this range of elements would not compete with another range of elements
	//when run in parallel. phh 7/17/08

	//set element pointers
	BasicElement *e=0;
	e=&mesh->element[elNum];
	e->setFileManager(filemanagerptr);

}
//========================================================================
void OFA_Assembler::assignElementMaterialPointers(Array<Material> &mList, int elNum)
{
	BasicElement *e = 0;
	e = &mesh->element[elNum];
	int matNum=e->getMaterialGroup();
	Material *matPoint = &(mList[matNum]);
	e->setMaterial(matPoint);

}


//========================================================================

//======================================================================
void OFA_Assembler::resetThreadElementLists()
{
	for (int i=0;i<BETA_num_threads_for_assembly;i++) {
		threadWorkspace *thread = &(threadWorkspaceList[i]);
		
		for(int j=0;j<thread->elementList.getNum();j++) {

			thread->elementList[j] = thread->secondaryElementList[j];

		}
	}
}
//========================================================================