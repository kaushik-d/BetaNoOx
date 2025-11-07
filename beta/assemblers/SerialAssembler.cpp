#include "stdafx.h"

#include "SerialAssembler.hpp"
#include "models/BasicModel.hpp"
#include <omp.h>

extern int verboseFlag;

void SerialAssembler::assembleInitialF()
{
// This is no good for incremental analysis !!!! ..... fix !!


// Current displacement state is required to properly calculate the
// nonlinear B matrix.... also, what about strain calculation?


	equations->zeroMatrix();

	BasicElement *e=0;
	for(int i=0;i<mesh->numElements;i++) {
		e=&mesh->element[i];
		e->calculateDofList();
		//element[i]->update("I", NULL);
		e->update("I", zeroVector);        // changed 6/98 ... make sure zeroVector has been set to zero!!!!
        double *Fe = e->getInitialForceVector();
		equations->addToLoadVector(e->calculateDofList(),
				              e->getNumDof(),
				              Fe);
	}	
}
//========================================================================

//========================================================================
void SerialAssembler::assembleK(double *theSolution)
{

	// Assemble coefficient matrix
	equations->zeroMatrix();  //This zeroes K and pseudoLoadVector
	BasicElement *e=0;

	for(int i=0;i<mesh->numElements;i++) {
		e=&mesh->element[i];
		int *doflist=e->calculateDofList();
		e->update("K", theSolution);
		Matrix *Ke = e->getStiffnessMatrix();

		equations->addMatrix((*Ke), 
				                 e->calculateDofList(),
				                 e->getNumDof(),
				                 NULL);

		if(verboseFlag == Max){
			BETA_OUT <<"Element #" << e->elementNumber << endl;
//			Ke->print("K_e", e->getNumDof(),e->getNumDof(),&BETA_OUT);
			Ke->printWithRowColumnNumbering("K_e", e->getNumDof(),e->getNumDof(), doflist, doflist, &BETA_OUT);
		}
	}
}
//========================================================================

//========================================================================
// Only use this routine for servicing "iterate"
void SerialAssembler::assembleKandF(double *theSolution)   // Move to nonlinear model class
{
// This assembles the forces Fe into the pseudoLoadVector ... not the real loadVector ====>>
//  do not use for distributed loads

	equations->zeroMatrix();
	BasicElement *e=0;

	for(int i=0;i<mesh->numElements;i++) {
		e=&mesh->element[i];
		e->update("KI", theSolution);
		Matrix *Ke = e->getStiffnessMatrix();
		double *Fe = e->getInitialForceVector();

		equations->addMatrix((*Ke), 
				                 e->calculateDofList(),
				                 e->getNumDof(),
				                 Fe);
		if(verboseFlag == Max){
			int *doflist=e->calculateDofList();
			BETA_OUT <<"Element #" << e->elementNumber << endl;
//			Ke->print("K_e", e->getNumDof(),e->getNumDof(),&BETA_OUT);
			Ke->printWithRowColumnNumbering("K_e", e->getNumDof(),e->getNumDof(), doflist, doflist, &BETA_OUT);
		}
	}
}
//========================================================================

void SerialAssembler::makeAssignments(int totalNumDof)
{
	assignElementWorkspacePointers(threadWorkspaceList[0].eWorkspace);
	assignElementFilemanagerPointers(threadWorkspaceList[0].filemanager);
	assignElementMaterialPointers(threadWorkspaceList[0].materialList);
}
//========================================================================

void SerialAssembler::assignElementWorkspacePointers(ElementWorkspace * elemWorkspace)
{
	int numEl = mesh->numElements;
	for(int i=0; i<numEl; i++){
		BasicElement *e=0;
		e=&mesh->element[i];
		e->setPointers(elemWorkspace);
	}

}

//========================================================================
void SerialAssembler::assignElementFilemanagerPointers(FileManager *filemanagerptr)
{
	int numEl = mesh->numElements;
	for(int i=0; i<numEl; i++){
		BasicElement *e=0;
		e=&mesh->element[i];
		e->setFileManager(filemanagerptr);
	}

}
//========================================================================
void SerialAssembler::assignElementMaterialPointers(Array<Material> &mList)
{
	int numEl = mesh->numElements;
	for(int i=0; i<numEl; i++){
		BasicElement *e=0;
		e=&mesh->element[i];
		int matNum=e->getMaterialGroup();
		Material *matPoint = &(mList[matNum]);
		e->setMaterial(matPoint);
	}

}


//========================================================================
