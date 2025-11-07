#include "stdafx.h"

#include <cstdlib>  
#include "mesh/BasicMesh.hpp"

void numberOfElementsAttachedToNodes(int numNodes, 
                             int numElem, ElementGroup &element, bool averageOnlyActive,
                            int **&numberAttachedInGroup, int MaxNumberOfMaterialsForContouring);

//==========================================================================

void averageDataValues(int numNodes, int numElem, ElementGroup &element, 
	bool averageOnlyActive, int MaxNumberOfMaterialsForContouring)
{
	/****************************************************************
		Determine number of elements connected to each node &
		Average Within a Group
	*****************************************************************/

int *numberAttached;
int currentGroup,i,j,elem,nnpe,n,materialGroup;

double **globalNodalValues=0;
globalNodalValues = (double **)malloc((MaxNumberOfMaterialsForContouring) * sizeof(double*));
if(globalNodalValues==0) EXIT_BETA("Memory Allocation Error!");
for(i=0; i<MaxNumberOfMaterialsForContouring; i++){ 
	globalNodalValues[i] = (double *)malloc((numNodes+1) * sizeof(double));
	if(globalNodalValues[i]==0) EXIT_BETA("Memory Allocation Error!");
	for(j=0; j<numNodes;j++) globalNodalValues[i][j] = 0.;
}

int **numberAttachedInGroup;
numberAttachedInGroup = (int **)malloc((MaxNumberOfMaterialsForContouring) * sizeof(int*));
if(numberAttachedInGroup==0) EXIT_BETA("Memory Allocation Error!");
for(i=0; i<MaxNumberOfMaterialsForContouring; i++)
{ 
	numberAttachedInGroup[i] = (int *)malloc((numNodes+1) * sizeof(int));
	if(numberAttachedInGroup[i]==0) EXIT_BETA("Memory Allocation Error!");
}

if(MaxNumberOfMaterialsForContouring <1 ) 
	EXIT_BETA("MaxNumberOfMaterialsForContouring is less that 1. Not an acceptable value!");
numberOfElementsAttachedToNodes(numNodes, numElem, element,averageOnlyActive,
                            numberAttachedInGroup, MaxNumberOfMaterialsForContouring);

//************************************
// Sum the nodal values within groups
//************************************
BasicElement *elemObj;		
		for(elem =0; elem< numElem; elem++) {
			elemObj = &element[elem];
			if(!averageOnlyActive || elemObj->activeElementFlag) {
					materialGroup = elemObj->getMaterialNumber();
					nnpe = elemObj->numNodesPerElement;
					for(i = 0; i<nnpe;i++) {
						//n = elemObj->connectivity[i];
						n = elemObj->node[i].nodeNum;
						globalNodalValues[materialGroup][n] += 
                                         elemObj->nodalValues[i];
					}//i
			}// if
		}// elem

/*
Element *elemObj;		
		for(elem =0; elem< numElem; elem++) {
			elemObj = element[elem]
			if(!averageOnlyActive || element[elem].activeElementFlag) {
					materialGroup = element[elem].materialGroup;
					nnpe = element[elem].numNodesPerElement;
					for(i = 0; i<nnpe;i++) {
						n = element[elem].connectivity[i];
						globalNodalValues[materialGroup][n] += 
                                         element[elem].nodalValues[i];
					}//i
			}// if
		}// elem


*/


//******************************
// Divide by number of attachments
//********************************
double *gnv ;

for(currentGroup = 0; currentGroup< MaxNumberOfMaterialsForContouring; currentGroup++) {
	numberAttached = numberAttachedInGroup[currentGroup]; 
   gnv = globalNodalValues[currentGroup];

		for(i=0; i<numNodes; i++)
			if(numberAttached[i]>0)
				gnv[i] = gnv[i] /  numberAttached[i];
		
//*********************************************************
//Put these averaged values back in the element structure
//*********************************************************
		
		for( elem =0; elem< numElem; elem++) {
			if(!averageOnlyActive || element[elem].activeElementFlag) {
				materialGroup = element[elem].getMaterialNumber();
				if(materialGroup==currentGroup) {
					nnpe = element[elem].numNodesPerElement;
					for(i = 0; i<nnpe;i++) {
						//n = element[elem]->connectivity[i] ;
						n = element[elem].node[i].nodeNum;
						element[elem].nodalValues[i] = gnv[n];
					}// i
				}// if
			}// if
		}//elem
}

	for(i=0; i<MaxNumberOfMaterialsForContouring; i++){ 
		free(numberAttachedInGroup[i]);
		free(globalNodalValues[i]);
	}
	free(numberAttachedInGroup);
	free(globalNodalValues);
}

//================================================================
void numberOfElementsAttachedToNodes(int numNodes, 
                             int numElem, ElementGroup &element, bool averageOnlyActive,
                            int **&numberAttachedInGroup, int MaxNumberOfMaterialsForContouring)
{
int i,j;
int elem,nnpe,n,materialGroup;

for(i=0;i<MaxNumberOfMaterialsForContouring;i++){
for(j=0; j< numNodes; j++) {numberAttachedInGroup[i][j] = 0;}
}
		
for( elem =0; elem< numElem; elem++) {
	if(!averageOnlyActive || element[elem].activeElementFlag) {
		materialGroup = element[elem].getMaterialNumber();
					nnpe = element[elem].numNodesPerElement;
					for(i = 0; i<element[elem].numNodesPerElement;i++) {
						//n = element[elem]->connectivity[i];
						n = element[elem].node[i].nodeNum;
						numberAttachedInGroup[materialGroup][n]++;
					}//loop i
			}// if
		}// loop elem
}//end of numberOfElementsAttachedToNodes   