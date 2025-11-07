#include "stdafx.h"
#include "BoundaryNodeList.hpp"
#include "elements/BasicElement.hpp"
//====================================================================================
int getNumFaces(BasicElement * ele){
	switch (ele->elementType) {
	case Q2D8N : 
		return 8;
		break;
	case Q2D4N : 
		return 4;
		break;
	case Q3D8N : 
		return 8;
		break;
	case Q3D4N : 
		return 4;
		break;
	case H3D20N : 
		return 6; //have to think about this more...
		break;
	case T3D10N : 
		return 4;
		break;
	default: return 0;
	}
}
//====================================================================================
void getFaceNodes(BasicElement * ele, int faceID, FaceNodes * faceNodes)
{
	int i;

	faceNodes->face=faceID;

	switch (ele->elementType) {
	case Q2D8N : 
		faceNodes->numNodesPerFace=2;
		faceNodes->faceType=1;  // line
		for (i=0; i<faceNodes->numNodesPerFace; i++) 
	        faceNodes->node[i]=&ele->node[cQ2D8N[faceID][i]];
		break;
	case Q2D4N : 
		faceNodes->numNodesPerFace=2;
		faceNodes->faceType=1;  // line
		for (i=0; i<faceNodes->numNodesPerFace; i++) 
	        faceNodes->node[i]=&ele->node[cQ2D4N[faceID][i]];
		break;
	case Q3D8N : 
		faceNodes->numNodesPerFace=2;
		faceNodes->faceType=1;  // line
		for (i=0; i<faceNodes->numNodesPerFace; i++) 
	        faceNodes->node[i]=&ele->node[cQ3D8N[faceID][i]];
		break;
	case Q3D4N : 
		faceNodes->numNodesPerFace=2;
		faceNodes->faceType=1;  // line
		for (i=0; i<faceNodes->numNodesPerFace; i++) 
	        faceNodes->node[i]=&ele->node[cQ3D4N[faceID][i]];
		break;
	case H3D20N : 
		faceNodes->numNodesPerFace=8;
		faceNodes->faceType=2;  // 2D face
		for (i=0; i<faceNodes->numNodesPerFace; i++) 
	        faceNodes->node[i]=&ele->node[cQ3D20N[faceID][i]];
		break;
	case T3D10N:
		faceNodes->numNodesPerFace=6;
		faceNodes->faceType=2;
		for (i=0; i<faceNodes->numNodesPerFace; i++) 
	        faceNodes->node[i]=&ele->node[cT3D10N[faceID][i]];
		break;
	default: break;
	}
}
//====================================================================================
BasicElement * findNeighborElement(BasicElement * ele, FaceNodes * faceNodes)
{
	int i, j, k, minMeet=1000000;
	Node *n, *nj;
	BasicElement * e;
	bool match=true, hasEle;

	for (i=0; i<faceNodes->numNodesPerFace; i++) {
		if (faceNodes->node[i]->numElementsAtNode < minMeet) {
			minMeet=faceNodes->node[i]->numElementsAtNode;
			n=faceNodes->node[i];
		}
	}
	for (i=0; i<minMeet; i++) {
		e=n->eList[i];
		if (e==ele) continue;
		match=true;
		for (j=0; j<faceNodes->numNodesPerFace; j++) {
			nj=faceNodes->node[j];
			if (n==nj) continue;
			hasEle=false;
			for (k=0; k<nj->numElementsAtNode; k++) {
				hasEle = hasEle || (e == nj->eList[k]);
			}
            match = match && hasEle;
		}
		if (match) return e; 
	}
	return 0;
}
//====================================================================================
BoundaryNodeList * getMeshBoundaryNodes(ElementGroup *ele, int numElements, int numNodes)
{
	
	int i,j,k;
	BasicElement *ei;
	FaceNodes fNodes;
	
	Node **nDump;
	nDump = new Node * [numNodes];

	bool isIn;
	int lNodes=0, l, numFaces;

	for (i=0; i<numElements; i++) {
		ei=&(*ele)[i];

		if (!ei->activeElementFlag) continue; 

		numFaces=getNumFaces(ei);
		for (j=0; j<numFaces; j++) {
			getFaceNodes(ei, j, &fNodes);
			if (findNeighborElement(ei, &fNodes)==0) {// the face is a boundary face
				for (k=0; k<fNodes.numNodesPerFace; k++) {
					isIn=false;
					for (l=0; l<lNodes; l++) {
						if (nDump[l]==fNodes.node[k]) {
							isIn=true;
							break;
						}
					}
					if (!isIn) {
						nDump[lNodes]=fNodes.node[k];
						lNodes++;
					}
				}
			}
		}
	}

	if (lNodes==0) return 0;
	
	BoundaryNodeList * bNodeList = new BoundaryNodeList();
	bNodeList->numNodesInList=lNodes;
	bNodeList->nList=new Node *[lNodes];
	for (l=0; l<lNodes; l++) bNodeList->nList[l]=nDump[l];

	delete [] nDump;
	return bNodeList;
}
