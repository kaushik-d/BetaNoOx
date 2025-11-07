#include "stdafx.h"

#include "mesh/BasicMesh.hpp"
#include "utility/utility.h"

//===========================================================
void readQuadPoints(int numEl,vector<Node> &qpList,istream *is)
{
	int dumi;
	double dum = 0.0;
	int count =0;
	int matNum;
	int numQPPerEl;
	for(int i=0; i<numEl; i++){
		(*is)>>dumi;
		(*is)>>matNum;
		(*is)>>numQPPerEl;

		for(int j=0;j<numQPPerEl;j++){
		qpList[count].nodeNum = count ;
		//numdof will be used for mat# of a quad point
		qpList[count].setNumDof(matNum); 

		(*is)>>qpList[count].x;
		(*is)>>qpList[count].y;
		(*is)>>qpList[count].z;
		count++;
		(*is)>>dum;		//this entry is the volume
		}				//associated with the quad point
	}


}
//===========================================================

