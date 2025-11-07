#include "stdafx.h"

#include <iomanip>
#include "BasicMesh.hpp"  

void BasicMesh::printNodalData(double *data, char * label, ostream *outStream )
{
#define FORMAT setw(30)<<setprecision(16)
	(*outStream)<<label<<endl;
	
	for(int i=0;i<numNodes;i++) {
		(*outStream) << i <<"\t";
		(*outStream)<<setiosflags(ios::scientific);
		for(int j=0;j<node[i].getNumDof();j++) {
			(*outStream) << FORMAT<< data[node[i].getFirstDof()+j]  ;
		}
		(*outStream) << endl;
#undef FORMAT
	}
}
//========================================================================

void BasicMesh::printNodalData(double *data, char * label )
{
#define FORMAT setw(15)<<setprecision(6)
	BETA_OUT<<label<<endl;
	
	for(int i=0;i<numNodes;i++) {
		BETA_OUT << i <<"\t";
		BETA_OUT<<setiosflags(ios::scientific);
		for(int j=0;j<node[i].getNumDof();j++) {
			BETA_OUT << FORMAT<< data[node[i].getFirstDof()+j]  ;
		}
		BETA_OUT << endl;
#undef FORMAT
	}
}
//========================================================================
void BasicMesh::printNodalData(double *data, char *label, double tolerance )
{
#define FORMAT setw(15)<<setprecision(6)
	int flag2 = 0;
	BETA_OUT<<label<<endl;
	for(int i=0;i<numNodes;i++) {
		//Check for non-zero data for node
		int flag1 = 0;
		
		for(int j=0;j<node[i].getNumDof();j++) {
			if( fabs(data[node[i].getFirstDof()+j] ) > tolerance )
					{ flag1 = flag2 = 1  ; } 
		}
		
		if(flag1==1)
				{
			BETA_OUT << i <<"\t";
			BETA_OUT<<setiosflags(ios::scientific);
			for(int j=0;j<node[i].getNumDof();j++) {
				BETA_OUT << FORMAT<< data[node[i].getFirstDof()+j]  ;
			}
			BETA_OUT << endl;
		}
#undef FORMAT
	}
	if(flag2==0)
			{BETA_OUT<<"There were no values larger than "<<tolerance<<endl;}
}
//========================================================================
