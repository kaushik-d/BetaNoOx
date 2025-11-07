#include "stdafx.h"
#include "mesh/BasicMesh.hpp"
#include "math/equation/equation.hpp"


void ReadMPCFile(string filename, BasicMesh &mesh, Equations &equations)
{
	ifstream mpcfile(filename.c_str());
	if(!mpcfile.good()){
		cerr << "error with file in ReadMPCFile()" << endl;
		exit(1);
	}
	cout << "Reading MPC File..." << endl;

	char buffer[256];
	mpcfile.getline(buffer,256);

	int slavenn,slavedof,slaveglobaldof,numMasters,masterNN[8];
	double masterWeight[8],diff;
	
	AdditionalEquation *eqn;
	Node *sn,*mn;


	int i;
	while (!mpcfile.eof()){
		mpcfile >> slavenn ;
		if(slavenn < 0) break;
		mpcfile >> slavedof >> slaveglobaldof >> numMasters;
		for(i=0;i<numMasters;i++){
			mpcfile >> masterNN[i] >> masterWeight[i];
		}
		mpcfile >> diff;	

		sn=&mesh.node[slavenn];

		eqn = new AdditionalEquation;
		eqn->setSlave(sn->getFirstDof() + slavedof );
		for(i=0;i<numMasters;i++){
			mn=&mesh.node[masterNN[i]];
			eqn->addAdditionalTerm(mn->getFirstDof()+ slavedof, masterWeight[i]);
		}
		eqn->setConstant(diff);
		
		//cout <<  slavenn << " " <<  slavedof << endl;
		equations.addAdditionalEquation(eqn);
	}
	mpcfile.close();
	cout << "Finished reading MPC File" << endl;
}