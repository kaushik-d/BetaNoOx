#include "stdafx.h"

#include "utility/utility.h"
#include "elemWorkspace.hpp"
#include "math/matrix.hpp"

//--------------------------------------------------------------------
ElementWorkspace::ElementWorkspace(void):
extrapMatrices(0)
{
	allocated=false;
	equations=0;
	dofList=0;
	mesh=0;

	quadPointStresses=0;
	rotatedStresses=0;
	extrapolatedStresses=0;

	quadPointStrains=0;
	rotatedStrains=0;
	extrapolatedStrains=0;

	quadPointFluxes=0;
	rotatedFluxes=0;
	extrapolatedFluxes=0;

	quadPointGradients=0;
	rotatedGradients=0;
	extrapolatedGradients=0;

	b=0;
	b1=0;
	b2=0;
	Fe=0;
	initialFe=0;
	Ke=0;
	//FRe=0;
	StoreAllMe=false;
	Me=0;
	MeList=0;
	Ke_original=0;
	F_original=0;

	//set limits for the array sizes
	//SetLimits();
}
//--------------------------------------------------------------------
ElementWorkspace::~ElementWorkspace(void)
{
	int i;
    for(int i=0;i<extrapMatrices.size();++i){
        if(extrapMatrices.at(i)){
            delete extrapMatrices.at(i);
        }
    }
    extrapMatrices.clear();
	if(quadPointStresses){ // FIX 
		for(i=0;i<WorkspaceMaxQuadPoints;i++){
			delete [] quadPointStresses[i];
			delete [] rotatedStresses[i];
			delete [] extrapolatedStresses[i];

			delete [] quadPointStrains[i];
			delete [] rotatedStrains[i];
			delete [] extrapolatedStrains[i];

			delete [] quadPointFluxes[i];
			delete [] rotatedFluxes[i];
			delete [] extrapolatedFluxes[i];

			delete [] quadPointGradients[i];
			delete [] rotatedGradients[i];
			delete [] extrapolatedGradients[i];
		}
		delete [] quadPointStresses;
		delete [] rotatedStresses;
		delete [] extrapolatedStresses;

		delete [] quadPointStrains;
		delete [] rotatedStrains;
		delete [] extrapolatedStrains;

		delete [] quadPointFluxes;
		delete [] rotatedFluxes;
		delete [] extrapolatedFluxes;

		delete [] quadPointGradients;
		delete [] rotatedGradients;
		delete [] extrapolatedGradients;
	}

	if(quadPointStresses){ // FIX 
		for(i=0;i<WorkspaceMaxNumberOfStrains;i++){
			delete [] b[i];
			delete [] b1[i];
			delete [] b2[i];

		}
		delete [] b;
		delete [] b1;
		delete [] b2;
	}

	if(Fe)			delete [] Fe;
	if(initialFe)	delete [] initialFe;
	if(F_original)	delete [] F_original;  

	//if(FRe)		delete [] FRe;

	if(Ke)		delete Ke;
	if(Me)		delete Me;  
	if(MeList)		delete MeList;  
	if(Ke_original)		delete Ke_original; 
	if(dofList) delete [] dofList;
}
//--------------------------------------------------------------------
void ElementWorkspace::PrintLimits(ostream *os)
{
	*os << "WorkspaceMaxQuadPoints		=" <<WorkspaceMaxQuadPoints << endl;
	*os << "WorkspaceMaxStresses		=" <<WorkspaceMaxStresses << endl;
	*os << "WorkspaceMaxFluxes			=" <<WorkspaceMaxFluxes << endl;
	*os << "WorkspaceMaxNodes			=" <<WorkspaceMaxNodes << endl;
	*os << "WorkspaceMaxNumberOfStrains	=" <<WorkspaceMaxNumberOfStrains << endl;
	*os << "WorkspaceMaxNumberDof		=" <<WorkspaceMaxNumberDof << endl;
	*os << endl;
}
//--------------------------------------------------------------------

