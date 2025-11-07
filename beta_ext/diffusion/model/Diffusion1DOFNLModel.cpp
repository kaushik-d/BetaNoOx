#include "stdafx.h"

#include <ctime>
#include <string>
#include "Diffusion1DOFNLModel.hpp"

#include "factory/Factory.hpp"
#include "utility/formWriter.hpp"
#include "math/gaussQuadrature/gaussPointList.hpp"
#include "mesh/MeshUtility.h"

#include "../element/TransientElementWorkspace.hpp"
#include "../element/HeatTransferElement3DTransient.hpp"
#include "../element/HeatTransferElement1DTransient.hpp"
//=========================================================================
extern Factory	*factory;
#define FORMAT setw(15)<<setprecision(10)
//-------------------------------------------------------------------
CreateErrorHandler(Diffusion1DOFNLModel);
//--------------------------------------------------------------------
void Diffusion1DOFNLModel::DoTransientNonLinearAnalysis(istream * inStream)
{
	bool DoIterate=true;

	checkModel();
	
	mainTimeCount=0.0;
	mainTimeStepCount=0;

	InitializeOptionalOutput();

	BETA_OUT << "solution:" << endl;
	BETA_OUT << "\n mainTimeStepCount mainTimeCount solution " << endl;

	int numTimeSegments=TimeSegmentList.getNumElements();
	for (int iTimeSegment=0; iTimeSegment<numTimeSegments; iTimeSegment++) {//loop over timesegments
		int numTimeSteps=TimeSegmentList[iTimeSegment].getNumElements();
		for (int iCount=0; iCount<numTimeSteps; iCount++) {//loop over time steps
			timeStepSize=TimeSegmentList[iTimeSegment][iCount];

			AnalysisParameters[0] = time_integration_alpha*timeStepSize;
			AnalysisParameters[1] = (1-time_integration_alpha)*timeStepSize;
			
			mainTimeCount+=timeStepSize;
			mainTimeStepCount++;
			cout << "Current Time= " << mainTimeCount << " seconds." << endl;

			// everytime clear rooms in equation xt: 102200
			equations.zeroLoadVector();
			equations.zeroMatrix(); //xtang 990525
			
			// Set up applied loadVector; Send loads to equation object
			for(int i=0; i<(int)loads.size(); i++){
	            if (loads[i]->isIncremental())
	                loads[i]->setRate(0.01);
				else
					loads[i]->setRate(1.0);  // uniform loading
				loads[i]->applyTo(equations,&mesh->node,mesh->numNodes, this);
			}

			// apply additional equations: important for dummy load !
			equations.applyAdditionalEquations(); 

			//linear solution:
			sumVector(u_s1,		u_s,		delta_u,totalNumDof);
			assembleK_F_s_alpha_family(u_s, delta_u, u_s1, AnalysisParameters, timeStepSize); 

//			displayTime("Time till solve", BETA_OUT);
			double * deltaSolution=equations.solve(u_s1);
//			displayTime("Time for solve", BETA_OUT);

			//equations.addIncrementalDisplacements(deltaSolution);
			copyVector(deltaSolution, delta_u, totalNumDof); //copy deltaU to local array
			sumVector(u_s1,		u_s,		delta_u,totalNumDof);

			TimeStepUpdate();//update Phi

			BETA_OUT << "\n" << mainTimeStepCount << " " << mainTimeStepCount << " " << FORMAT << u_s1[1];

			if(DoIterate) TransientAnalysisIterate();

			WriteOptionalOutput();

			//BETA_OUT << iCount << " " << count << " " << u_s1[1] << endl;

			initializeVector( delta_u, 0.0, totalNumDof);
			copyVector( u_s1, u_s, totalNumDof);
		}//end TimeStep loop
		BETA_OUT << endl;
	}//end TimeSegment loop
	BETA_OUT << endl;

	FinalizeOptionalOutput();
} // end of DoTransientNonLinearAnalysis
//========================================================================
void Diffusion1DOFNLModel::TransientAnalysisIterate()
{
	int iIter=0;
	int nIter = 50;  //fix this later
    double * deltaSolution;

	//for calculating residual
    double maxResidual;
    int    dofForMaxResidual;
    double * mappedResidualVector=0; // = new double[totalNumDof]; //do not modify mappedResidualVector
	//only for 'reading', because mappedResidualVector points to equation object's data

	getGlobalForcesAndResidual(internalForce,
					AnalysisParameters, timeStepSize,
					u_s, delta_u, u_s1, 
					mappedResidualVector,maxResidual,dofForMaxResidual );

//	BETA_OUT << " (" << mappedResidualVector[1] << ") ";

	while ( (maxResidual > maxAllowableResidual) && (iIter <= nIter) ){

		assembleK_alpha_family(u_s, delta_u, u_s1, AnalysisParameters, timeStepSize); 

		deltaSolution=equations.incrementalSolutionUsingResidual(u_s1);

		equations.addIncrementalDisplacements(deltaSolution);
		sumVector(delta_u,	delta_u,	deltaSolution,totalNumDof);
		sumVector(u_s1,		u_s,		delta_u,totalNumDof);

		BETA_OUT << " " << u_s1[1];
		//mesh->printNodalData(u_s1, "displacements", &BETA_OUT );

		getGlobalForcesAndResidual(internalForce,
						AnalysisParameters, timeStepSize,
						u_s, delta_u, u_s1, 
						mappedResidualVector,maxResidual,dofForMaxResidual );

//				BETA_OUT << " (" << mappedResidualVector[1] << ") ";
    
		iIter++;
	}//end of iterate loop

	cout << "Interations taken: " << iIter << endl;
	if(iIter > nIter){ FatalError("Exceeded allowed number of iterations");}
}
//========================================================================
