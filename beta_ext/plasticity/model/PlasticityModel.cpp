#include "stdafx.h"

#include <ctime>
#include <string>
#include "PlasticityModel.hpp"

#include "factory/Factory.hpp"
#include "utility/formWriter.hpp"
#include "math/gaussQuadrature/gaussPointList.hpp"
#include "mesh/MeshUtility.h"

#include "elements/3D/ElasticityElement3D.hpp"
#include "plasticity/element/plasticElement.hpp"
//=========================================================================
//extern FileManager filemanager;
extern Factory	*factory;
void deleteFiles(string files);
//-------------------------------------------------------------------
CreateErrorHandler(PlasticityModel);
//--------------------------------------------------------------------
PlasticityModel::PlasticityModel(void)
{
	NodalLCS_stressRequired=false;//-DG_Oct24_2006
	stepIncrementsNodalLCS_stress=1;//-DG_Oct24_2006
	NodalGCS_stressRequired=false ;//-DG_Oct24_2006
	stepIncrementsNodalGCS_stress=1;//-DG_Oct24_2006
	Quad_stressRequired=false ;//-DG_Oct24_2006
	stepIncrementsQuad_stress=1;//-DG_Oct24_2006
	Damage_Factors_Required=false; //-DG_March29_2007
	stepIncrementsDamage_Factors=1; //-DG_March29_2007
}
//--------------------------------------------------------------------
PlasticityModel::~PlasticityModel(void)
{
}
//--------------------------------------------------------------------
int PlasticityModel::processCommand(char** localTokenList,ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;
	char newFileName[80];

	char token[200];
	int  foundMatch=0;

strcpy(token,localTokenList[0]);
BETA_OUT<<"\nCommand processor = 'PlasticityModel::processCommand'          "<<token<<endl;

	if( COMPARE(token,"openFile")==0 ) {   //This is for this pass only
		//Should be 3 tokens
		strcpy(newFileName,localTokenList[1]);
		currentStream = openFileForInput(newFileName);
		strcpy(token,localTokenList[2]);
		BETA_OUT<<"Input stream for command ***"<<token<<"*** = "<<newFileName<<endl;
		BETA_OUT<<"Stream(openFile) = "<<currentStream<<endl;
		BETA_OUT<<"Stream(openFile) = "<<(*currentStream)<<endl;  
	}

	//----------- Insert new commands below this line ----------------------
	IS_COMMAND(DoPlasticAnalysis)  
	//----------- Insert new commands above this line ----------------------

	if(currentStream!=originalStream) {
		BETA_OUT<<"Close alternate input stream"<<endl;
		(*currentStream).close();
		delete currentStream;
		currentStream=originalStream;
		BETA_OUT<<"Automatically reverts to default stream "<<currentStream
			<<" upon return to readCommands"<<endl;
	}

	if(foundMatch ==0){	
		foundMatch=DamageModel::processCommand(localTokenList,currentStream, numToken);
	}

	return foundMatch;
}
//======================================================================
void PlasticityModel::ProcessSetElementProperty(string command, istream * inStream)
{

	int intCommand=0 ;
	enum {
		SET_ELEMENT_ISV_gauss,  // +xtang 991031 
	};

	ChangeToUpper(command)
	if(command=="SETELEMENTISV_GAUSS")
		{intCommand=SET_ELEMENT_ISV_gauss;}
	else {
		return DamageModel::ProcessSetElementProperty(command, inStream);
	}

int i,first,last,increment;
int isvType;
int matID;

while (2==2) {				
	(*inStream)>>first>>last>>increment;
 	    //BETA_OUT<<"First, last, increment="
		//	<<first<<" "<<last<<" "<<increment<<endl;
	if(first<0) return;

	//Error checking...
	if(first < 0 || first    > mesh->numElements || first > last
				 || increment<1 )
	  {
		BETA_OUT<<"Fatal error...input is bad"<<endl;
	    BETA_OUT<<"First, last, increment="
			<<first<<" "<<last<<" "<<increment<<endl;
	    exit(1);}				


	BasicElement *e=0;

switch(intCommand) {
		//xtang 122000	
		case SET_ELEMENT_ISV_gauss:
			{
			(*inStream) >> isvType;	
			BETA_OUT << "\tisvType = " << isvType << " for materials  "
  			<< first <<" to "<<last << " by " << increment << '\n';				
			for (int jj=first; jj<=last; jj+=increment) {  // material 
				for(i=0;i<mesh->numElements;i++) {
					e=&mesh->element[i];
					matID=e->getMaterialNumber();
					if (matID == jj) {
						((PlasticElement *) e)->isv=((PlasticMaterial *) (&materialList[matID]) )->createISVSpace(isvType,27); // hardwired !!! 3x3x3=gauss points
					}
				}
			}
			}
			break;

	default:
		BETA_OUT<<"Match not found"<<endl;
		exit(1);
}

}//end of while (2==2)
} //end of ProcessSetElementProperty
//======================================================================
void PlasticityModel::DoPlasticAnalysis(istream * inStream)
{
    
    int i;
    double maxResidual;
    int    dofForMaxResidual;
    double * mappedResidualVector = new double[totalNumDof];
    double * totalForce = new double[totalNumDof];

	// xtang: note that model.currentSolution is a unmapped vector,
	//        equ.solution is generally a unmapped vector

	double * mappedCurrentSolution = new double[totalNumDof];  

    

//    double * currentSolutionInEquations;
    double * deltaSolution;

	// double * dq=equations.getResultantVector();
    // double * df=equations.getLoadVector();

	

    int iIter;
    int count=1;

    ofstream  os, os1;

    // for sparse doc.
    //os.open("sparse.dat");
    //os.close();

    double aRate=1.0/double(MaxAllowableLoadSteps);

	
	

	char aName[200];
	//ofstream ost;

 //   FormWriter fr(OUT);

	allocateForAnalysis();
    
    // initializing solution vectors
    solution = new double[totalNumDof];
	double * dump = new double[totalNumDof];
	double * dq=equations.getResultantVector();
    double * df=equations.getLoadVector();

	for (i=0;i<totalNumDof;i++) {solution[i]=totalForce[i]=0.0;}

    // "Load loop" begins here
	// xtang 990720
	PlasticElement *e=(PlasticElement*)(&mesh->element[0]);
//	PlasticElement::openOutputFiles(filemanager, OutputOptions);

	PlasticElement::eVolDist->seekp(0,ios::beg);
	//setPrintStressFlag(0);

	// xtang 052801: The internal force (resultantLoad) is obtained by using stresses 
	// stored at ISV.
	for(i=0; i<mesh->numElements; i++){
		e=(PlasticElement*)(&mesh->element[i]);
		e->update("A", currentSolution);
	} // update ISV (get stress increment only !!!
    

unloading:

	do { 
        // everytime clear rooms in equation xt: 102200
		equations.zeroLoadVector();
		equations.zeroMatrix(); //xtang 990525
		equations.zeroSolution();  
		equations.zeroResultantVector();  
		for (i=0;i<totalNumDof;i++) {mappedCurrentSolution[i]=0.0;}


        // Set up applied loadVector; Send loads to equation object
		for( i=0; i<(int)loads.size(); i++){
            if(loads[i] == NULL) break;
            if (loads[i]->isIncremental()) { 
                loads[i]->setRate(aRate);
            }
            loads[i]->applyTo(equations,&mesh->node,mesh->numNodes, this);
		}
        
        // apply additional equations: important for dummy load !
        equations.applyAdditionalEquations(); 

        // switch Primary BC
//        currentSolutionInEquations = equations.getSolution();  // unmapped
//        copyVector(currentSolutionInEquations, currentSolution, totalNumDof);  
        copyVector(equations.getSolution(), currentSolution, totalNumDof);  

        iIter=0;

		// xtang 052801: The internal force (resultantLoad) is obtained by using stresses 
		// stored at ISV.
		for(i=0; i<mesh->numElements; i++){
			e=(PlasticElement*)(&mesh->element[i]);
			e->update("V", currentSolution);
		} // update ISV (get stress increment only !!!

		// get initial solution for the current load increment;
		assembler->assembleKandF(currentSolution); // K_tangent and initial Force
        getGlobalForces(internalForce,currentSolution, // d_q
           mappedResidualVector,maxResidual,dofForMaxResidual );

        deltaSolution=equations.incrementalSolutionUsingResidual(currentSolution); // mapped
//		for (i=0;i<totalNumDof;i++) {mappedCurrentSolution[i]+=deltaSolution[i];} // mapped

        equations.addIncrementalDisplacements(deltaSolution);  // d_a+=delt_a : unmapped
        
//		currentSolutionInEquations = equations.getSolution();  //  unmapped
//      copyVector(currentSolutionInEquations,currentSolution,totalNumDof); // unmapped
        copyVector(equations.getSolution(),currentSolution,totalNumDof); // unmapped
       
        do {  // iteration begins here
			for(i=0; i<mesh->numElements; i++){
				e=(PlasticElement*)(&mesh->element[i]);
				e->update("V", currentSolution);
			} // update ISV (get stress increment only !!!

            assembler->assembleKandF(currentSolution); 
            getGlobalForces(internalForce,currentSolution, 
                mappedResidualVector,maxResidual,dofForMaxResidual );
            
            if(maxResidual< (2.0*maxAllowableResidual)) { break;}

            deltaSolution=equations.incrementalSolutionUsingResidual(currentSolution);
            equations.addIncrementalDisplacements(deltaSolution);
            
//            copyVector(currentSolution,solution,totalNumDof); // U_r after solving
//          currentSolutionInEquations = equations.getSolution();
//          copyVector(currentSolutionInEquations,currentSolution,totalNumDof);  // U_r+1 after soving
            copyVector(equations.getSolution(),currentSolution,totalNumDof);  // U_r+1 after soving
			iIter++;

		} while (iIter <MaxAllowableIterations); 

		for (i=0;i<totalNumDof;i++) {
			solution[i]+=currentSolution[i]; // unmapped 
			totalForce[i]+=internalForce[i]; // xi*(df[i]); unmapped
		}

		// update element ISV

	    updateElementISV();


        BETA_OUT<<iIter<< "  "<<maxResidual<<"  ";

		// how to output?????????????
        //outVolumeStr11(internalForce, currentSolution, volS1, volS2);  // volS1=volS2=1 by default xtang 990713

        outVolumeAverageISV(); // 052801

		// xtang 052901: print displacement
		sprintf(aName,"disp%d.dat",count);
		ofstream *osdg;
		//ost.open(aName);
		osdg = filemanager->OpenOutputStream(aName);
		mesh->printNodalData(solution, "displacements", osdg);
		//osdg->close();
		filemanager->CloseOutputStream(osdg);

		// xtang 052901: print element stress
		sprintf(aName,"%d.LCS.stress",count);
		//ost.open(aName);
		osdg = filemanager->OpenOutputStream(aName);
		printElementStress("stress", osdg);
		//osdg->close();
		filemanager->CloseOutputStream(osdg);

		//JV071107 - why is this block here ??? 
		//it is printing the same thing as the previous block !!!
		/*
		// xtang 060901
		sprintf(aName,"q_stress%d.dat",count);
		//ost.open(aName);
		osdg = filemanager.OpenOutputStream(aName);
		printElementStress("stress", osdg);
		osdg->close();
		*/

//====
//quick fix for printing quad stresses.. will need to implement a neat output routine 
		{
//			if(count==1 || count==stepIncrementsQuad_stress*kkkk)
//			{
				char cdumpDG[500];
				char outBase[] = {'Q','u','a','d','s','t','r','e','s','s','\0'};
				ofstream *osdg;
				sprintf(cdumpDG,"%d.%s",count,outBase);  // like "1.o.Quad.stress"
				//osdg = new ofstream(); 
				//osdg->open(cdumpDG);
				osdg = filemanager->OpenOutputStream(cdumpDG);
				
				GaussPointList  gaussPointList;
				int totalNumIPs;
								
				equations.zeroResultantVector();	
				for(i=0; i<mesh->numElements; i++)
				{
					e=(PlasticElement*)(&mesh->element[i]);
					int group = e->getMaterialNumber();		
					e->update("F", currentSolution);
					e->getQuadraturePoints(gaussPointList);
					totalNumIPs=gaussPointList.totalNumIPs;
					e->printQuadStresses(totalNumIPs, *osdg);
				}
				//osdg->close();
				filemanager->CloseOutputStream(osdg);
				//if(count==1 && stepIncrementsNodalGCS_stress==1){kkkk=kkkk+1;}
				//kkkk++;
//			}
		}
//====

        count++;

    } while (count <= MaxAllowableLoadSteps);

	if (MaxAllowableIterations < 0) {
		MaxAllowableIterations=-MaxAllowableIterations;
		aRate=-aRate;
		goto unloading;
	}
                
    mesh->printNodalData(currentSolution, "Nodal Displacements" );    
    
	// +xtang 990720
    // delete "_stiff.dat" and "_quadinfo.dat"
	
	//eMF->close();
	//eQuadInfo->close();
	deleteFiles("_stiff.dat");
	deleteFiles("_quadInfo.dat");


    // keeep these two
	//eVolDist->close();
	//((fstream*)PlasticElement::eVolDist)->close();
//	PlasticElement::closeOutputFiles(filemanager, OutputOptions);
	//eDamage->close();
	// -xtang 990720

    delete [] mappedResidualVector; //xt062901
    delete [] totalForce; //xt062901
	delete [] mappedCurrentSolution; //xt062901
    delete [] solution; //xt062901
	delete [] dump; //xt062901
	
} // end of DoPlasticAnalysis

#undef F1
#undef F2
//=========================================================
void PlasticityModel::updateElementISV()
{
    int i;
    for(i=0;i<mesh->numElements;i++) { mesh->element[i].update("U",0);}
}
//=========================================================
void PlasticityModel::outVolumeAverageISV()
{
#define maxStresses 6
#define FORMAT setw(15)<<setprecision(6)
    
	double *volAvgStress, 
		   *volAvgStrain, 
		   *volAvgPlasticStrain,
		   effectiveStress=0.0,
		   effectivePlasticStrain=0.;

    double Totalvolume=0., volume=0.0;

    double *TotalvolAvgStress, 
		   *TotalvolAvgStrain,
	       *TotalvolAvgPlasticStrain, 
		   TotalEffectiveStress=0.0,
		   TotalEffectivePlasticStrain=0.0;

    int i, is;

    volAvgStress = new double [maxStresses];
    volAvgStrain = new double [maxStresses];
    volAvgPlasticStrain = new double [maxStresses];
 
	TotalvolAvgStress = new double [maxStresses];
    TotalvolAvgStrain = new double [maxStresses];
    TotalvolAvgPlasticStrain = new double [maxStresses];

    for(i=0;i<maxStresses;i++){
        TotalvolAvgStress[i]=0.0;
        TotalvolAvgStrain[i]=0.0;
        TotalvolAvgPlasticStrain[i]=0.0;
    }

	PlasticElement* e=0;
    // Calculate Volume
    for(i=0;i<mesh->numElements;i++) {
		e=(PlasticElement*)(&mesh->element[i]);
        
        e->calculateDofList();
        volume=0.;
        ((PlasticElement * ) e)->getVolAvgISV(
			volAvgStress,volAvgStrain, volAvgPlasticStrain, 
			effectiveStress, effectivePlasticStrain, volume);      

        Totalvolume       += volume;
		TotalEffectiveStress += effectiveStress;
		TotalEffectivePlasticStrain += effectivePlasticStrain;

        for(is=0; is<maxStresses; is++){
           TotalvolAvgStress[is] += volAvgStress[is];
           TotalvolAvgStrain[is] += volAvgStrain[is];
           TotalvolAvgPlasticStrain[is] += volAvgPlasticStrain[is];
        }
    }
    #define FF1 setw(7)<<setprecision(4)<<setiosflags(ios::fixed)

    for(is=0; is<maxStresses; is++){
		BETA_OUT<<FORMAT<<TotalvolAvgStrain[is]/Totalvolume;
		BETA_OUT<<FORMAT<<TotalvolAvgStress[is]/Totalvolume; 
		BETA_OUT<<FORMAT<<TotalvolAvgPlasticStrain[is]/Totalvolume; 
	}

	BETA_OUT<<FORMAT<<TotalEffectivePlasticStrain/Totalvolume; 
	BETA_OUT<<FORMAT<<TotalEffectiveStress/Totalvolume;
    BETA_OUT<<endl;

    delete [] volAvgStress; //xt062901
    delete [] volAvgStrain; //xt062901
    delete [] volAvgPlasticStrain; //xt062901
 
	delete [] TotalvolAvgStress; //xt062901
    delete [] TotalvolAvgStrain; //xt062901
    delete [] TotalvolAvgPlasticStrain; //xt062901

}
#undef FORMAT
//====================================================
void PlasticityModel::printElementStress(char * label, ostream *outStream )
{
	(*outStream) <<label<<endl;
	int iele;
	PlasticElement* e=0;
	for (iele=0; iele<mesh->numElements; iele++) {
		e=(PlasticElement*)(&mesh->element[iele]);
		((PlasticElement *) e)->printElementalStresses(outStream);
	}
}
