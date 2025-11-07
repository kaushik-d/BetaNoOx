#include "stdafx.h"

#include "utility/utility.h"
#include "utility/FileManager.hpp"
#include "SingleFieldMacroElement.hpp"

#include "utility/excepts.hpp"
#include "math/matrix.hpp"
#include "math/equation/equation.hpp"
#include "factory/Factory.hpp"

extern int verboseFlag;
extern Factory	*factory;

void deriv3D(int numberOfInterp, double xi,  double eta, double zeta, double * p_S_Lx1, double * p_S_Lx2, double * p_S_Lx3);
void shape3D(int numberOfInterp, double xi,  double eta, double zeta, double * S);
bool readScript(Factory *factory, BasicModel* &modelObject, string filename, ifstream *& inStream);
//=====================================================================
SingleFieldMacroElement::SingleFieldMacroElement(void)
{
}
//=====================================================================
SingleFieldMacroElement::~SingleFieldMacroElement()
{
}
//=====================================================================
void SingleFieldMacroElement::CalculateK(double *uGlobal,string filename)
{
	BasicMesh *subdomains=model->mesh;

	IntegrationDataType   idt;
	idt.id=new InterpAndSolDerivatives(numNodesPerElement,numNodesPerElement); idt.id->allocate();

	idt.elemDisp[0] = idt.id->nodalValues_u1; 
	idt.elemDisp[1] = idt.id->nodalValues_u2;
	idt.elemDisp[2] = idt.id->nodalValues_u3;     
	idt.id->numberOfInterp = numNodesPerElement;
	idt.uGlobal=uGlobal;

  subdomains->PrintMesh(filename,false); //ksm 6/29/10

	extractNodalCoordinates(idt.id->xCoor,idt.id->yCoor,idt.id->zCoor);
	if( uGlobal != NULL) {extractSolution( uGlobal, idt.elemDisp);} 
	getQuadraturePoints(idt.gaussPointList);

	initializeArrays();

	double integrationFactor;
	double volume =0.;
	double weight;

	int numberofsubdomains=subdomains->numElements;
	for (int subdomain=0; subdomain < numberofsubdomains; subdomain++){
		ElasticityElement3D* element=(ElasticityElement3D*)(&subdomains->element[subdomain]);

		IntegrationDataType   SDidt;
		SDidt.id=new InterpAndSolDerivatives(element->numNodesPerElement,element->numNodesPerElement); SDidt.id->allocate();
		SDidt.id->numberOfInterp = element->numNodesPerElement;

		element->extractNodalCoordinates(SDidt.id->xCoor,SDidt.id->yCoor,SDidt.id->zCoor);

		material=element->getMaterial();//JV092205 the macro element has 'dynamic' material properties
		setMaterialGroup(element->getMaterialNumber());

		getQuadraturePoints(SDidt.gaussPointList);
		int SDtotalNumIPs = SDidt.gaussPointList.totalNumIPs;

		for(int ip=0;ip<SDtotalNumIPs;ip++)
		{
			element->interpolation(ip, SDidt.gaussPointList, SDidt.id);
			element->jacobian(SDidt.id);
			weight = SDidt.gaussPointList.weightList[ip];

			double xi=0.;
			double eta=0.;
			double zeta=0.;
			for(int kk=0; kk < element->numNodesPerElement; kk++){
				xi   += SDidt.id->S[kk] * SDidt.id->xCoor[kk];
				eta  += SDidt.id->S[kk] * SDidt.id->yCoor[kk];
				zeta += SDidt.id->S[kk] * SDidt.id->zCoor[kk];
			}

			deriv3D(idt.id->numberOfInterp, xi, eta, zeta, idt.id->p_S_Lx1, idt.id->p_S_Lx2, idt.id->p_S_Lx3);
			shape3D(idt.id->numberOfInterp, xi, eta, zeta, idt.id->S);
			jacobian(idt.id);

			integrationFactor = SDidt.id->detJ * idt.id->detJ * weight * getIntegrationFactor();
			volume += integrationFactor;


            //element->getOrientationAtIP()->interpolate_rotation(SDidt.id,element->getElementNodeRotations());
            element->getOrientationAtIP()->interpolate_rotation_by_angle(SDidt.id,element->getElementNodeRotations());
            //element->getOrientationAtIP()->interpolate_old_beta(SDidt.id,element->getElementNodeRotations(),material->getMaterialRotation());
			ElasticMaterial* material = (ElasticMaterial*)BasicElement::material;
			material->rotateMaterialToGlobal(element->getOrientationAtIP());

			if(idt.uGlobal!=NULL){derivativesOfFields(idt.id);}
			calculateBmatrix(idt.id,b);

			Matrix *cMat = material->getPointerToGlobalCmat();
			double *stress	= bag->quadPointStresses[ip];
			addToKe(cMat,idt.id,integrationFactor,stress);
		}  // end of integration loop
	} // end of submesh loop !!!
	checkForNegativeDiagonal(); 
	if(MasterKe==0){
		MasterKe=new Matrix(*Ke,numDof,numDof);
		MasterKe->print("MasterKe SDI = ",&BETA_OUT);
//		MasterKe->convertSymmetricToGeneral();
//		MasterKe->printInMaple("MasterKe SDI Maple= ",numDof,numDof,&BETA_OUT);
	}

}
//==========================================================
void SingleFieldMacroElement::update(char * action, double *uGlobal)
{
	if(strstr(action,"K") !=NULL) {
		if(model)
			Ke=modelOwner->MasterKe;
		else
			EXIT_BETA("you cannot call SingleFieldMacroElement::update() before the macroelement has been defined.");
	}else if(strstr(action,"I") !=NULL) {
		for(int i=0;i<numDof;i++) initialFe[i] = 0.;
	}else{
		cout << "calling SingleFieldMacroElement::update() with '"<< action << "'"<< endl;
		//EXIT_BETA("SingleFieldMacroElement::update() can currently calculate only the stiffness matrix");
    cout<<"WARNING: Incorrect stress output for macroelements...see SingleFieldMacroElement.cpp"<<endl; //ksm6-29-10
	}

}
//==========================================================
void SingleFieldMacroElement::readSpecialCommand(istream &inStream, ElementGroup *element, int numElements, char * command)
{
BETA_OUT<<"SingleFieldMacroElement::readSpecialCommand"<<endl;

if (COMPARE(command, "DefineMacroElementGroups")==0 ) {
	//calculate masterKe using first element in the ElementGroup
	SingleFieldMacroElement* master=(SingleFieldMacroElement*)(&(*element)[0]);

	BasicModel *modelObject =0; 
	ifstream *subdomainstream=0;
	string filename;
	inStream >> filename;
	readScript(factory, modelObject, filename, subdomainstream);
	displayTime("Time to finish running Subdomain script", *modelObject->filemanager->outStream);
	if(subdomainstream) delete subdomainstream;

	master->model=modelObject;
	master->modelOwner=master;
	modelObject->createAssembler();
	modelObject->allocateForAnalysis(BETA_DO_NOT_ALLOCATE_EQUATIONS);
	master->setPointers(modelObject->assembler->threadWorkspaceList[0].eWorkspace);//using the subdomain model's workspace to calculate the master Ke
	master->calculateDofList();
	master->CalculateK(NULL, filename);

	//linking all the other elements to the subdomain model as well
	for(int i=1;i<numElements;i++) { 
		SingleFieldMacroElement* e=(SingleFieldMacroElement*)(&(*element)[i]);
		e->model=modelObject;
		e->modelOwner=master;
	}
	return;
}
// If no match, try super class
BETA_OUT<<"No match... try super class"<<endl;
MacroElement::readSpecialCommand(inStream,element,numElements, command);
}//readSpecialCommand
//==========================================================
