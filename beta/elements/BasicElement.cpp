#include "stdafx.h"

#include "utility/utility.h"
#include "utility/FileManager.hpp"
#include "BasicElement.hpp"

#include "utility/excepts.hpp"
#include "math/matrix.hpp"
#include "math/equation/equation.hpp"

extern int verboseFlag;
//=====================================================================
BasicElement::BasicElement(void):
isv(0)
{
	dofList = NULL;
	numDof=0;
	materialGroup= -1;
	activeElementFlag=true;
	nodalValues=0;
	bag=0;
	equations=0;
	Ke=0;
	Fe=0;
	initialFe=0;
	material=0;
	elementType=0;
}

BasicElement::~BasicElement()
{
	node.clear();
	if(nodalValues) delete [] nodalValues; nodalValues=0;

	if(requiresParallelESAConsideration){
		if(Ke) delete Ke; Ke=0;
		if(Fe) delete [] Fe; Fe=0;
		if(initialFe) delete [] initialFe; initialFe=0;
		if(dofList) delete [] dofList; dofList=0;
	}
    if(!isv.empty()){
        if(isv[0]){
            material->deallocateISVs(isv);
        }
    }
}
//---------------------------------------------------------------------
//WARNING: doflist is shared (dynamic)- the data is stored in the elementWorkspace
int *BasicElement::calculateDofList()
{
	int i,j,k;

// Count numDof
	numDof = 0;
	for(i=0;i<numNodesPerElement;i++)
		numDof+=node[i].getNumDof();

	k=0;
	for(i=0;i<numNodesPerElement;i++)
		for(j=0;j<node[i].getNumDof();j++)
			dofList[k++]=node[i].getFirstDof()+j;

			//BETA_OUT<< "DOF list  "<<numDof<<"\n";
			//for(int iii=0; iii<numDof; iii++)
			 //   { BETA_OUT<< dofList[iii]<<" ";}
			//	 BETA_OUT<<endl;
	return dofList;
}
//=====================================================================
void BasicElement::setElementType(int nDims)    // xtang 042902: hardwired...fix later!
{
	switch(nDims) {
		case -2:if      (numNodesPerElement==8)  { elementType=Q2D8N; }
				else if (numNodesPerElement==4)  { elementType=Q2D4N; }
				else { elementType=DEF2D; }
				break;
		case 1: 
				if      (numNodesPerElement==2)  { elementType=L1D2N; }
				else { elementType=DEF1D; }
				break;
		case 2: if      (numNodesPerElement==8)  { elementType=Q2D8N; }
				else if (numNodesPerElement==4)  { elementType=Q2D4N; }
                else if (numNodesPerElement==6)  { elementType=T2D6N; }
                else if (numNodesPerElement==3)  { elementType=T2D3N; }
				else { elementType=DEF2D; }
				break;
		case 3: if      (numNodesPerElement==20) { elementType=H3D20N; }
				else if (numNodesPerElement==8)  { elementType=H3D8N; }
				else if (numNodesPerElement==10) { elementType=T3D10N; }
				else if (numNodesPerElement==4)  { elementType=T3D4N; }
				else if (numNodesPerElement==15) { elementType=W3D15N; }
				else if (numNodesPerElement==6)  { elementType=W3D6N; }

				else { elementType=DEF3D; }
				break;
		default:break;
	}
}
//=====================================================================
void BasicElement::specifyNonZeroLocation()
{
	if(equations==0) {
		cout << "equations pointer not set!" <<endl;
		cout << "must be set before calling specifyNonZeroLocation()!" <<endl;
		
		exit(1);
	}

	int *map,elemRank;

 	map = calculateDofList();
	elemRank = getNumDof();
	equations->specifyNonZeroLocation(map,elemRank);
}
//=====================================================================
void BasicElement::checkForNegativeDiagonal()
{
  //Check for negative diagonals
	 for(int i=0; i<numDof; i++){
		 if((*Ke)[i][i]<=0) 
			{	BETA_OUT<<"BasicElement::bad diagonal in K - element No." << elementNumber <<endl;
				cerr <<"BasicElement::bad diagonal in K - element No." << elementNumber <<endl;
                ofstream fout("BadKe.txt");
                char BadElString[50];
                sprintf(BadElString,"Bad Ke - Element # %d",elementNumber);
                Ke->print(BadElString,&fout);
                fout.close();
				exit(1);}
            }
}
//======================================================================
void BasicElement::extractSolution( double * uGlobal, double ** elementDisp)
{int i;

   if(uGlobal !=NULL) 
     {
      int firstDof;
      int numberOfDof;
      for(i=0 ; i<numNodesPerElement;i++) {
	      firstDof = node[i].getFirstDof();
         numberOfDof = node[i].getNumDof();              
         switch (numberOfDof){
            case 1:
		{elementDisp[0][i] = uGlobal[firstDof];   break;}
            case 2:
		{elementDisp[0][i] = uGlobal[firstDof];
		 elementDisp[1][i] = uGlobal[++firstDof]; break;}
            case 3:
		{elementDisp[0][i] = uGlobal[firstDof];
		 elementDisp[1][i] = uGlobal[++firstDof];
       elementDisp[2][i] = uGlobal[++firstDof];
			 break;}
          case 4:
		{elementDisp[0][i] = uGlobal[firstDof];
		 elementDisp[1][i] = uGlobal[++firstDof];
       elementDisp[2][i] = uGlobal[++firstDof];
elementDisp[3][i] = uGlobal[++firstDof];
		//BETA_OUT<<elementDisp[0][i]<<"  "
		//	 <<elementDisp[1][i]<<"  "
       //   <<elementDisp[2][i]<<endl;
			 break;}
            default:
                {
					BETA_OUT<<"Case not implemented yet"<<endl;
					BETA_OUT<<"Check node #"<<node[i].nodeNum <<". Incorrect NumDofs." << endl;
					exit(1);
                }
            }//end of switch
	  }
     }
}
//=========================================================================
void BasicElement::outputResults()
{
 BETA_OUT<<"Virtual function ....outputResults"<<endl;
}
//=========================================================================
void BasicElement::print(ostream &)
{
 //BETA_OUT<<"connectivity"<<endl;
 for(int i=0; i<numNodesPerElement; i++)
	 {
	  BETA_OUT<<node[i].nodeNum<<" ";
	 }
	 BETA_OUT<<endl;
}
//=====================================================================
void BasicElement::setConnectivity(int nnpe, int * conn, NodeGroup *nodeList)
{
	numNodesPerElement=nnpe;
	if( !node.empty() ) node.clear();
	node.reserve(numNodesPerElement);
	char arrayname[256];
	sprintf(arrayname,"Element %d: Connectivity List",elementNumber);
	node.SetName(arrayname);

	int nodeOffset=0;
	if( (*nodeList)[0].nodeNum==1) nodeOffset=1;

	for(int i=0;i<numNodesPerElement;i++) {
		int nodeNum=conn[i];
		node.push_back( (*nodeList)[nodeNum-nodeOffset] );
		if(i==0 || i%2 == 0) node[i].isMidNode=false; else node[i].isMidNode=true; //JV010803 for interfacing with GeomPack++
	}

}
//=====================================================================
void BasicElement::read(istream &istrm, NodeGroup *nodeList, int nnpe)
{
	int i,nodeNum;
	int nodeOffset=0;
	if( (*nodeList)[0].nodeNum==1) nodeOffset=1;

	numNodesPerElement=nnpe;
	if( !node.empty() ) node.clear();
	node.reserve(numNodesPerElement);

	char arrayname[256];
	sprintf(arrayname,"Element %d: Connectivity List",elementNumber);
	node.SetName(arrayname);

	for(i=0;i<numNodesPerElement;i++) {
		istrm >> nodeNum;
		node.push_back( (*nodeList)[nodeNum-nodeOffset] );
	}
}
//=====================================================================
void BasicElement::readBinary(FILE *fp, NodeGroup *nodeList, int nnpe)
{
	const int SIZEOF_INT	=sizeof(int);

	int i,nodeNum;
	int nodeOffset=0;
	if( (*nodeList)[0].nodeNum==1) nodeOffset=1;

	numNodesPerElement=nnpe;
	if( !node.empty() ) node.clear();
	node.reserve(numNodesPerElement);

	char arrayname[256];
	sprintf(arrayname,"Element %d: Connectivity List",elementNumber);
	node.SetName(arrayname);

	for(i=0;i<numNodesPerElement;i++) {
		fread(&nodeNum, SIZEOF_INT, 1, fp);
		node.push_back( (*nodeList)[nodeNum-nodeOffset] );
	}
}
//=====================================================================
void BasicElement::setMidNodes(int numDims)
//JV010803 initially written for interfacing with GeomPack++
{
	switch(abs(numDims)){
		case 1:
			//handle only linear and quadratic 1D elements
			if(numNodesPerElement==2){
				node[0].isMidNode=false; 
				node[1].isMidNode=false; 
			}
			else{
				node[0].isMidNode=false; 
				node[1].isMidNode=true; 
				node[2].isMidNode=false; 
			}
			break;
		case 2:
		case 3:
			//Adds possibility for linear hex mesh
			for(int i=0;i<numNodesPerElement;i++) {
				node[i].isMidNode=false;
				if(numNodesPerElement == 20){
					if(i!=0 || i%2 != 0) 
						node[i].isMidNode=true; 
				}
			}
			break;
		default:
			cout << "A " << numDims << " dimensional mesh?! What have you been smoking? :)" << endl;
			exit(1);
	}
}
//=====================================================================
void BasicElement::readAllElementConnectivity(istream &istrm,
															 NodeGroup *node,
															 ElementGroup *element,
															 int numElements)
{
 for(int i=0; i< numElements; i++)
	{
	     int number;
		 int nnpe;
	     istrm >> number >> nnpe;
         (*element)[number].setElementNumber(number);
	     (*element)[number].read(istrm,node, nnpe);
	}
}
//====================================================================
void  BasicElement::setAnalysisType(istream *inStream, ostream &outStream) 
{ 
	//figure out what to do with this function!!
	//maybe analysisType should not be static
	//you cannot write to anything other than a static ostream or cout !!!

	(*inStream) >> BasicElement::analysisType;
	 outStream <<"Analysis type: "; 
	 // Eventually change this... use enumeration
	if(analysisType==0) (outStream) << ">>LINEAR<<"<<endl;
	if(analysisType==1) (outStream) << ">>GEOMETRIC_NONLINEAR<<"<<endl;
}
//==========================================================
void BasicElement::extractNodalCoordinates(double *xCoor,
				                           double *yCoor,
										   double *zCoor)
{int i; 
  for(i=0;i<numNodesPerElement;i++) {
      xCoor[i] = node[i].x;
      yCoor[i] = node[i].y;
      zCoor[i] = node[i].z;
   }
}
//=========================================================================
void BasicElement::setPointers(ElementWorkspace * space)
{
	if(!requiresParallelESAConsideration){
	Fe  = space->Fe;
	initialFe = space->initialFe;
	Ke = space->Ke;
	dofList = space->dofList; }

	equations		= space->equations;
}
//==========================================================================
void BasicElement::setEquationsPointer(Equations *equationsPointer)
{
	equations = equationsPointer;
}
//==========================================================================
double * BasicElement::getForceReVector(int Nincre)
{
	double temp=0.;
	double *tempi=&temp;
	return tempi;
}
//==========================================================================
void BasicElement::allocateISVs()
{
    if(!material){
        std::cout<<"!!!ERROR!!!: BasicElement::allocateISVs" << std::endl;
        std::cout<<"Cannot allocate ISVs if no material is assigned to element." << std::endl;
        exit(1);
    }
    material->allocateISVs(isv,1);
}
//==========================================================================
bool BasicElement::checkForAllocatedISVs()const
{
    if(!isv.empty()){
        return !(isv[0]==NULL);//Will be NULL if not allocated.
    } else {
        return false;
    }
}
//==========================================================================
void BasicElement::openOutputFiles(FileManager *fm, int flags)
{}
//=========================================================================
void BasicElement::closeOutputFiles(FileManager *fm, int flags)
{}
//=========================================================================
void BasicElement::SetActiveFlag(bool flag)
{
	int i;
	activeElementFlag=flag;
	for(i=0;i<numNodesPerElement;i++){
		node[i].isActive=flag;
	}
}
//=========================================================================
bool BasicElement::setNodalValues(double constantValue)
{
	if(nodalValues==0) return false;
	for(int i=0; i<numNodesPerElement; i++){
		nodalValues[i]=constantValue;
	}
	return true;
}
//=========================================================================
bool BasicElement::setNodalValues(vector<double> &list)
{
	if(nodalValues==0) return false;
	if(numNodesPerElement != list.size()) return false;
	for(int i=0; i<numNodesPerElement; i++){
		nodalValues[i]=list[i];
	}
	return true;
}

//=========================================================================
void BasicElement::AssignElementsToThreads(Array<threadWorkspace> &tList)
{
	//Required for parallel assembly
	if(equations==0) {
		cout << "equations pointer not set!" <<endl;
		cout << "must be set before calling specifyNonZeroLocation()!" <<endl;
		
		exit(1);
	}

	int *map,elemRank;

	map = calculateDofList();
	elemRank = getNumDof();
	equations->assignElementsToThreads(map,elemRank,tList,elementNumber); //bco use numthreads object....need to give element access to thread workspace or find work around
}
//========================================================================
void BasicElement::allocateElementalDataMemory(int WorkspaceMaxNumberDof)
{
	//Required for parallel assembly (ESA)
	Ke = new Matrix(WorkspaceMaxNumberDof,WorkspaceMaxNumberDof);
	Fe = new double [WorkspaceMaxNumberDof];
	initialFe = new double [WorkspaceMaxNumberDof];
	dofList = new int [WorkspaceMaxNumberDof];
}

//========================================================================
void BasicElement::releaseElementalDataMemory()
{

	if(Ke) delete Ke;						Ke = NULL;
	if(Fe) delete [] Fe;					Fe = NULL;
	if(initialFe) delete [] initialFe;		initialFe = NULL;
	if(dofList) delete [] dofList;			dofList = NULL;
}
//========================================================================
void BasicElement::getMasterNodeCoords(double* const MasterCoordsArray)
{
    //Assumes that MasterCoordsArray is appropriately sized for the current element
    //Fills in the array and sets the ElementOrder value

    // To explain the constness - you can change the data that MasterCoordsArray is pointing to,
    // but you can't modify the pointer itself (can't change where it is pointing)

    double* MasterCoordsPointer;
    //Right now the master coordinates have to be hard-coded in :-(
    //1D elements
    double E1D2N[2] = {-1.0,
                        1.0};

    double E1D3N[3] = {-1.0,
                        0.0,
                        1.0};
    //2D elements
    double E2D3N[6] = {0.0,0.0,
                       1.0,0.0,
                       0.0,1.0};

    double E2D6N[12] = {0.0,0.0,
                        0.5,0.0,
                        1.0,0.0,
                        0.5,0.5,
                        0.0,1.0,
                        0.0,0.5};

    double E2D4N[8] = {-1.0,-1.0,
                        1.0,-1.0,
                        1.0,1.0,
                       -1.0,1.0}; //Check

    double E2D8N[16] = {-1.0,-1.0,
                         0.0,-1.0,
                         1.0,-1.0,
                         1.0,0.0,
                         1.0,1.0,
                         0.0,1.0,
                         -1.0,1.0,
                         -1.0,0.0}; //Check
    //3D elments
    double E3D4N[12]   =  {0, 0, 0,
                           1, 0, 0,
                           0, 1, 0,
                           0, 0, 1};

    double E3D10N[30]  = {0,    0,    0,
                          1,    0,    0,
                          0,    1,    0,
                          0,    0,    1,
                          0.5,  0,    0,
                          0.5,  0.5,  0,
                          0,    0.5,  0,
                          0,    0,    0.5,
                          0.5,  0,    0.5,
                          0,    0.5,  0.5};
    
    double E3D6N[18]  = {0,  0, -1,
                         1,  0, -1,
                         0,  1, -1,
                         0,  0,  1,
                         1,  0,  1,
                         0,  1,  1};
    
    
    double E3D15N[45]  = {0,    0,   -1,
                          0.5,  0,   -1,
                          1,    0,   -1,
                          0.5,  0.5, -1,
                          0,    1,   -1,
                          0,    0.5, -1,
                          0,    0,    0,
                          1,    0,    0,
                          0,    1,    0,
                          0,    0,    1,
                          0.5,  0,    1,
                          1,    0,    1,
                          0.5,  0.5,  1,
                          0,    1,    1,
                          0,    0.5,  1};

    double E3D8N[24]   =  {-1,-1,-1,
                           -1, 1,-1,
                           -1, 1, 1,
                           -1,-1, 1,
                            1,-1,-1,
                            1, 1,-1,
                            1, 1, 1,
                            1,-1, 1};

    double E3D20N[60] = {-1,-1,-1, 
                         -1, 0,-1,
                         -1, 1,-1,
                         -1, 1, 0,
                         -1, 1, 1,
                         -1, 0, 1,
                         -1,-1, 1,
                         -1,-1, 0,
                          0,-1,-1,
                          0, 1,-1,
                          0, 1, 1,
                          0,-1, 1,
                          1,-1,-1,
                          1, 0,-1,
                          1, 1,-1,
                          1, 1, 0,
                          1, 1, 1,
                          1, 0, 1,
                          1,-1, 1,
                          1,-1, 0};    
    #define n -0.3333333333333
    #define p  0.3333333333333
    double E3D32N[96] =   {-1,-1,-1,
                           -1, n,-1,
                           -1, p,-1,
                           -1, 1,-1,
                           -1, 1, n,
                           -1, 1, p,
                           -1, 1, 1,
                           -1, p, 1,
                           -1, n, 1,
                           -1,-1, 1,
                           -1,-1, p,
                           -1,-1, n,
                            n,-1,-1,
                            n, 1,-1,
                            n, 1, 1,
                            n,-1, 1,
                            p,-1,-1,
                            p, 1,-1,
                            p, 1, 1,
                            p,-1, 1,
                            1,-1,-1,
                            1, n,-1,
                            1, p,-1,
                            1, 1,-1,
                            1, 1, n,
                            1, 1, p,
                            1, 1, 1,
                            1, p, 1,
                            1, n, 1,
                            1,-1, 1,
                            1,-1, p,
                            1,-1, n};
    #undef n
    #undef p
    int NumDims = 0;
    switch(elementType){
        //1D
        case L1D2N:// 1-D Linear 2-node linear element
            MasterCoordsPointer = E1D2N;
            NumDims = 1;
            break;
        case L1D3N:// 1-D Linear 3-node linear element
            MasterCoordsPointer = E1D3N;
            NumDims = 1;
            break;
        //2D
        case Q2D4N:// quadrilateral 4-node linear element
            MasterCoordsPointer = E2D4N;
            NumDims = 2;
            break;
        case Q2D8N:// quadrilateral 8-node quadratic element
            MasterCoordsPointer = E2D8N;
            NumDims = 2;
            break;
        case T2D3N:// triangular 2D 3-node element (linear)
            MasterCoordsPointer = E2D3N;
            NumDims = 2;
            break;
        case T2D6N:// triangular 2D 6-node element (quadratic)
            MasterCoordsPointer = E2D6N;
            NumDims = 2;
            break;
        //3D
        case T3D4N:// 3D tetrahedral 4-node  (linear) element
            MasterCoordsPointer = E3D4N;
            NumDims = 3;
            break;
        case T3D10N:// 3D tetrahedral 10-node  (quadratic) element
            MasterCoordsPointer = E3D10N;
            NumDims = 3;
            break;
        case W3D6N:// 3D wedge 6-node  (linear) element
            MasterCoordsPointer = E3D6N;
            NumDims = 3;
            break;
        case W3D15N:// 3D wedge 15-node  (quadratic) element
            MasterCoordsPointer = E3D15N;
            NumDims = 3;
            break;
        case H3D8N:// 3D hexahedral 8-node (linear) element
            MasterCoordsPointer = E3D8N;
            NumDims = 3;
            break;
        case H3D20N:// 3D hexahedral 20-node (quadratic) element
            MasterCoordsPointer = E3D20N;
            NumDims = 3;
            break;
        default:
            cout << "No master nodes defined for Element Type " << elementType << endl;
            exit(1);
            break;
    }
    memcpy(MasterCoordsArray,MasterCoordsPointer,sizeof(double)*numNodesPerElement*NumDims);
    //for(int i=0;i<numNodesPerElement;++i){
    //    for(int j=0;j<NumDims;++j){
    //        MasterCoordsArray[i][j] = MasterCoordsPointer[i*NumDims + j];
    //    }
    //}
}

int BasicElement::getElementOrder()
{
    int ElementOrder = 0;
    switch(elementType){
        //1D
        case L1D2N:// 1-D Linear 2-node linear element
            ElementOrder = 1;
            break;
        case L1D3N:// 1-D Linear 3-node linear element
            ElementOrder = 2;
            break;
        //2D
        case Q2D4N:// quadrilateral 4-node linear element
            ElementOrder = 1;
            break;
        case Q2D8N:// quadrilateral 8-node quadratic element
            ElementOrder = 2;
            break;
        case T2D3N:// triangular 2D 3-node element (linear)
            ElementOrder = 1;
            break;
        case T2D6N:// triangular 2D 6-node element (quadratic)
            ElementOrder = 2;
            break;
        //3D
        case T3D4N:// 3D tetrahedral 4-node  (linear) element
            ElementOrder = 1;
            break;
        case T3D10N:// 3D tetrahedral 10-node  (quadratic) element
            ElementOrder = 2;
            break;
        case W3D6N:// 3D wedge 6-node  (linear) element
            ElementOrder = 1;
            break;
        case W3D15N:// 3D wedge 15-node  (quadratic) element
            ElementOrder = 2;
            break;
        case H3D8N:// 3D hexahedral 8-node (linear) element
            ElementOrder = 1;
            break;
        case H3D20N:// 3D hexahedral 20-node (quadratic) element
            ElementOrder = 2;
            break;
        default:
            cout << "Order not defined for element type " << elementType << " (El # " << elementNumber << ")" << endl;
            exit(1);
            break;
    }
    return ElementOrder;
}
//=========================================================================
//                             Static Data                        
//=========================================================================
int BasicElement::analysisType=0;
