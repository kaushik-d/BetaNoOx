#include "stdafx.h"
#include "utility/utility.h"

#include "isoElement.hpp"
#include "InterpAndSolDerivatives.hpp"
#include "materials/Material.hpp"

#include "math/gaussQuadrature/hexQuadPoint.hpp"
#include "math/gaussQuadrature/tetQuadPoint.hpp"
#include "math/gaussQuadrature/wedgeQuadPoint.hpp"
#include "math/gaussQuadrature/lineQuadPoint.hpp"
#include "math/gaussQuadrature/triQuadPoint.hpp"
#include "math/gaussQuadrature/quadQuadPoint.hpp"

//===========
void deriv1D( int numberOfInterp, double zeta, double * p_S_Lx1);
void shape1D( int numberOfInterp, double zeta, double * S  );
//===========
void deriv2D(int numberOfInterp,
              double xi,  double eta, 				  
				  double * p_S_Lx1, double * p_S_Lx2);

void  shape2D(int numberOfInterp,
              double xi,  double eta,				  
				  double * S);

//============
void deriv3D(int numberOfInterp,
              double xi,  double eta, double zeta,				  
				  double * p_S_Lx1, double * p_S_Lx2, double * p_S_Lx3);

void  shape3D(int numberOfInterp,
              double xi,  double eta, double zeta,				  
				  double * S);
//============
extern	int			verboseFlag;
//===================================================================
int IsoElement::getTotalNumberOfIPs() const
{
switch(numberOfIntegrationDirections){
case 1:
	return integrationOrder ;
	break;	
case 2:
    switch(numNodesPerElement){
    case 4:
    case 8:
	    return integrationOrder * integrationOrder;
        break;
    case 3:
    case 6:
        switch(integrationOrder){
        case 1:
            return 1;
            break;
        case 2:
            return 3;
            break;
        case 3:
            return 4;
            break;
        case 4:
            return 7;
            break;
        default:
            banner("Specified Integration order not defined for triangles", BETA_OUT); exit(1);
        }
	break;
    }
case 3:
	switch(numNodesPerElement){
	case 8:
	case 20:
		//Quadrilateral elements (tensor product)
		return integrationOrder * integrationOrder * integrationOrder;
		break;
	case 4:
	case 10:
		//Tetrahedral Elements
		switch(integrationOrder){
		case 1:
			return 1;
			break;
		case 2:
			return 4;
			break;
		case 3:
			return 15;
			break;
		default:
			banner("Specified Integration order not defined for tetrahedrons", BETA_OUT); exit(1);
		}//end of switch over tetrahedron elements
		break;
	case 6:
	case 15:
		//Wedge Elements
		switch(integrationOrder){
		case 1:
			return 2;
			break;
		case 2:
			return 9;
			break;
		case 3:
			return 18;
			break;
		default:
			banner("Specified Integration order not defined for wedges", BETA_OUT); exit(1);
		}//end of switch over wedge elements
	default:
		banner("Integration poitns not defined for given element type", BETA_OUT); exit(1);
	}//end over switch for 3D element types
	break;
	

default:
	banner("Case not defined", BETA_OUT); exit(1);
}//end Of Switch
}//end of getTotalNumberOfIPs()
//===================================================================
void IsoElement::setIntegrationOrder(int i)
{
    if(checkForAllocatedISVs()){
        std::cout<<"!!!ERROR!!!- IsoElement::setIntegrationOrder:" <<std::endl;
        std::cout<<"Cannot change integration order after ISVs are allocated!" << std::endl;
        exit(1);
    }
    integrationOrder = i;
}
//===================================================================
void IsoElement::getQuadraturePoints(GaussPointList &gpL)
{

switch(numberOfIntegrationDirections){

case 1:
		LineQuadPoints::getQuadraturePointsForLine( 
										integrationOrder,&(gpL.xiList),
										&(gpL.weightList)  );                              
		gpL.totalNumIPs = integrationOrder  ;
break;

case 2:
    switch(numNodesPerElement){
    case 4:
    case 8:
        QuadQuadPoints::getQuadraturePointsForQuad( 
										integrationOrder,&(gpL.xiList),
										&(gpL.etaList),
										&(gpL.weightList)  );

		gpL.totalNumIPs = integrationOrder * integrationOrder ;
        break;
    case 3:
    case 6:
        TriQuadPoints::getQuadraturePointsForTri( 
										integrationOrder,&(gpL.xiList),
										&(gpL.etaList),
										&(gpL.weightList)  );

		gpL.totalNumIPs = getTotalNumberOfIPs();
        break;
    }

break;

case 3:
	switch(numNodesPerElement){
	case 8:
	case 20:
		//Hexahedron elements
		HexQuadPoints::getQuadraturePointsForHex( 
										integrationOrder,
										&(gpL.xiList),  &(gpL.etaList),
										&(gpL.zetaList),&(gpL.weightList)  );

		gpL.totalNumIPs = integrationOrder * integrationOrder * integrationOrder;
		break;
	case 4:
	case 10:
		//Tetrahedron elements
		TetQuadPoints::getQuadraturePointsForTet( 
										integrationOrder,
										&(gpL.xiList),  &(gpL.etaList),
										&(gpL.zetaList),&(gpL.weightList)  );

		gpL.totalNumIPs = getTotalNumberOfIPs();
		break;
	case 6:
	case 15:
		//Wedge elements
		WedgeQuadPoints::getQuadraturePointsForWedge( 
										integrationOrder,
										&(gpL.xiList),  &(gpL.etaList),
										&(gpL.zetaList),&(gpL.weightList)  );

		gpL.totalNumIPs = getTotalNumberOfIPs();
		break;
	default:
		banner("Unknown Element Type", BETA_OUT); exit(1);
	}//end of switch over element type
		
break;

default:
	banner("Case not defined", BETA_OUT); exit(1);
}//end Of Switch
}//end of getQuadraturePoints()

//====================================================================
void IsoElement::printQuadPoints(int totalNumIPs, ostream &SO)
{
    //Outputs Quad Point locations and associated volume.
    double x,y,z;

	GaussPointList	   gaussPointList;
    InterpAndSolDerivatives interpData(numNodesPerElement,numNodesPerElement); interpData.allocate();

    interpData.numberOfInterp = numNodesPerElement;
    extractNodalCoordinates(interpData.xCoor,interpData.yCoor,interpData.zCoor);

    getQuadraturePoints(gaussPointList); 

	double weight,vol;
	int ip,i;

	for(ip=0;ip<totalNumIPs;ip++){
		interpolation(ip, gaussPointList, &interpData);
		jacobian(&interpData);

		x=y=z=0.0;
		for (i=0; i<numNodesPerElement; i++){
			x=x+interpData.S[i]*interpData.xCoor[i];
			y=y+interpData.S[i]*interpData.yCoor[i];
			z=z+interpData.S[i]*interpData.zCoor[i];
		}
		weight = gaussPointList.weightList[ip];
		vol = interpData.detJ * weight * getIntegrationFactor();

        SO << DOUBLE_FORMAT << x << "\t" << y << "\t" << z << "\t" << vol << endl;
    }
}

//=============================================================
void IsoElement::allocateISVs()
{
    material->allocateISVs(isv,getTotalNumberOfIPs());
}
//=============================================================
void IsoElement::jacobian(InterpAndSolDerivatives *id_ptr) 
{
	InterpAndSolDerivatives &id=*id_ptr;

    double px1,px2, px3, py1,py2,py3, pz1,pz2,pz3;
	 double dett1,dett2,dett3,dett4,dett5,dett6,dett7,dett8,dett9;
	 double invDetJ;
	 double aji11, aji12, aji13,  aji21, aji22, aji23,
		      aji31, aji32, aji33;
    int i,j;
double p_x_xi;
int numberOfInterp = id.numberOfInterp;

double *S=id.S, *p_S_Lx1=id.p_S_Lx1, *p_S_Lx2=id.p_S_Lx2,*p_S_Lx3=id.p_S_Lx3,
		        *p_S_x1=id.p_S_x1,    *p_S_x2=id.p_S_x2,     *p_S_x3=id.p_S_x3;


double *xCoor=id.xCoor, *yCoor=id.yCoor, *zCoor=id.zCoor;

//====================================
switch(numberOfIntegrationDirections){

case 1:
double invDetj;
p_x_xi = 0.;

for(i=0; i< numberOfInterp; i++){p_x_xi += p_S_Lx1[i] * xCoor[i];}
// Determinant of J and its inverse
             id.detJ = p_x_xi;
             invDetj = 1/id.detJ;

// Global derivatives of shape functions
 for(i=0; i< numberOfInterp; i++){p_S_x1[i] = 0.;}
 for(i=0; i< numberOfInterp; i++){p_S_x1[i] = invDetj * p_S_Lx1[i] ;}	
break;

case 2:

double J[2][2], invJ[2][2];
              J[0][0] = 0.;
              J[0][1] = 0.;
              J[1][0] = 0.;
              J[1][1] = 0.;
double temp;

 for(i=0; i< numberOfInterp; i++)
     {
      J[0][0] += p_S_Lx1[i] * xCoor[i];
      J[0][1] += p_S_Lx1[i] * yCoor[i];
      J[1][0] += p_S_Lx2[i] * xCoor[i];
      J[1][1] += p_S_Lx2[i] * yCoor[i];
     }
// Determinant of J and its inverse
             id.detJ = J[0][0]*J[1][1]-
                    J[0][1]*J[1][0];
             temp = 1/id.detJ;
             invJ[0][0]= J[1][1] * temp;
             invJ[0][1]=-J[0][1] * temp;
             invJ[1][0]=-J[1][0] * temp;
             invJ[1][1]= J[0][0] * temp;
// Global derivatives of shape functions
 for(i=0; i< numberOfInterp; i++)
     {
      p_S_x1[i] = 0.;
      p_S_x2[i] = 0.;
     }

 for(i=0; i< numberOfInterp; i++)
     {
      p_S_x1[i] = invJ[0][0] * p_S_Lx1[i] + invJ[0][1] * p_S_Lx2[i];
      p_S_x2[i] = invJ[1][0] * p_S_Lx1[i] + invJ[1][1] * p_S_Lx2[i];
     }
break;

case 3:
	 px1=px2=px3=py1=py2=py3=pz1=pz2=pz3=0.;

	for(i=0; i< numberOfInterp; i++)
     {
         px1 += p_S_Lx1[i]   * xCoor[i];
	      px2 += p_S_Lx2[i]  * xCoor[i];	
	      px3 += p_S_Lx3[i] * xCoor[i];	
	      py1 += p_S_Lx1[i]   * yCoor[i];
	      py2 += p_S_Lx2[i]  * yCoor[i];
	      py3 += p_S_Lx3[i] * yCoor[i];
	      pz1 += p_S_Lx1[i]   * zCoor[i];
	      pz2 += p_S_Lx2[i]  * zCoor[i];
	      pz3 += p_S_Lx3[i] * zCoor[i];
	  }

		dett1=  py2 * pz3 - py3 * pz2 ; 
		dett2=  px2 * pz3 - px3 * pz2 ;
      dett3=  px2 * py3 - px3 * py2 ;
     	dett4=  py1 * pz3 - py3 * pz1 ;
     	dett5=  px1 * pz3 - px3 * pz1 ;
     	dett6=  px1 * py3 - px3 * py1 ;
     	dett7=  py1 * pz2 - py2 * pz1 ;
     	dett8=  px1 * pz2 - px2 * pz1 ;
     	dett9=  px1 * py2 - px2 * py1 ;

      id.detJ =  px1 * dett1 - py1 * dett2 + pz1 * dett3;

      invDetJ = 1./id.detJ;

// Form inverse of jacobian
           aji11  =   dett1 * invDetJ; 
           aji12  =  -dett4 * invDetJ; 
           aji13  =   dett7 * invDetJ; 
           aji21  =  -dett2 * invDetJ; 
           aji22  =   dett5 * invDetJ; 
           aji23  =  -dett8 * invDetJ;
           aji31  =   dett3 * invDetJ; 
           aji32  =  -dett6 * invDetJ; 
           aji33  =   dett9 * invDetJ; 

//...Derivatives of shape functions in term of global coordintes
   for(j=0; j< numberOfInterp; j++)
     {
      p_S_x1[j] = aji11 * p_S_Lx1[j] + aji12 * p_S_Lx2[j] + aji13 * p_S_Lx3[j];
      p_S_x2[j] = aji21 * p_S_Lx1[j] + aji22 * p_S_Lx2[j] + aji23 * p_S_Lx3[j];
      p_S_x3[j] = aji31 * p_S_Lx1[j] + aji32 * p_S_Lx2[j] + aji33 * p_S_Lx3[j];
	   }
break;

default:
	banner("Case not defined", BETA_OUT); exit(1);
}//end Of Switch

//JV120205
//calculating global coordinates of current quadrature point
id.x=0;
id.y=0;
id.z=0;
 for(i=0; i< numberOfInterp; i++){
	 id.x += S[i] * xCoor[i];
	 id.y += S[i] * yCoor[i];
	 id.z += S[i] * zCoor[i];
 }

}  //end of Jacobian

//==============================================================
void IsoElement::interpolation(int integPointNumber,
											  GaussPointList &gpList,
											  InterpAndSolDerivatives *id_ptr)
{
	InterpAndSolDerivatives &id=*id_ptr;

switch(numberOfIntegrationDirections){

case 1:
		id.numberOfInterp = numNodesPerElement;

		deriv1D(id.numberOfInterp, gpList.xiList[integPointNumber],id.p_S_Lx1);
		shape1D(id.numberOfInterp, gpList.xiList[integPointNumber],id.S  );
break;

case 2:
		id.numberOfInterp = numNodesPerElement;

		deriv2D(id.numberOfInterp, gpList.xiList[integPointNumber],
											  gpList.etaList[integPointNumber],
											  id.p_S_Lx1, id.p_S_Lx2);
		shape2D(id.numberOfInterp, gpList.xiList[integPointNumber],
											  gpList.etaList[integPointNumber], id.S);
break;

case 3:
		id.numberOfInterp = numNodesPerElement;


		deriv3D(id.numberOfInterp, gpList.xiList[integPointNumber],
											 gpList.etaList[integPointNumber],
											 gpList.zetaList[integPointNumber],
											 id.p_S_Lx1, id.p_S_Lx2, id.p_S_Lx3);
		shape3D(id.numberOfInterp, gpList.xiList[integPointNumber],
											 gpList.etaList[integPointNumber],
											 gpList.zetaList[integPointNumber], id.S);
break;

default:
	banner("Case not defined", BETA_OUT); exit(1);
}//end Of Switch
}
//====================================================================
IsoElement::IsoElement(void):
BasicElement()
{ 
    extrapMatrix = NULL;
	numberOfIntegrationDirections=0;
	setIntegrationOrder(3);

    ElementNodeRotations = 0;
    OrientationAtIP = 0;
    OrientationAtIP = new Rotation();
}
//====================================================================
IsoElement::~IsoElement(void)
{ 
    if (ElementNodeRotations)
        delete [] ElementNodeRotations;
    if(OrientationAtIP)
        delete OrientationAtIP;
	//if(traceFlag==1){cout<< "Deleting IsoElement element"<<endl;}
}
//======================================================================
void IsoElement::readSpecialCommand(istream &inStream, ElementGroup *element,
								int numElements, char * command)
{
	int i;
	BETA_OUT<<"IsoElement::readSpecialCommand"<<endl;

//================================================================
	if (COMPARE(command, "SetIntegrationOrder")==0 ) {
		int first, last, increment, quadOrder;
		IsoElement *p_elem;

		inStream>>first>>last>>increment;
		inStream>>quadOrder;

		BETA_OUT << "\tQuadrature order = " << quadOrder << " for elements "
		<< first <<" to "<<last << " by " << increment << '\n';

		for(i=first;i<numElements;i +=increment) { 
			p_elem = ((IsoElement*) &(*element)[i]);
			p_elem->setIntegrationOrder(quadOrder);
		}// i
		return;
	}//SetIntegrationOrder
//================================================================
    if (COMPARE(command, "SetElementMaterialAngles")==0 ) {
        //One rotation axis, but rotation can be different at each node
        IsoElement * ThisElement;
        int elementNumber, RotationAxis;
        double * Angles = 0;
        bool AxisWarningOutput = false;
        bool MaterialDefinitionWarningOutput = false;
        for(i=0;i<numElements;++i){
            ThisElement = ((IsoElement*) &(*element)[i]);
            inStream >> elementNumber >> RotationAxis;
            if(ThisElement->elementNumber!=i) {
		        BETA_OUT<<" i, elemNumber "<<i<<"  "<<elementNumber<<endl;
		        banner("Incorrect element order...", BETA_OUT); exit(1);
		    }

            ThisElement = (IsoElement*) &(*element)[i];
            Angles = new double[ThisElement->numNodesPerElement];
            for(int j=0;j<ThisElement->numNodesPerElement;++j){
                inStream >> Angles[j];
            }

            if(RotationAxis == 1){
                // If the material axis is something other than x2, there are 2
                // likely reasons.  Treat them in the following way:
                // 1 - No material rotation:
                //   If this is the case, then the material has probably been defined
                //   with the 2 direction as the fiber direction.  Go ahead and rotate
                //   about the 1 axis, but print a warning as this way of defining
                //   material properties leads to confusion and needs to end
                // 2 - With material rotation:
                //   The rotations are defined to accomodate the former order of
                //   rotations, first the material rotation, then the mangles rotation.
                //   The new approach does mangles, then material (to restore the ability
                //   to handle braids).  Switch the axis for mangles rotation to 1 and
                //   use the negative angle, and print a warning because in the future
                //   rotations need to be defined differently.
                if(ThisElement->material->getMaterialRotation().no_rotation()){
                    if(MaterialDefinitionWarningOutput){
                        cout << "!!!WARNING!!! - It appears that you are rotating about the" << endl;
                        cout << "x1 axis because the material is defined so that the x2 axis is" << endl;
                        cout << "the fiber direction.  In the future, define the material so" << endl;
                        cout << "that x1 is the fiber direction and apply a rotation using the" << endl;
                        cout << "readAngles command, and apply a mangles rotation about x2" << endl;
                    }
                }else{
                    if(!AxisWarningOutput){
                    
                        cout << "!!!WARNING!!! - Mangles should always be about the x2 axis." << endl;
                        cout << "Mangles rotations about the 1 axis are being adjusted so" << endl;
                        cout << "that they are about the x2 axis. If your material is defined" << endl;
                        cout << "with the local 1  axis as the fiber direction with a rotation" << endl;
                        cout << "about x3, then things will work okay.  If the local x2 axis" << endl;
                        cout << "is the material's fiber direction then your model is broken." << endl;
                        AxisWarningOutput = true;
                    }
                    // If the old mangles rotation was about the 1 axis, rotating it 90 degrees about z 
                    // results in a rotation about the negative 2 axis (or a negative rotation about the 
                    // positive 2 axis)
                    RotationAxis = 2;
                    for(int j=0;j<ThisElement->numNodesPerElement;++j){
                        Angles[j] = -Angles[j];
                    }
                }
            }
            ThisElement->readRotationsFromMangles(RotationAxis,Angles);
            delete [] Angles;
            Angles = 0;
        }// End loop over elements to assgn mangles

        if(verboseFlag>=Min){
		    ofstream* matangles=filemanager->OpenOutputStream("MatAnglesOutput");
		    (*matangles) << "Coordinate System rotation data for each element" << endl;

		    for(i=0;i<numElements;i++) { 
			    (*matangles) << "Element: " << i << endl;
                ThisElement = (IsoElement*) &(*element)[i];
                for(int j=0;j<ThisElement->numNodesPerElement;++j){
                    (*matangles) << "\tNode " << ThisElement->node[j].nodeNum << "\t";
                    ThisElement->ElementNodeRotations[j].output_rotation(matangles);
                    (*matangles) << endl;
                }
		    }// i
            filemanager->CloseOutputStream(matangles);
	    }else{
		    BETA_OUT<<"Set verboseFlag to at least Min for echo of material angles info in file: MatAnglesOutput"<<endl;
	    }
		return;
	}//SetElementMaterialAngles
	else{// If no match, try super class
		EXIT_BETA("IsoElement::readSpecialCommand : Cannot recognize this command - " << command);
	}

};//readSpecialCommand

//======================================================================
void IsoElement::readRotationsFromMangles(const int &Axis, const double* Angles){
    // If the material has already been assigned, then ElementNodeRotations has already
    // been allocated, and the rotation for the material has already been applied.
    // Apply the mangles rotation before the material rotaiton.  If the material has not
    // yet been assigned, then give an error.  You can't simply apply the material rotation
    // when a material is assigned because materials get assigned to elements multiple
    // times in an analysis
    
    if(material == 0){ //Material must be assigned to the material before mangles is read
        cout << "!!!ERROR!!! IsoElement::readRotationsFromMangles-" << endl;
        cout << "Element material must be set before mangles is read" << endl;
        cout << "to allow for rotations to be combined correctly." << endl;
        exit(1);
    }    
    ElementNodeRotations = new Rotation[numNodesPerElement];
    for(int i=0;i<numNodesPerElement;++i){
        // General Note about Rotation Class and Angles:
        // Angles are input into Beta in terms of rotating a material about
        // a coordinate system.  The rotations class is defined in terms of 
        // rotating a coordinate system about a material.  One is the opposite
        // of the other, so when constructing rotation objects, input the opposite
        // of the angle that is used to input the rotation into Beta.

        // Mangles, then Material
        ElementNodeRotations[i] = Rotation(Axis,-Angles[i]);
        ElementNodeRotations[i].add_after_self(material->getMaterialRotation());
    }    
}
//======================================================================
void IsoElement::AllocateRotations()
{
    if(ElementNodeRotations != 0){
        cout << "WARNING - IsoElement::AllocateRotations()" << endl;
        cout << "ElementNodeRotations already allocated.  Skipping..." << endl;
        return;
    }
    ElementNodeRotations = new Rotation[numNodesPerElement];
    for(int i=0;i<numNodesPerElement;++i){
        ElementNodeRotations[i] = material->getMaterialRotation();
    }
}
//======================================================================
void IsoElement::readElementOrientation(istream* &IS)
{
    int ElNum,NumRotations;
    vector<int> RotationAxes;
    *IS >> ElNum;
    if(ElNum != elementNumber){
        cout << "ERROR:  IsoElement::readElementOrientation" << endl;
        cout << "Incorrect element numbering" << endl;
        exit(1);
    }
    *IS >> NumRotations;
    RotationAxes.resize(NumRotations);
    for(int i=0;i<NumRotations;++i){
        *IS >> RotationAxes.at(i);
    }
    double RotationAngle;
    ElementNodeRotations = new Rotation[numNodesPerElement];
    for(int i=0;i<numNodesPerElement;++i){
        // The rotation stored at the nodes is a single rotation which rotates the local coordinate system
        // to the global coordinate system.  The input euler angles, however, represent a series of rotations
        // going from the global coordinate system to the local material coordinates at a node.  Therefore, 
        // to generate the single rotation which will be stored, the sign of each input rotation must be 
        // switched, and the order must be reversed.

        // Note, one could read in all the angles and then loop over them in reverse, using 
        // Rotation::add_after_self to generate the composite rotation.  The code below simply loops
        // through the input rotations in order and uses Rotation::add_before_self, which yields an
        // equivalent result.
        for(int j=0;j<NumRotations;++j){
            *IS >> RotationAngle;
            ElementNodeRotations[i].add_before_self(Rotation(RotationAxes.at(j),-RotationAngle));
        }
    }
}
//======================================================================
void IsoElement::printQuadOrientations(ostream &SO)
{
    // Outputs the principle vectors of the local coordinate system expressed in the global (primed) 
    // coordinate system for each quad point
    GaussPointList	   gaussPointList;
    InterpAndSolDerivatives interpData(numNodesPerElement,numNodesPerElement); interpData.allocate();

    interpData.numberOfInterp = numNodesPerElement;
    extractNodalCoordinates(interpData.xCoor,interpData.yCoor,interpData.zCoor);

    getQuadraturePoints(gaussPointList);

    SO << elementNumber << " " << gaussPointList.totalNumIPs << endl;

    for(int ip=0;ip<gaussPointList.totalNumIPs;++ip){
        interpolation(ip, gaussPointList, &interpData);
        //OrientationAtIP->interpolate_rotation(&interpData,ElementNodeRotations);
        OrientationAtIP->interpolate_rotation_by_angle(&interpData,ElementNodeRotations);
        //OrientationAtIP->interpolate_old_beta(&interpData,ElementNodeRotations,material->getMaterialRotation());
        Matrix DCM(3,3);
        OrientationAtIP->get_direction_cosine_matrix(DCM);
        for(int i=0;i<3;++i){
            SO << "   ";
            for(int j=0;j<3;++j){
                SO << " " << DCM(j,i);
            }
        }
        SO << endl;
    }
}
//======================================================================
double IsoElement::CalculateVolume()
{
	double vol=0.0;
	GaussPointList gaussPointList;
	InterpAndSolDerivatives       interpData(numNodesPerElement,numNodesPerElement); interpData.allocate();
	getQuadraturePoints(gaussPointList);
	int totalNumIPs = gaussPointList.totalNumIPs;
	double integrationFactor, weight;

	double *elemDisp[3];
	elemDisp[0] = interpData.nodalValues_u1; 
	elemDisp[1] = interpData.nodalValues_u2; 
	elemDisp[2] = interpData.nodalValues_u3; 

	interpData.numberOfInterp = numNodesPerElement;
	extractNodalCoordinates(interpData.xCoor,interpData.yCoor,interpData.zCoor);

	for(int ip=0;ip<totalNumIPs;ip++){
		interpolation(ip, gaussPointList, &interpData);
		jacobian(&interpData);
		weight = gaussPointList.weightList[ip];
		integrationFactor = interpData.detJ * weight * getIntegrationFactor();
		vol += integrationFactor;
	}
	return vol;
}
//======================================================================
void IsoElement::calculateElementVolume(IntegrationDataType &idt)
{
	double &elementVolume           = bag->elementVolume;
	elementVolume += idt.integrationFactor;  // this can also be area
}
//======================================================================
void IsoElement::getExtrapolationMatrix(int LessOrder)
{
    //Creates an extrapolation matrix in the workspace if an appropriate one doesn't exist and
    //points extrapMatrix to it.
    int ElementOrder = getElementOrder();
    if(extrapMatrix == NULL || extrapMatrix->getDesiredOrder() != (ElementOrder-LessOrder)){
        GaussPointList GaussPoints;
        getQuadraturePoints(GaussPoints);
        for(int i=0;i<bag->extrapMatrices.size();++i){
            ExtrapolationMatrix* Mat = (ExtrapolationMatrix*) bag->extrapMatrices[i];
            if(Mat->getDimensionality() == numberOfIntegrationDirections &&
               Mat->getNumKnownCoords() == GaussPoints.totalNumIPs && 
               Mat->getNumUnknownCoords() == numNodesPerElement &&
               Mat->getDesiredOrder() == (ElementOrder - LessOrder) ){
                extrapMatrix = Mat;
                break;
            }
        }
        //If no appropriate matrix was found, then create one
        if(extrapMatrix == NULL){
            bag->extrapMatrices.push_back(NULL);
            double *MasterCoords = new double[numNodesPerElement*numberOfIntegrationDirections];
            getMasterNodeCoords(MasterCoords);
            bag->extrapMatrices.back() = new ExtrapolationMatrix(numNodesPerElement,
                                                                 numberOfIntegrationDirections,
                                                                 GaussPoints,
                                                                 MasterCoords,
                                                                 ElementOrder-LessOrder,
                                                                 TensorProduct);
            extrapMatrix = (ExtrapolationMatrix*) bag->extrapMatrices.back();
            delete [] MasterCoords;
        }
    }
}
//=========================================================================
void IsoElement::extrapolateQuadValuesToNodes(const double* const* valsAtQuadPoints, const int &numVals, double** &valsAtNodes)
{
    //Assumes that valsAtNodes has already been sized appropriately to hold the extrapolated values
    //(it should be #nodes x #vals)
    //valsAtQuadPoints should be a #quad points x #vals array
    getExtrapolationMatrix(0);
    int NumIPs = getTotalNumberOfIPs();
    for(int i=0;i<numNodesPerElement;++i){
        for(int j=0;j<numVals;++j){
            valsAtNodes[i][j] = 0.0;
            for(int k=0; k<NumIPs; ++k){
                valsAtNodes[i][j] += (*extrapMatrix)[i][k] * valsAtQuadPoints[k][j];
            }
        }
    }
}
//=========================================================================
void IsoElement::extrapolateQuadValuesToNodes(const double* valsAtQuadPoints, const int &numVals, double* &valsAtNodes)
{
    //For single value extrapolation
    //Assumes that valsAtNodes has already been sized appropriately to hold the extrapolated values
    //(it should be #nodes long)
    //valsAtQuadPoints should be a #quad points 
    getExtrapolationMatrix(0);
    int NumIPs = getTotalNumberOfIPs();
    for(int i=0;i<numNodesPerElement;++i){
        valsAtNodes[i] = 0.0;
        for(int j=0; j<NumIPs; ++j){
            valsAtNodes[i] += (*extrapMatrix)[i][j] * valsAtQuadPoints[j];
        }
    }
}
//=========================================================================
void getDoubleMinMax(double *x,int length, double &min, double &max);
//=========================================================================
bool IsoElement::SearchForPointInsideElement(double *pointCoord, double *localCoord)
{
	//See beta_docs/searchInsideElementForPoint.docx for documentation
	bool isPointInside = false;
	int numNodes = getNumNodesPerElement();

	double tol = 1.0e-5;
	double NRtol = 1.0e-8;
	double xi=0.0; double eta = 0.0; double zeta = 0.0;

	double *x,*y,*z = 0;
	x = new double [numNodes];
	y = new double [numNodes];
	z = new double [numNodes];
	double xp = pointCoord[0];
	double yp = pointCoord[1];
	double zp = pointCoord[2];

	extractNodalCoordinates(x,y,z);

	if(elementType == Q2D8N){ //for 8 Noded Quad (2D)
			
		//Filter out points that are far away from element
		bool insideBoundingBox = SearchElementBoundingBoxForPoint(pointCoord,tol);
		if(insideBoundingBox == true){
				isPointInside = SearchElementNodesForPoint(pointCoord,localCoord,tol);
				if(isPointInside == false){
					isPointInside = InverseInterpolationForPointQuad(pointCoord,localCoord,NRtol);
			}
		}
	}
	else if(elementType == H3D20N){ //for 20 Noded Hex (3D)
		bool insideBoundingBox = SearchElementBoundingBoxForPoint(pointCoord,tol);
		if (insideBoundingBox == true)
		{		
			isPointInside = SearchElementNodesForPoint(pointCoord,localCoord,tol);
			if(isPointInside ==false){
				isPointInside = InverseInterpolationForPointHex(pointCoord,localCoord,NRtol);
			}
		}
	}
		

	else {
		cout<<"IsoElement::SearchForPointInsideElement - Search not implemented for this element type. Exiting."<<endl;
		exit(1);
	}


	delete [] x; delete [] y; delete [] z;
	return isPointInside;
	
	}
//=========================================================================	
bool IsoElement::SearchElementNodesForPoint(double *pointCoord, double *localCoord, double tol)
{
bool pointFound = false;
int numNodes;
int numDim;
double local[100];

	//Define local coordinate of nodes for element types
	//(hopefully this will get cleaned up as we define these elsewhere in the code later)
	//double E2D8N[16] = {-1.0,-1.0,
	//					0.0,-1.0,
	//					1.0,-1.0,
	//					1.0,0.0,
	//					1.0,1.0,
	//					0.0,1.0,
	//					-1.0,1.0,
	//					-1.0,0.0};
	//	double E3D20N[60] = {-1,-1,-1,
 //                        -1, 0,-1,
 //                        -1, 1,-1,
 //                        -1, 1, 0,
 //                        -1, 1, 1,
 //                        -1, 0, 1,
 //                        -1,-1, 1,
 //                        -1,-1, 0,
 //                         0,-1,-1,
 //                         0, 1,-1,
 //                         0, 1, 1,
 //                         0,-1, 1,
 //                         1,-1,-1,
 //                         1, 0,-1,
 //                         1, 1,-1,
 //                         1, 1, 0,
 //                         1, 1, 1,
 //                         1, 0, 1,
 //                         1,-1, 1,
 //                         1,-1, 0};

	if(elementType == Q2D8N){
		numNodes = 8;
		numDim = 2;
        getMasterNodeCoords(local);
		//local = E2D8N;
		
	}
	else if(elementType == H3D20N){
		numNodes = 20;
		numDim = 3;
        getMasterNodeCoords(local);
		//local = E3D20N;
	}
	else{
		cout<<"IsoElement::SearchElementNodesForPoint - element type not implemented"<<endl;
		exit(1);}
				
	double *x,*y,*z;
	double xp,yp,zp;
	x = new double [numNodes];
	y = new double [numNodes];
	z = new double [numNodes];
	extractNodalCoordinates(x,y,z);
	xp = pointCoord[0]; yp=pointCoord[1]; zp=pointCoord[2];
	
	for(int i=0;i<numNodes;i++){
	if (fabs(xp-x[i]) < tol && fabs(yp-y[i]) < tol && fabs(zp-z[i]) < tol)
		{
//			cout << "The point is coincident with node " << i << endl;
			if(numDim==2){
				localCoord[0] = local[i*2];
				localCoord[1] = local[i*2+1];
				localCoord[2] = 0.0;
			}
			else if(numDim == 3){
				localCoord[0] = local[i*3];
				localCoord[1] = local[i*3+1];
				localCoord[2] = local[i*3+2];
			}
			pointFound=true;
			break;
		}
	}
	
	delete [] x; 	delete [] y;	delete [] z;

return pointFound;
}
//=========================================================================
bool IsoElement::SearchElementBoundingBoxForPoint(double *pointCoord, double tol)
{
	//See beta_docs/searchInsideElementForPoint.docx for documentation
	bool isInsideBoundingBox = false;
	double xp = pointCoord[0]; double yp = pointCoord[1]; double zp = pointCoord[2];
	int numNodes = getNumNodesPerElement();
	double *x; double *y; double *z;
	x = new double [numNodes];
	y = new double [numNodes];
	z = new double [numNodes];
	extractNodalCoordinates(x,y,z);

	double min_x, min_y, min_z, max_x, max_y, max_z;
	min_x = max_x = min_y = max_y = min_z = max_z = 0.0;
	getDoubleMinMax(x,numNodes,min_x,max_x);
	getDoubleMinMax(y,numNodes,min_y,max_y);
	getDoubleMinMax(z,numNodes,min_z,max_z);
	
	//resize the bounding box to allow for the possibility of badly shaped elements
	double fac = 0.1;
	double xfac = fac*fabs(max_x - min_x);
	double yfac = fac*fabs(max_y - min_y);
	double zfac = fac*fabs(max_z - min_z);
	
	min_x = min_x - xfac; max_x = max_x + xfac;
	min_y = min_y - yfac; max_y = max_y + yfac;
	min_z = min_z - zfac; max_z = max_z + zfac;
	//Checks to see if points lie on edge of bounding box within tol
	if(fabs(xp - min_x)<tol){
		xp = min_x;}
	else if(fabs(xp-max_x)<tol){
		xp = max_x;}
	if(fabs(yp - min_y)<tol){
		yp = min_y;}
	else if(fabs(yp-max_y)<tol){
		yp = max_y;}
	if(fabs(zp - min_z)<tol){
		zp = min_z;}
	else if(fabs(zp-max_z)<tol){
		zp = max_z;}

	if (xp <= max_x && xp >= min_x && yp <= max_y && yp >= min_y && zp<= max_z && zp>=min_z)
	{
		isInsideBoundingBox = true;
	}

delete [] x; delete [] y; delete [] z;

return isInsideBoundingBox;
}
//=========================================================================
bool IsoElement::InverseInterpolationForPointQuad(double *pointCoord,double *localCoord,double tol)
{
//See beta_docs/searchInsideElementForPoint.docx for documentation
bool isPointInside = false;
double *x; double *y; double *z;
int numNodes = getNumNodesPerElement();
x = new double [numNodes];
y = new double [numNodes];
z = new double [numNodes];
extractNodalCoordinates(x,y,z);
double xp = pointCoord[0];
double yp = pointCoord[1];
double xi=0.0; double eta=0.0;

//Newton Raphson Solver
double *S; double *Sxi; double *Seta;
S = new double [numNodes];
Sxi = new double [numNodes];
Seta = new double [numNodes];
double Residual[2];
int maxIteration = 20;
int iterationCount = 0;
double error=2*tol;
double dxi,deta;

	while(error>tol && iterationCount < maxIteration)
		{
		//calculate residual
		shape2D(numNodes,xi,eta,S);
		deriv2D(numNodes,xi,eta,Sxi,Seta);
		Residual[0] = xp;  Residual[1] = yp;
		for (int i=0;i<numNodes;i++){
			Residual[0] -= x[i]*S[i]; 
			Residual[1] -= y[i]*S[i];
		}

		error = sqrt(Residual[0]*Residual[0] + Residual[1]*Residual[1]);
		if(error<=tol)
		{
			break;
		}

		//calculate jacobian
		double J[2][2], invJ[2][2];
              J[0][0] = 0.;
              J[0][1] = 0.;
              J[1][0] = 0.;
              J[1][1] = 0.;
		double temp;
			
		for(int i=0; i<numNodes; i++)
		{
		J[0][0] += Sxi[i] * x[i];
		J[0][1] += Seta[i] * x[i];
		J[1][0] += Sxi[i] * y[i];
		J[1][1] += Seta[i] * y[i];
		}
			// Determinant of J and its inverse
        double detJ = J[0][0]*J[1][1]-
                       J[0][1]*J[1][0];
        temp = 1/detJ;
        invJ[0][0]= J[1][1] * temp;
        invJ[0][1]=-J[0][1] * temp;
        invJ[1][0]=-J[1][0] * temp;
        invJ[1][1]= J[0][0] * temp;

		dxi = invJ[0][0]*Residual[0] + invJ[0][1]*Residual[1];
		deta = invJ[1][0]*Residual[0] + invJ[1][1]*Residual[1];

		xi += dxi;
		eta += deta;

		++iterationCount;
		}

		if (fabs(xi-1.0)<tol){
			xi = 1.0;
		}
		else if (fabs(xi+1.0)<tol){
			xi = -1.0;
		}
		if (fabs(eta-1.0)<tol){
			eta = 1.0;
		}
		else if (fabs(eta+1.0)<tol){
			eta = -1.0;
		}
		//check local coordinate for inside/outside element
		if(xi<=1.0 && xi>=-1.0 && eta<=1.0 && eta>=-1.0)
		{
			isPointInside = true;
		}
		else{
			isPointInside = false;
			if(iterationCount >=maxIteration){
					BETA_OUT<<"Element #: "<<elementNumber;
					BETA_OUT<<"  MAXIMUM ITERATION COUNT EXCEEDED IN SEARCH FOR POINT"<<endl;}
			}

localCoord[0] = xi; localCoord[1] = eta; localCoord[2] = 0.0;
delete [] x; delete [] y; delete [] z;
delete [] S; delete [] Sxi; delete []Seta;
return isPointInside;
}
//=========================================================================
bool IsoElement::InverseInterpolationForPointHex(double *pointCoord,double* localCoord,double tol)
{
	//See beta_docs/searchInsideElementForPoint.docx for documentation
	bool isPointInside = false;
	double *x; double *y; double *z;
	int numNodes = getNumNodesPerElement();
	x = new double [numNodes];
	y = new double [numNodes];
	z = new double [numNodes];
	extractNodalCoordinates(x,y,z);
	double xp = pointCoord[0];
	double yp = pointCoord[1];
	double zp = pointCoord[2];
	double xi=0.0; double eta=0.0; double zeta=0.0;

	double *S; double *Sxi; double *Seta; double *Szeta;
	S = new double [numNodes];
	Sxi = new double [numNodes];
	Seta = new double [numNodes];
	Szeta = new double [numNodes];
	double Residual[3];
	int maxIteration = 20;
	int iterationCount = 0;
	double error=2*tol;
	double dxi,deta,dzeta;

			while(error>tol && iterationCount < maxIteration)
			{
			//calculate residual
			shape3D(numNodes,xi,eta,zeta,S);
			deriv3D(numNodes,xi,eta,zeta,Sxi,Seta,Szeta);
			Residual[0] = xp;  Residual[1] = yp; Residual[2] = zp;
			for (int i=0;i<numNodes;i++){
				Residual[0] -= x[i]*S[i]; 
				Residual[1] -= y[i]*S[i];
				Residual[2] -= z[i]*S[i];
			}

			error = sqrt(Residual[0]*Residual[0] + Residual[1]*Residual[1] + Residual[2]*Residual[2]);
			if(error<=tol)
			{
				break;
			}

			//calculate jacobian
			double J[3][3], invJ[3][3];
              J[0][0] = 0.;
              J[0][1] = 0.;
			  J[0][2] = 0.;
              J[1][0] = 0.;
              J[1][1] = 0.;
			  J[1][2] = 0.;
			  J[2][0] = 0.;
              J[2][1] = 0.;
			  J[2][2] = 0.;
			double temp;
			
			for(int i=0; i<numNodes; i++)
			{
			J[0][0] += Sxi[i] * x[i];
			J[0][1] += Seta[i] * x[i];
			J[0][2] += Szeta[i] * x[i];
			J[1][0] += Sxi[i] * y[i];
			J[1][1] += Seta[i] * y[i];
			J[1][2] += Szeta[i] * y[i];
			J[2][0] += Sxi[i] * z[i];
			J[2][1] += Seta[i] * z[i];
			J[2][2] += Szeta[i] * z[i];
			}
			// Determinant of J and its inverse
             double detJ = J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][0]
			              +J[0][2]*J[1][0]*J[2][1] - J[0][0]*J[1][2]*J[2][1]
						  -J[0][1]*J[1][0]*J[2][2] - J[0][2]*J[1][1]*J[2][0];
             temp = 1/detJ;
             invJ[0][0]= (J[1][1]*J[2][2]-J[1][2]*J[2][1])* temp;
             invJ[0][1]= (J[0][2]*J[2][1]-J[0][1]*J[2][2]) * temp;
			 invJ[0][2]= (J[0][1]*J[1][2]-J[0][2]*J[1][1]) * temp;
             invJ[1][0]= (J[1][2]*J[2][0] - J[1][0]*J[2][2])* temp;
             invJ[1][1]= (J[0][0]*J[2][2] - J[0][2]*J[2][0]) * temp;
			 invJ[1][2]= (J[0][2]*J[1][0] - J[0][0]*J[1][2]) * temp;
			 invJ[2][0]= (J[1][0]*J[2][1] - J[1][1]*J[2][0])* temp;
             invJ[2][1]= (J[0][1]*J[2][0] - J[0][0]*J[2][1]) * temp;
			 invJ[2][2]= (J[0][0]*J[1][1] - J[0][1]*J[1][0]) * temp;

			 dxi = invJ[0][0]*Residual[0] + invJ[0][1]*Residual[1] + invJ[0][2]*Residual[2];
			 deta = invJ[1][0]*Residual[0] + invJ[1][1]*Residual[1] + invJ[1][2]*Residual[2];
			 dzeta = invJ[2][0]*Residual[0] + invJ[2][1]*Residual[1] + invJ[2][2]*Residual[2];

			 xi += dxi;
			 eta += deta;
			 zeta += dzeta;

			 ++iterationCount;
			}

			if (fabs(xi-1.0)<tol){
				xi = 1.0;
			}
			else if (fabs(xi+1.0)<tol){
				xi = -1.0;
			}
			if (fabs(eta-1.0)<tol){
				eta = 1.0;
			}
			else if (fabs(eta+1.0)<tol){
				eta = -1.0;
			}
			if (fabs(zeta-1.0)<tol){
				zeta =1.0;
			}
			else if (fabs(zeta+1.0)<tol){
				zeta = -1.0;
			}
			//check local coordinate for inside/outside element
			if(xi<=1.0 && xi>=-1.0 && eta<=1.0 && eta>=-1.0)
			{
				isPointInside = true;
			}
			else{
				isPointInside = false;
				if(iterationCount >=maxIteration){
					BETA_OUT<<"Element #: "<<elementNumber;
					BETA_OUT<<"  MAXIMUM ITERATION COUNT EXCEEDED IN SEARCH FOR POINT"<<endl;}

			}

delete [] x; delete [] y; delete []z;
delete [] S; delete [] Sxi; delete [] Seta; delete [] Szeta;

localCoord[0] = xi; localCoord[1] = eta; localCoord[2] = zeta;
return isPointInside;
}
//=========================================================================

