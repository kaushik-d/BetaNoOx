#include"stdafx.h"

#include"Q3DDamageElement.hpp"
#include "math/gaussQuadrature/quadQuadPoint.hpp"
#include "elements/InterpAndSolDerivatives.hpp" 


Q3DDamageElement::Q3DDamageElement(void):
DMElasticElement()
{   
    numberOfIntegrationDirections=2;
    setIntegrationOrder(2);
}

//==================================================================

Q3DDamageElement::~Q3DDamageElement(void)
{  
}

//==================================================================

void Q3DDamageElement::calculateBmatrix(InterpAndSolDerivatives *id,
                                        double **b)
{
    int numberOfInterp = id->numberOfInterp;

    //Order: xx yy zz xy yz xz (Engineering shear)
    double *p_S_x2=id->p_S_x1,     
           *p_S_x3=id->p_S_x2  ;  //Note shifting of coor.
    double *S = id->S;
    int i;
    int c1 = 0;
    int c2 = 1;
    int c3 = 2;
    int c4 = 3;

    for( i=0; i< numberOfInterp; i++, c1+=4, c2+=4, c3+=4, c4+=4)
    {
        b[0][c1] = 0.       ; b[0][c2] = 0.       ; b[0][c3] = 0.        ;b[0][c4] = S[i];
        b[1][c1] = 0.       ; b[1][c2] = p_S_x2[i]; b[1][c3] = 0.        ;b[1][c4] = 0.;
        b[2][c1] = 0.       ; b[2][c2] = 0.       ; b[2][c3] = p_S_x3[i] ;b[2][c4] = 0.;

        b[3][c1] = p_S_x2[i]; b[3][c2] = 0.       ; b[3][c3] = 0.        ;b[3][c4] = 0.;
        b[4][c1] = 0.       ; b[4][c2] = p_S_x3[i]; b[4][c3] = p_S_x2[i] ;b[4][c4] = 0.;
        b[5][c1] = p_S_x3[i]; b[5][c2] = 0.       ; b[5][c3] = 0.        ;b[5][c4] = 0.;
    }
}//end of linearBmat

//==================================================================

void Q3DDamageElement::calculateStrain(InterpAndSolDerivatives *id,
                                       double *strain)
{//engineering shear
     strain[0] = id->p_u1_x1 ;
     strain[1] = id->p_u2_x2 ;
     strain[2] = id->p_u3_x3 ;
     strain[3] = id->p_u1_x2 ;
     strain[4] = id->p_u2_x3 + id->p_u3_x2;
     strain[5] = id->p_u1_x3 ;
     
     //BETA_OUT<<"Strains"<<endl;
     //for(int i=0; i<6; i++){BETA_OUT<<strain[i]<<"  ";}BETA_OUT<<endl;
}

//==================================================================

void Q3DDamageElement::derivativesOfFields(InterpAndSolDerivatives *id)
{
    int i;
    int numberOfInterp = id->numberOfInterp;
    double *S=id->S, *p_S_Lx1=id->p_S_Lx1, 
                    *p_S_Lx2=id->p_S_Lx2;

    double *p_S_x2=id->p_S_x1,     *p_S_x3= id->p_S_x2;  //Note the shift in coordinates

    double *xDisp=id->nodalValues_u1, *yDisp=id->nodalValues_u2,
           *zDisp=id->nodalValues_u3,
           *u4   =id->nodalValues_u4;

    double p_u1_x1;
    double p_u1_x2, p_u1_x3;
    double p_u2_x2, p_u2_x3;
    double p_u3_x2, p_u3_x3;

    p_u1_x1 =           0.;
    p_u1_x2 = p_u1_x3 = 0.;
    p_u2_x2 = p_u2_x3 = 0.;
    p_u3_x2 = p_u3_x3 = 0.;

    for(i=0; i< numberOfInterp; i++)
    {
        p_u1_x1 += S[i]      * u4[i];
        p_u1_x2 += p_S_x2[i] * xDisp[i];
        p_u1_x3 += p_S_x3[i] * xDisp[i];
        p_u2_x2 += p_S_x2[i] * yDisp[i];
        p_u2_x3 += p_S_x3[i] * yDisp[i];
        p_u3_x2 += p_S_x2[i] * zDisp[i];
        p_u3_x3 += p_S_x3[i] * zDisp[i];
    }
    id->p_u1_x1 = p_u1_x1;
    id->p_u1_x2 = p_u1_x2;
    id->p_u1_x3 = p_u1_x3;
    id->p_u2_x2 = p_u2_x2;
    id->p_u2_x3 = p_u2_x3;
    id->p_u3_x2 = p_u3_x2;
    id->p_u3_x3 = p_u3_x3;
}  // end of derivativesOfFields

//==================================================================

void Q3DDamageElement::update(char * action, double *uGlobal)
{
    IntegrationDataType   idt;
    idt.id=new InterpAndSolDerivatives(numNodesPerElement,numNodesPerElement); idt.id->allocate();

    idt.elemDisp[0] = idt.id->nodalValues_u1; 
    idt.elemDisp[1] = idt.id->nodalValues_u2;
    idt.elemDisp[2] = idt.id->nodalValues_u3;  
    idt.elemDisp[3] = idt.id->nodalValues_u4;  //jdw  This is the only change in base class method
    idt.id->numberOfInterp = numNodesPerElement;
    idt.uGlobal=uGlobal;

    extractNodalCoordinates(idt.id->xCoor,idt.id->yCoor,idt.id->zCoor);
    if( uGlobal != 0) {extractSolution( uGlobal, idt.elemDisp);} 
    getQuadraturePoints(idt.gaussPointList); 

    //Will be allocated if mangles has been called.
    if(ElementNodeRotations == 0){AllocateRotations();}

    loopOverIntegrationPoints(action, idt);
}