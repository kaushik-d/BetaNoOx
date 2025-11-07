#include "stdafx.h"
using namespace std;
#undef BETA_OUT
#define BETA_OUT cout


//=====================================================
//                 2-node line
//=====================================================
void shapeLine2( double * S, double zeta)
{
      S[0] = 1.0/2.0-zeta/2;
      S[1] = 1.0/2.0+zeta/2;
}
void derivLine2( double * Szeta, double zeta)
{
      Szeta[0] = -1.0/2.0;
      Szeta[1] = 1.0/2.0;
}
//=====================================================
//                 3-node line
//=====================================================

void shapeLine3( double * S, double zeta)
{
 double t1;
      t1 = zeta*zeta;
      S[0] = -zeta/2+t1/2;
      S[1] = 1.0-t1;
      S[2] = zeta/2+t1/2;
}
void derivLine3( double * Szeta, double zeta)
{
      Szeta[0] = -1.0/2.0+zeta;
      Szeta[1] = -2.0*zeta;
      Szeta[2] = 1.0/2.0+zeta;
}

//=====================================================
//                 4-node line
//=====================================================

void shapeLine4( double * S, double zeta)
{double t1,t2;
      t1 = zeta*zeta;
      t2 = t1*zeta;
      S[0] = -1.0/16.0+zeta/16+9.0/16.0*t1-9.0/16.0*t2;
      S[1] = 9.0/16.0-27.0/16.0*zeta-9.0/16.0*t1+27.0/16.0*t2;
      S[2] = 9.0/16.0+27.0/16.0*zeta-9.0/16.0*t1-27.0/16.0*t2;
      S[3] = -1.0/16.0-zeta/16+9.0/16.0*t1+9.0/16.0*t2;
}
void derivLine4( double * Szeta, double zeta)
{double t1;  
      t1 = zeta*zeta;
      Szeta[0] = 1.0/16.0+9.0/8.0*zeta-27.0/16.0*t1;
      Szeta[1] = -27.0/16.0-9.0/8.0*zeta+81.0/16.0*t1;
      Szeta[2] = 27.0/16.0-9.0/8.0*zeta-81.0/16.0*t1;
      Szeta[3] = -1.0/16.0+9.0/8.0*zeta+27.0/16.0*t1;
}
//=====================================================
//                 5-node line
//=====================================================

void shapeLine5( double * S, double zeta)
{double t1,t2,t3;
      t1 = zeta*zeta;
      t2 = t1*zeta;
      t3 = t1*t1;
      S[0] = zeta/6-t1/6-2.0/3.0*t2+2.0/3.0*t3;
      S[1] = -4.0/3.0*zeta+8.0/3.0*t1+4.0/3.0*t2-8.0/3.0*t3;
      S[2] = 1.0-5.0*t1+4.0*t3;
      S[3] = 4.0/3.0*zeta+8.0/3.0*t1-4.0/3.0*t2-8.0/3.0*t3;
      S[4] = -zeta/6-t1/6+2.0/3.0*t2+2.0/3.0*t3;
}
void derivLine5( double * Szeta, double zeta)
{double t1,t2; 
      t1 = zeta*zeta;
      t2 = t1*zeta;
      Szeta[0] = 1.0/6.0-zeta/3-2.0*t1+8.0/3.0*t2;
      Szeta[1] = -4.0/3.0+16.0/3.0*zeta+4.0*t1-32.0/3.0*t2;
      Szeta[2] = -10.0*zeta+16.0*t2;
      Szeta[3] = 4.0/3.0+16.0/3.0*zeta-4.0*t1-32.0/3.0*t2;
      Szeta[4] = -1.0/6.0-zeta/3+2.0*t1+8.0/3.0*t2;
}

//=====================================================
//                 3-node triangle
//=====================================================
void shapeTri3( double * S, double zeta, double eta)
{
    S[0] = 1 - zeta - eta;
    S[1] = zeta;
    S[2] = eta;
}

void derivTri3(double * Szeta      ,
               double * Seta       ,
               double zeta, double eta)
{
    Szeta[0] = -1;
    Szeta[1] =  1;
    Szeta[2] =  0;

    Seta[0]  = -1;
    Seta[1]  =  0;
    Seta[2]  =  1;
}

//=====================================================
//                 6-node triangle
//=====================================================
void shapeTri6( double * S, double xi, double eta)
{
    double t1 = xi*xi;
    double t2 = eta*eta;
    double t3 = eta*xi;
    S[0] = 1-3*xi-3*eta+4*t3+2*t1+2*t2;
    S[1] = 4*xi-4*t3-4*t1;
    S[2] = -xi+2*t1;
    S[3] = 4*t3;
    S[4] = -eta+2*t2;
    S[5] = 4*eta-4*t3-4*t2;
}

void derivTri6(double * Sxi        ,
               double * Seta       ,
               double xi, double eta)
{
    Sxi[0] = -3+4*eta+4*xi;
    Sxi[1] = 4-4*eta-8*xi;
    Sxi[2] = -1+4*xi;
    Sxi[3] = 4*eta;
    Sxi[4] = 0;
    Sxi[5] = -4*eta;

    Seta[0]  = -3+4*eta+4*xi;
    Seta[1]  = -4*xi;
    Seta[2]  = 0;
    Seta[3]  = 4*xi;
    Seta[4]  = -1+4*eta;
    Seta[5]  = 4-4*xi-8*eta;
}

//=====================================================
//                 4-node quad
//=====================================================
void shapeQuad4( double * S, double zeta, double eta)
{
		S[0] = .25 * (1-zeta) * (1-eta);
		S[1] = .25 * (1+zeta) * (1-eta);
		S[2] = .25 * (1+zeta) * (1+eta);
		S[3] = .25 * (1-zeta) * (1+eta);
}

void derivQuad4(double * Szeta      ,
					  double * Seta       ,
					  double zeta, double eta)
{
 Szeta[0] = -.25*(1-eta);
 Szeta[1] =  .25*(1-eta);
 Szeta[2] =  .25*(1+eta);
 Szeta[3] = -.25*(1+eta);

 Seta[0]  = -.25*(1-zeta);
 Seta[1]  = -.25*(1+zeta);
 Seta[2]  =  .25*(1+zeta);
 Seta[3]  =  .25*(1-zeta);
}

//=====================================================
//                 8-node quad
//=====================================================

void shapeQuad8( double * S, double xi, double eta)
{
S[0]= (-1 + eta)*(1 + eta - eta*xi - xi*xi)/4;
S[1]= (-1 + eta)*(-1 + xi*xi)/2;
S[2]= (-1 + eta)*(1 + eta + eta*xi - xi*xi)/4;
S[3]= (1 - eta*eta)*(1 + xi)/2;
S[4]= (1 + eta)*(-1 + eta + eta*xi + xi*xi)/4;
S[5]= (1 + eta)*(1 - xi*xi)/2;
S[6]= (1 + eta)*(-1 + eta - eta*xi + xi*xi)/4;
S[7]= (-1 + eta*eta)*(-1 + xi)/2;
}

void derivQuad8(double * Szeta      ,
					double * Seta       ,
					double zeta, double eta)
{
		 Szeta[0] =  -((1-eta)*(1-zeta))/4 - (1 - eta)*(-1 - eta - zeta)/4 ;
		 Szeta[1] =   -((1 - eta)*zeta) ;
		 Szeta[2] =   (1-eta)*(1 + zeta)/4 + (1 - eta)*(-1 - eta + zeta)/4  ;
		 Szeta[3] =   (1 - eta*eta)/2  ;
		 Szeta[4] =   (1+eta)*(1 + zeta)/4 + (1 + eta)*(-1 + eta + zeta)/4 ;
		 Szeta[5] =   -((1 + eta)*zeta) ;
		 Szeta[6] =   -((1+eta)*(1-zeta))/4 - (1+eta)*(-1 + eta - zeta)/4  ;
		 Szeta[7] =   -(1 - eta*eta)/2  ;
		 Seta[0] =   -((1-eta)*(1-zeta))/4 - (1- zeta)*(-1 - eta - zeta)/4 ;
		 Seta[1] =   -(1 - zeta*zeta)/2  ;
		 Seta[2] =   -((1-eta)*(1+zeta))/4 - (1+ zeta)*(-1 - eta + zeta)/4  ;
		 Seta[3] =   -(eta*(1 + zeta))  ;
		 Seta[4] =   (1+eta)*(1+zeta)/4 + (1+ zeta)*(-1 + eta + zeta)/4   ;
		 Seta[5] =   (1 - zeta*zeta)/2  ;
		 Seta[6] =   (1+eta)*(1 - zeta)/4 + (1 - zeta)*(-1 + eta - zeta)/4  ;
		 Seta[7] =   -(eta*(1 - zeta))  ;
}
//======================================================

//======================================================
void shapeHex8( double * S, double xi, double eta, double zeta)
{
double t1, t2, t3, t4;
      t1 = xi*eta;
      t2 = xi*zeta;
      t3 = eta*zeta;
      t4 = t1*zeta;
      S[0] = 1.0/8.0-xi/8-eta/8+t1/8-zeta/8+t2/8+t3/8-t4/8;
      S[1] = 1.0/8.0-xi/8+eta/8-t1/8-zeta/8+t2/8-t3/8+t4/8;
      S[2] = 1.0/8.0-xi/8+eta/8-t1/8+zeta/8-t2/8+t3/8-t4/8;
      S[3] = 1.0/8.0-xi/8-eta/8+t1/8+zeta/8-t2/8-t3/8+t4/8;
      S[4] = 1.0/8.0+xi/8-eta/8-t1/8-zeta/8-t2/8+t3/8+t4/8;
      S[5] = 1.0/8.0+xi/8+eta/8+t1/8-zeta/8-t2/8-t3/8-t4/8;
      S[6] = 1.0/8.0+xi/8+eta/8+t1/8+zeta/8+t2/8+t3/8+t4/8;
      S[7] = 1.0/8.0+xi/8-eta/8-t1/8+zeta/8+t2/8-t3/8-t4/8;

}

void derivHex8(double * Sxi, double * Seta, double * Szeta,
	        double xi,  double eta, double zeta)
{
double t1;
      t1 = eta*zeta;
      Sxi[0] = -1.0/8.0+eta/8+zeta/8-t1/8;
      Sxi[1] = -1.0/8.0-eta/8+zeta/8+t1/8;
      Sxi[2] = -1.0/8.0-eta/8-zeta/8-t1/8;
      Sxi[3] = -1.0/8.0+eta/8-zeta/8+t1/8;
      Sxi[4] = 1.0/8.0-eta/8-zeta/8+t1/8;
      Sxi[5] = 1.0/8.0+eta/8-zeta/8-t1/8;
      Sxi[6] = 1.0/8.0+eta/8+zeta/8+t1/8;
      Sxi[7] = 1.0/8.0-eta/8+zeta/8-t1/8;
      t1 = xi*zeta;
      Seta[0] = -1.0/8.0+xi/8+zeta/8-t1/8;
      Seta[1] = 1.0/8.0-xi/8-zeta/8+t1/8;
      Seta[2] = 1.0/8.0-xi/8+zeta/8-t1/8;
      Seta[3] = -1.0/8.0+xi/8-zeta/8+t1/8;
      Seta[4] = -1.0/8.0-xi/8+zeta/8+t1/8;
      Seta[5] = 1.0/8.0+xi/8-zeta/8-t1/8;
      Seta[6] = 1.0/8.0+xi/8+zeta/8+t1/8;
      Seta[7] = -1.0/8.0-xi/8-zeta/8-t1/8;
      t1 = xi*eta;
      Szeta[0] = -1.0/8.0+xi/8+eta/8-t1/8;
      Szeta[1] = -1.0/8.0+xi/8-eta/8+t1/8;
      Szeta[2] = 1.0/8.0-xi/8+eta/8-t1/8;
      Szeta[3] = 1.0/8.0-xi/8-eta/8+t1/8;
      Szeta[4] = -1.0/8.0-xi/8+eta/8+t1/8;
      Szeta[5] = -1.0/8.0-xi/8-eta/8-t1/8;
      Szeta[6] = 1.0/8.0+xi/8+eta/8+t1/8;
      Szeta[7] = 1.0/8.0+xi/8-eta/8-t1/8;

}


//======================================================
void shapeHex20( double * S, double xi, double eta, double zeta)
{
double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,
       t15,t16,t18,t24,t29,t31,t33,t35,
		 t20,t22;

//S =============================================
      t1 = xi*xi;
      t2 = eta*eta;
      t3 = zeta*zeta;
      t4 = xi*t2;
      t5 = t1*eta;
      t6 = eta*t3;
      t7 = t2*zeta;
      t8 = xi*t3;
      t9 = t1*zeta;
      t10 = xi*eta;
      t11 = t10*t3;
      t12 = t10*zeta;
      t13 = t4*zeta;
      t14 = t5*zeta;
      t15 = -1.0/4.0+t1/8+eta/8+t2/8+t3/8+xi/8+zeta/8-t4/8-t5/8-t6/8-t7/8-t8/8-
t9/8+t11/8-t12/8+t13/8+t14/8;
      t16 = xi*zeta;
      t18 = -1.0/4.0+t1/8-eta/8+t2/8+t3/8+xi/8+zeta/8-t4/8+t5/8+t6/8-t7/8-t8/8-
t9/8-t11/8+t12/8+t13/8-t14/8;
      t20 = -1.0/4.0+t1/8-eta/8+t2/8+t3/8+xi/8-zeta/8-t4/8+t5/8+t6/8+t7/8-t8/8+
t9/8-t11/8-t12/8-t13/8+t14/8;
      t22 = -1.0/4.0+t1/8+eta/8+t2/8+t3/8+xi/8-zeta/8-t4/8-t5/8-t6/8+t7/8-t8/8+
t9/8+t11/8+t12/8-t13/8-t14/8;
      t24 = eta*zeta;
      t29 = -1.0/4.0+t1/8+eta/8+t2/8+t3/8-xi/8+zeta/8+t4/8-t5/8-t6/8-t7/8+t8/8-
t9/8-t11/8+t12/8-t13/8+t14/8;
      t31 = -1.0/4.0+t1/8-eta/8+t2/8+t3/8-xi/8+zeta/8+t4/8+t5/8+t6/8-t7/8+t8/8-
t9/8+t11/8-t12/8-t13/8-t14/8;
      t33 = -1.0/4.0+t1/8-eta/8+t2/8+t3/8-xi/8-zeta/8+t4/8+t5/8+t6/8+t7/8+t8/8+
t9/8+t11/8+t12/8+t13/8+t14/8;
      t35 = -1.0/4.0+t1/8+eta/8+t2/8+t3/8-xi/8-zeta/8+t4/8-t5/8-t6/8+t7/8+t8/8+
t9/8-t11/8-t12/8+t13/8-t14/8;
      S[0] = t15;
      S[1] = 1.0/4.0-t2/4-xi/4-zeta/4+t4/4+t7/4+t16/4-t13/4;
      S[2] = t18;
      S[3] = 1.0/4.0+eta/4-t3/4-xi/4-t10/4-t6/4+t8/4+t11/4;
      S[4] = t20;
      S[5] = 1.0/4.0-t2/4-xi/4+zeta/4+t4/4-t7/4-t16/4+t13/4;
      S[6] = t22;
      S[7] = 1.0/4.0-eta/4-t3/4-xi/4+t10/4+t6/4+t8/4-t11/4;
      S[8] = 1.0/4.0-t1/4-eta/4-zeta/4+t5/4+t24/4+t9/4-t14/4;
      S[9] = 1.0/4.0-t1/4+eta/4-zeta/4-t5/4-t24/4+t9/4+t14/4;
      S[10] = 1.0/4.0-t1/4+eta/4+zeta/4-t5/4+t24/4-t9/4-t14/4;
      S[11] = 1.0/4.0-t1/4-eta/4+zeta/4+t5/4-t24/4-t9/4+t14/4;
      S[12] = t29;
      S[13] = 1.0/4.0-t2/4+xi/4-zeta/4-t4/4+t7/4-t16/4+t13/4;
      S[14] = t31;
      S[15] = 1.0/4.0+eta/4-t3/4+xi/4+t10/4-t6/4-t8/4-t11/4;
      S[16] = t33;
      S[17] = 1.0/4.0-t2/4+xi/4+zeta/4-t4/4-t7/4+t16/4-t13/4;
      S[18] = t35;
      S[19] = 1.0/4.0-eta/4-t3/4+xi/4-t10/4+t6/4-t8/4+t11/4;


}

void derivHex20(double * Sxi, double * Seta, double * Szeta,
                double xi,  double eta, double zeta)
{

double  t1, t2, t3,t4, t5, t6, t7, t8;

//Sxi =============================================
      t1 = eta*eta;
      t2 = xi*eta;
      t3 = zeta*zeta;
      t4 = xi*zeta;
      t5 = eta*t3;
      t6 = eta*zeta;
      t7 = t1*zeta;
      t8 = t2*zeta;
      Sxi[0] = xi/4+1.0/8.0-t1/8-t2/4-t3/8-t4/4+t5/8-t6/8+t7/8+t8/4;
      Sxi[1] = -1.0/4.0+t1/4+zeta/4-t7/4;
      Sxi[2] = xi/4+1.0/8.0-t1/8+t2/4-t3/8-t4/4-t5/8+t6/8+t7/8-t8/4;
      Sxi[3] = -1.0/4.0-eta/4+t3/4+t5/4;
      Sxi[4] = xi/4+1.0/8.0-t1/8+t2/4-t3/8+t4/4-t5/8-t6/8-t7/8+t8/4;
      Sxi[5] = -1.0/4.0+t1/4-zeta/4+t7/4;
      Sxi[6] = xi/4+1.0/8.0-t1/8-t2/4-t3/8+t4/4+t5/8+t6/8-t7/8-t8/4;
      Sxi[7] = -1.0/4.0+eta/4+t3/4-t5/4;
      Sxi[8] = -xi/2+t2/2+t4/2-t8/2;
      Sxi[9] = -xi/2-t2/2+t4/2+t8/2;
      Sxi[10] = -xi/2-t2/2-t4/2-t8/2;
      Sxi[11] = -xi/2+t2/2-t4/2+t8/2;
      Sxi[12] = xi/4-1.0/8.0+t1/8-t2/4+t3/8-t4/4-t5/8+t6/8-t7/8+t8/4;
      Sxi[13] = 1.0/4.0-t1/4-zeta/4+t7/4;
      Sxi[14] = xi/4-1.0/8.0+t1/8+t2/4+t3/8-t4/4+t5/8-t6/8-t7/8-t8/4;
      Sxi[15] = 1.0/4.0+eta/4-t3/4-t5/4;
      Sxi[16] = xi/4-1.0/8.0+t1/8+t2/4+t3/8+t4/4+t5/8+t6/8+t7/8+t8/4;
      Sxi[17] = 1.0/4.0-t1/4+zeta/4-t7/4;
      Sxi[18] = xi/4-1.0/8.0+t1/8-t2/4+t3/8+t4/4-t5/8-t6/8+t7/8-t8/4;
      Sxi[19] = 1.0/4.0-eta/4-t3/4+t5/4;

//Seta ============================================
      t1 = xi*eta;
      t2 = xi*xi;
      t3 = zeta*zeta;
      t4 = eta*zeta;
      t5 = xi*t3;
      t6 = xi*zeta;
      t7 = t1*zeta;
      t8 = t2*zeta;
      Seta[0] = 1.0/8.0+eta/4-t1/4-t2/8-t3/8-t4/4+t5/8-t6/8+t7/4+t8/8;
      Seta[1] = -eta/2+t1/2+t4/2-t7/2;
      Seta[2] = -1.0/8.0+eta/4-t1/4+t2/8+t3/8-t4/4-t5/8+t6/8+t7/4-t8/8;
      Seta[3] = 1.0/4.0-xi/4-t3/4+t5/4;
      Seta[4] = -1.0/8.0+eta/4-t1/4+t2/8+t3/8+t4/4-t5/8-t6/8-t7/4+t8/8;
      Seta[5] = -eta/2+t1/2-t4/2+t7/2;
      Seta[6] = 1.0/8.0+eta/4-t1/4-t2/8-t3/8+t4/4+t5/8+t6/8-t7/4-t8/8;
      Seta[7] = -1.0/4.0+xi/4+t3/4-t5/4;
      Seta[8] = -1.0/4.0+t2/4+zeta/4-t8/4;
      Seta[9] = 1.0/4.0-t2/4-zeta/4+t8/4;
      Seta[10] = 1.0/4.0-t2/4+zeta/4-t8/4;
      Seta[11] = -1.0/4.0+t2/4-zeta/4+t8/4;
      Seta[12] = 1.0/8.0+eta/4+t1/4-t2/8-t3/8-t4/4-t5/8+t6/8-t7/4+t8/8;
      Seta[13] = -eta/2-t1/2+t4/2+t7/2;
      Seta[14] = -1.0/8.0+eta/4+t1/4+t2/8+t3/8-t4/4+t5/8-t6/8-t7/4-t8/8;
      Seta[15] = 1.0/4.0+xi/4-t3/4-t5/4;
      Seta[16] = -1.0/8.0+eta/4+t1/4+t2/8+t3/8+t4/4+t5/8+t6/8+t7/4+t8/8;
      Seta[17] = -eta/2-t1/2-t4/2-t7/2;
      Seta[18] = 1.0/8.0+eta/4+t1/4-t2/8-t3/8+t4/4-t5/8-t6/8+t7/4-t8/8;
      Seta[19] = -1.0/4.0-xi/4+t3/4+t5/4;

//Szeta ===========================================
      t1 = eta*zeta;
      t2 = eta*eta;
      t3 = xi*zeta;
      t4 = xi*xi;
      t5 = xi*eta;
      t6 = t5*zeta;
      t7 = xi*t2;
      t8 = t4*eta;
      Szeta[0] = zeta/4+1.0/8.0-t1/4-t2/8-t3/4-t4/8+t6/4-t5/8+t7/8+t8/8;
      Szeta[1] = -1.0/4.0+t2/4+xi/4-t7/4;
      Szeta[2] = zeta/4+1.0/8.0+t1/4-t2/8-t3/4-t4/8-t6/4+t5/8+t7/8-t8/8;
      Szeta[3] = -zeta/2-t1/2+t3/2+t6/2;
      Szeta[4] = zeta/4-1.0/8.0+t1/4+t2/8-t3/4+t4/8-t6/4-t5/8-t7/8+t8/8;
      Szeta[5] = 1.0/4.0-t2/4-xi/4+t7/4;
      Szeta[6] = zeta/4-1.0/8.0-t1/4+t2/8-t3/4+t4/8+t6/4+t5/8-t7/8-t8/8;
      Szeta[7] = -zeta/2+t1/2+t3/2-t6/2;
      Szeta[8] = -1.0/4.0+eta/4+t4/4-t8/4;
      Szeta[9] = -1.0/4.0-eta/4+t4/4+t8/4;
      Szeta[10] = 1.0/4.0+eta/4-t4/4-t8/4;
      Szeta[11] = 1.0/4.0-eta/4-t4/4+t8/4;
      Szeta[12] = zeta/4+1.0/8.0-t1/4-t2/8+t3/4-t4/8-t6/4+t5/8-t7/8+t8/8;
      Szeta[13] = -1.0/4.0+t2/4-xi/4+t7/4;
      Szeta[14] = zeta/4+1.0/8.0+t1/4-t2/8+t3/4-t4/8+t6/4-t5/8-t7/8-t8/8;
      Szeta[15] = -zeta/2-t1/2-t3/2-t6/2;
      Szeta[16] = zeta/4-1.0/8.0+t1/4+t2/8+t3/4+t4/8+t6/4+t5/8+t7/8+t8/8;
      Szeta[17] = 1.0/4.0-t2/4+xi/4-t7/4;
      Szeta[18] = zeta/4-1.0/8.0-t1/4+t2/8+t3/4+t4/8-t6/4-t5/8+t7/8-t8/8;
      Szeta[19] = -zeta/2+t1/2-t3/2+t6/2;

}

void shapeTet10( double * S, double xi, double eta, double zeta)
{
	double t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15;

// S ==============================================
	t4 = eta * zeta;
	t5 = 0.4e1 * t4;
	t6 = zeta * xi;
	t7 = 0.4e1 * t6;
	t8 = xi * eta;
	t9 = 0.4e1 * t8;
	t10 = xi*xi;
	t11 = 0.2e1 * t10;
	t12 = eta * eta;
	t13 = 0.2e1 * t12;
	t14 = zeta * zeta;
	t15 = 0.2e1 * t14;
	S[0] = 0.1e1 - 0.3e1 * xi - 0.3e1 * eta - 0.3e1 * zeta + t5 + t7 + t9 + t11 + t13 + t15;
	S[1] = -xi + t11;
	S[2] = -eta + t13;
	S[3] = -zeta + t15;
	S[4] = 0.4e1 * xi - 0.4e1 * t6 - 0.4e1 * t8 - 0.4e1 * t10;
	S[5] = t9;
	S[6] = 0.4e1 * eta - 0.4e1 * t4 - 0.4e1 * t8 - 0.4e1 * t12;
	S[7] = 0.4e1 * zeta - 0.4e1 * t4 - 0.4e1 * t6 - 0.4e1 * t14;
	S[8] = t7;
	S[9] = t5;

}

void derivTet10(double * Sxi, double * Seta, double * Szeta,
                double xi,  double eta, double zeta)
{

double  t1, t2, t3;

//Sxi =============================================
	t1 = 0.4e1 * zeta;
	t2 = 0.4e1 * eta;
	t3 = 0.4e1 * xi;
	Sxi[0] = -0.3e1 + t1 + t2 + t3;
	Sxi[1] = -0.1e1 + t3;
	Sxi[2] = 0.0e0;
	Sxi[3] = 0.0e0;
	Sxi[4] = 0.4e1 - t1 - t2 - 0.8e1 * xi;
	Sxi[5] = t2;
	Sxi[6] = -t2;
	Sxi[7] = -t1;
	Sxi[8] = t1;
	Sxi[9] = 0.0e0;

//Seta ============================================
	t1 = 0.4e1 * zeta;
	t2 = 0.4e1 * eta;
	t3 = 0.4e1 * xi;
	Seta[0] = -0.3e1 + t1 + t2 + t3;
	Seta[1] = 0.0e0;
	Seta[2] = -0.1e1 + t2;
	Seta[3] = 0.0e0;
	Seta[4] = -t3;
	Seta[5] = t3;
	Seta[6] = 0.4e1 - t1 - t3 - 0.8e1 * eta;
	Seta[7] = -t1;
	Seta[8] = 0.0e0;
	Seta[9] = t1;

//Szeta ===========================================
	t1 = 0.4e1 * zeta;
	t2 = 0.4e1 * eta;
	t3 = 0.4e1 * xi;
	Szeta[0] = -0.3e1 + t1 + t2 + t3;
	Szeta[1] = 0.0e0;
	Szeta[2] = 0.0e0;
	Szeta[3] = -0.1e1 + t1;
	Szeta[4] = -t3;
	Szeta[5] = 0.0e0;
	Szeta[6] = -t2;
	Szeta[7] = 0.4e1 - t2 - t3 - 0.8e1 * zeta;
	Szeta[8] = t3;
	Szeta[9] = t2;

}

void shapeWedge15( double * S, double xi, double eta, double zeta)
{
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,
           t15,t16,t17,t18,t19,t20,t22,t25,t31;

// S ===================================================
	t1 = zeta / 0.2e1;
	t2 = eta * zeta;
	t3 = 0.3e1 / 0.2e1 * t2;
	t4 = zeta * xi;
	t5 = 0.3e1 / 0.2e1 * t4;
	t6 = xi * eta;
	t7 = 0.2e1 * t6;
	t8 = xi*xi;
	t9 = eta * eta;
	t10 = zeta * zeta;
	t11 = t10 / 0.2e1;
	t12 = t8 * zeta;
	t13 = t9 * zeta;
	t14 = t10 * xi;
	t15 = t14 / 0.2e1;
	t16 = t10 * eta;
	t17 = t16 / 0.2e1;
	t18 = t6 * zeta;
	t19 = 0.2e1 * t18;
	t20 = -xi - eta - t1 + t3 + t5 + t7 + t8 + t9 + t11 - t12 - t13 - t15 - t17 - t19;
	t22 = t4 / 0.2e1;
	t25 = t2 / 0.2e1;
	t31 = -xi - eta + t1 - t3 - t5 + t7 + t8 + t9 + t11 + t12 + t13 - t15 - t17 + t19;
	S[0] = t20;
	S[1] = 0.2e1 * xi - 0.2e1 * t4 - 0.2e1 * t6 - 0.2e1 * t8 + 0.2e1 * t12 + 0.2e1 * t18;
	S[2] = -xi + t22 + t8 - t12 + t15;
	S[3] = 0.2e1 * t6 - 0.2e1 * t18;
	S[4] = -eta + t25 + t9 - t13 + t17;
	S[5] = 0.2e1 * eta - 0.2e1 * t2 - 0.2e1 * t6 - 0.2e1 * t9 + 0.2e1 * t13 + 0.2e1 * t18;
	S[6] = 0.1e1 - xi - eta - t10 + t14 + t16;
	S[7] = xi - t14;
	S[8] = eta - t16;
	S[9] = t31;
	S[10] = 0.2e1 * xi + 0.2e1 * t4 - 0.2e1 * t6 - 0.2e1 * t8 - 0.2e1 * t12 - 0.2e1 * t18;
	S[11] = -xi - t22 + t8 + t12 + t15;
	S[12] = 0.2e1 * t6 + 0.2e1 * t18;
	S[13] = -eta - t25 + t9 + t13 + t17;
	S[14] = 0.2e1 * eta + 0.2e1 * t2 - 0.2e1 * t6 - 0.2e1 * t9 - 0.2e1 * t13 - 0.2e1 * t18;

}


void derivWedge15(double * Sxi, double * Seta, double * Szeta,
                double xi,  double eta, double zeta)
{

double  t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, 
        t11, t12, t13, t14, t15, t16, t17, t18, t20, t22;

//Sxi =============================================
t1 = 0.3e1 / 0.2e1 * zeta;
t2 = 0.2e1 * eta;
t3 = 0.2e1 * xi;
t4 = zeta * xi;
t5 = 0.2e1 * t4;
t6 = zeta * zeta;
t7 = t6 / 0.2e1;
t8 = eta * zeta;
t9 = 0.2e1 * t8;
t11 = 0.2e1 * zeta;
t12 = 0.4e1 * xi;
t13 = 0.4e1 * t4;
t15 = zeta / 0.2e1;
t17 = eta - t8;
t18 = -0.1e1 + t6;
t22 = eta + t8;
Sxi[0] = -0.1e1 + t1 + t2 + t3 - t5 - t7 - t9;
Sxi[1] = 0.2e1 - t11 - t2 - t12 + t13 + t9;
Sxi[2] = -0.1e1 + t15 + t3 - t5 + t7;
Sxi[3] = 0.2e1 * t17;
Sxi[4] = 0.0e0;
Sxi[5] = -0.2e1 * t17;
Sxi[6] = t18;
Sxi[7] = -t18;
Sxi[8] = 0.0e0;
Sxi[9] = -0.1e1 - t1 + t2 + t3 + t5 - t7 + t9;
Sxi[10] = 0.2e1 + t11 - t2 - t12 - t13 - t9;
Sxi[11] = -0.1e1 - t15 + t3 + t5 + t7;
Sxi[12] = 0.2e1 * t22;
Sxi[13] = 0.0e0;
Sxi[14] = -0.2e1 * t22;


//Seta ============================================
	t1 = 0.3e1 / 0.2e1 * zeta;
	t2 = 0.2e1 * eta;
	t3 = 0.2e1 * xi;
	t4 = zeta * xi;
	t5 = 0.2e1 * t4;
	t6 = zeta * zeta;
	t7 = t6 / 0.2e1;
	t8 = eta * zeta;
	t9 = 0.2e1 * t8;
	t11 = -xi + t4;
	t12 = zeta / 0.2e1;
	t14 = 0.2e1 * zeta;
	t15 = 0.4e1 * eta;
	t16 = 0.4e1 * t8;
	t18 = -0.1e1 + t6;
	t20 = -xi - t4;
	Seta[0] = -0.1e1 + t1 + t2 + t3 - t5 - t7 - t9;
	Seta[1] = 0.2e1 * t11;
	Seta[2] = 0.0e0;
	Seta[3] = -0.2e1 * t11;
	Seta[4] = -0.1e1 + t12 + t2 - t9 + t7;
	Seta[5] = 0.2e1 - t14 - t3 - t15 + t16 + t5;
	Seta[6] = t18;
	Seta[7] = 0.0e0;
	Seta[8] = -t18;
	Seta[9] = -0.1e1 - t1 + t2 + t3 + t5 - t7 + t9;
	Seta[10] = 0.2e1 * t20;
	Seta[11] = 0.0e0;
	Seta[12] = -0.2e1 * t20;
	Seta[13] = -0.1e1 - t12 + t2 + t9 + t7;
	Seta[14] = 0.2e1 + t14 - t3 - t15 - t16 - t5;

//Szeta ===========================================
	t1 = 0.3e1 / 0.2e1 * eta;
	t2 = 0.3e1 / 0.2e1 * xi;
	t3 = xi*xi;
	t4 = eta * eta;
	t5 = zeta * xi;
	t6 = eta * zeta;
	t7 = xi * eta;
	t8 = 0.2e1 * t7;
	t10 = -xi + t3 + t7;
	t11 = xi / 0.2e1;
	t13 = eta / 0.2e1;
	t15 = -eta + t4 + t7;
	Szeta[0] = -0.1e1 / 0.2e1 + t1 + t2 + zeta - t3 - t4 - t5 - t6 - t8;
	Szeta[1] = 0.2e1 * t10;
	Szeta[2] = t11 - t3 + t5;
	Szeta[3] = -t8;
	Szeta[4] = t13 - t4 + t6;
	Szeta[5] = 0.2e1 * t15;
	Szeta[6] = -0.2e1 * zeta + 0.2e1 * t5 + 0.2e1 * t6;
	Szeta[7] = -0.2e1 * t5;
	Szeta[8] = -0.2e1 * t6;
	Szeta[9] = 0.1e1 / 0.2e1 - t1 - t2 + zeta + t3 + t4 - t5 - t6 + t8;
	Szeta[10] = -0.2e1 * t10;
	Szeta[11] = -t11 + t3 + t5;
	Szeta[12] = t8;
	Szeta[13] = -t13 + t4 + t6;
	Szeta[14] = -0.2e1 * t15;

}
//=========================================================

//       New format
//=========================================================
#include <stdlib.h>

void deriv3D(int numberOfInterp,
              double xi,  double eta, double zeta,				  
				  double * Sxi, double * Seta, double * Szeta)
{
//BETA_OUT<<"derivHex"<<xi<<eta<<zeta<<endl;

switch(numberOfInterp)                              
   {
    case 8  :derivHex8(Sxi,Seta,Szeta, xi,eta,zeta);             
              break;

    case 20 : derivHex20(Sxi,Seta,Szeta, xi,eta,zeta);             
              break;
    
	case 10 : derivTet10(Sxi,Seta,Szeta, xi,eta,zeta);
		      break;

	case 15 : derivWedge15(Sxi,Seta,Szeta, xi,eta,zeta);
		      break;

    default :BETA_OUT<<"Cannot find required interpolation functions"<<endl;
             BETA_OUT<<"This element has "<<numberOfInterp
                 <<" nodes"<<endl;    
             exit(1); 
   }

}

//=========================================================
void  shape3D(int numberOfInterp,
              double xi,  double eta, double zeta,				  
				  double * S)
{
switch(numberOfInterp)                              
   {
    case 8  : shapeHex8(S, xi,eta,zeta);             
              break;

    case 20 : shapeHex20(S, xi,eta,zeta);            
              break;
    
	case 10 : shapeTet10(S, xi,eta,zeta);
		      break;

	case 15 : shapeWedge15(S, xi,eta,zeta);
		      break;

    default :BETA_OUT<<"Cannot find required interpolation functions"<<endl;
             BETA_OUT<<"This element has "<<numberOfInterp
                 <<" nodes"<<endl;    
             exit(1); 
   }

}
//=========================================================
void deriv2D(int numberOfInterp,
              double xi,  double eta,				  
				  double * Sxi, double * Seta)
{

switch(numberOfInterp)                              
   {
    case 3  : derivTri3(Sxi,Seta, xi,eta);             
              break;

    case 4  : derivQuad4(Sxi,Seta, xi,eta);             
              break;

    case 6  : derivTri6(Sxi,Seta, xi,eta);             
              break;

    case 8  : derivQuad8(Sxi,Seta,xi,eta);             
              break;

    default :BETA_OUT<<"Cannot find required interpolation functions"<<endl;
             BETA_OUT<<"This element has "<<numberOfInterp
                 <<" nodes"<<endl;    
             exit(1); 
   }

}

//=========================================================
void  shape2D(int numberOfInterp,
              double xi,  double eta,				  
				  double * S)
{
switch(numberOfInterp)                              
   {
    case 3  : shapeTri3(S, xi,eta);             
              break;

    case 4  : shapeQuad4(S, xi,eta);             
              break;

    case 6  : shapeTri6(S, xi,eta);             
              break;

    case 8  : shapeQuad8(S, xi,eta);            
              break;

    default :BETA_OUT<<"Cannot find required interpolation functions"<<endl;
             BETA_OUT<<"This element has "<<numberOfInterp
                 <<" nodes"<<endl;    
             exit(1); 
   }

}

//=========================================================
//=========================================================
void deriv1D(int numberOfInterp,double xi,double * Sxi)
{

switch(numberOfInterp)                              
   {
    case 2  : derivLine2(Sxi,xi);             
              break;
    case 3  : derivLine3(Sxi,xi);             
              break;
    case 4  : derivLine4(Sxi,xi);             
              break;

    default :BETA_OUT<<"Cannot find required interpolation functions"<<endl;
             BETA_OUT<<"This element has "<<numberOfInterp
                 <<" nodes"<<endl;    
             exit(1); 
   }

}

//=========================================================
void  shape1D(int numberOfInterp,double xi,double * S)
{
switch(numberOfInterp)                              
   {
    case 2  : shapeLine2(S,xi);             
              break;
    case 3  : shapeLine3(S,xi);             
              break;
    case 4  : shapeLine4(S,xi);             
              break;

    default :BETA_OUT<<"Cannot find required interpolation functions"<<endl;
             BETA_OUT<<"This element has "<<numberOfInterp
                 <<" nodes"<<endl;    
             exit(1); 
   }

}
