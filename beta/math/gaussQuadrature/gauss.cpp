#include "stdafx.h"

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "gauss.hpp"
using namespace std;

double *GAUSS::coords[MAX_Gauss_ORDER]={0};
double *GAUSS::weights[MAX_Gauss_ORDER]={0};
bool GAUSS::calculated[MAX_Gauss_ORDER]={false};

GAUSS::GAUSS(int k)
{  
	double EPS=3e-15,z,z1,P1,P2,P3,PP,xl,xm;
	int i,j,m;

	//cout<<"Called gauss: "<<k<<endl;

	currentOrder=k-1;
	currentCoords=coords[currentOrder];
	currentWeights=weights[currentOrder];

	if(k > MAX_Gauss_ORDER || k < 0) 
	{
		printf("\nError: GAUSS::GAUSS - Order %d requested, maximum is %d\n\n",k,MAX_Gauss_ORDER);
		exit(1);
	}
	if(calculated[currentOrder]) return;
	cout << "Calculating Gauss quadrature order " << (currentOrder+1) << endl;
	currentCoords = coords[currentOrder] = new double[k];
	currentWeights = weights[currentOrder] = new double[k];
	m=(k+1)/2;
	xm=0;
	xl=1;
	for(i=0;i<m;i++)
	{
//		cout << i << endl;
		z=cos(M_PI*((double)i+1-.25)/((double)k+.5));
L1:
		P1=1;
		P2=0;
		for(j=0;j<k;j++)
		{
			P3=P2;
			P2=P1;
			P1=((2*j+1)*z*P2-(j)*P3)/(j+1);
		}
		PP=k*(z*P1-P2)/(z*z-1);
		z1=z;
		z=z1-P1/PP;
		if(fabs(z-z1)>EPS) goto L1;
		coords[currentOrder][i]=xm-xl*z;
		coords[currentOrder][k-i-1]=xm+xl*z;
		weights[currentOrder][i]=2*xl/((1-z*z)*PP*PP);
		weights[currentOrder][k-i-1]=weights[currentOrder][i];
	}
	calculated[currentOrder] = true;
}
void GAUSS::deleteAll()
{
	int i;
	for(i=0;i<MAX_Gauss_ORDER;i++) {
		if(coords[i]) delete [] coords[i];coords[i]=0;
		if(weights[i]) delete [] weights[i];weights[i]=0;
		calculated[i]=false;
	}
}
//==============================================================
//==============================================================

//==============================================================
//==============================================================

void getQuadraturePointsForHex( int integrationOrder,
									  double *xi, double *eta, double *zeta,
									  double *weight )
{
	GAUSS gauss(integrationOrder);
	double xiPoint, etaPoint;
	double weightXi, weightEta;
	int ip = -1;

	for(int i=0;i<integrationOrder;i++)
		{
		  xiPoint = gauss.getCoord(i);
		  weightXi= gauss.getWeight(i);
		for(int j=0;j<integrationOrder;j++)
		  {
			etaPoint = gauss.getCoord(j);
			weightEta= gauss.getWeight(j);

			for(int k=0;k<integrationOrder;k++)
				{
				ip++;
					 xi[ip] = xiPoint;
					eta[ip] = etaPoint;
				  zeta[ip] = gauss.getCoord(k);
				weight[ip] = weightXi * weightEta * gauss.getWeight(k);
				}
			 }
		}

 //cout<<"Integration points and weights "<<endl;
 //for(ip= 0; ip<integrationOrder*integrationOrder*integrationOrder; ip++)
 //{cout<<xi[ip]<<"  "<<eta[ip]<<"  "<<zeta[ip]<<"  "<<weight[ip]<<endl;  }

 }//end of  getQuadraturePointsForHex

//..............................................................

//===========================================================================
void getQuadraturePointsForQuad( int integrationOrder,
									  double *xi, double *eta, double *weight )
{
	GAUSS gauss(integrationOrder);
	double xiPoint;
	double weightXi;
	int ip = -1;

	for(int i=0;i<integrationOrder;i++)
		{
		  xiPoint = gauss.getCoord(i);
		  weightXi= gauss.getWeight(i);

			for(int k=0;k<integrationOrder;k++)
				{
				ip++;
					 xi[ip] = xiPoint;
				   eta[ip] = gauss.getCoord(k);
				weight[ip] = weightXi * gauss.getWeight(k);
				}
		}

// cout<<"Integration points and weights "<<endl;
// for(ip= 0; ip<integrationOrder*integrationOrder; ip++)
// {cout<<xi[ip]<<"  "<<eta[ip]<<"  "<<weight[ip]<<endl;  }

 }//end of  getQuadraturePointsForQuad

//..............................................................
//==============================================================
#include "triQuadPoint.hpp"

void TriQuadPoints::getQuadraturePointsForTri( int integrationOrder,
									  double **p_xi, double **p_eta,
									  double **p_weight )
{

   if(xi[integrationOrder]==NULL){
	   int length = 0;
	   switch (integrationOrder){
		   case 1:
			   length = 1;
			   break;
		   case 2:
			   length = 3;
			   break;
		   case 3:
			   length = 4;
			   break;
		   case 4:
			   length = 7;
			   break;
		   default:
			   cout << "Integration Order " << integrationOrder << " is not supported for triangles." << endl;
			   exit(1);
	   } 
	   xi[integrationOrder] = new double [ length ];
	   eta[integrationOrder] = new double [ length ];
	   weight[integrationOrder] = new double [ length ];
	   cout<<"Initialized order:TriQuadPoints::getQuadraturePointsForTri "
		     <<integrationOrder<<endl;
   }
	else
	{
 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_eta)    =    eta[integrationOrder];
	 (*p_weight) = weight[integrationOrder];
	 return;
	}

   switch(integrationOrder){
	   double val1, val2, weightVal;
	   case 1:
		   xi[integrationOrder][0]     = 1.0/3;
		   eta[integrationOrder][0]    = 1.0/3;
		   weight[integrationOrder][0] = 0.5;
           break;

	   case 2:
		   val1 = 0.5;
		   val2 = 0.0;
		   weightVal = 1.0/6;

		   xi[integrationOrder][0]     = val1;
		   eta[integrationOrder][0]    = val1;
		   weight[integrationOrder][0] = weightVal;

		   xi[integrationOrder][1]     = val1;
		   eta[integrationOrder][1]    = val2;
		   weight[integrationOrder][1] = weightVal;

		   xi[integrationOrder][2]     = val2;
		   eta[integrationOrder][2]    = val1;
		   weight[integrationOrder][2] = weightVal;
		   break;

	   case 3:
		   val1 = 0.2;
		   val2 = 0.6;
		   weightVal = 25.0/96;

		   xi[integrationOrder][0]     = val1;
		   eta[integrationOrder][0]    = val1;
		   weight[integrationOrder][0] = weightVal;

		   xi[integrationOrder][1]     = val1;
		   eta[integrationOrder][1]    = val2;
		   weight[integrationOrder][1] = weightVal;

		   xi[integrationOrder][2]     = val2;
		   eta[integrationOrder][2]    = val1;
		   weight[integrationOrder][2] = weightVal;
		   
		   
		   xi[integrationOrder][3]     = 1.0/3;
		   eta[integrationOrder][3]    = 1.0/3;
		   weight[integrationOrder][3] = -27.0/96;
		   break;

	   case 4:
           xi[integrationOrder][0]     = 1.0/3;
		   eta[integrationOrder][0]    = 1.0/3;
		   weight[integrationOrder][0] = 0.225*0.5;


		   val1 = 0.101286507323;
		   val2 = 0.797426985353;
		   weightVal = 0.125939180544 * 0.5;

		   xi[integrationOrder][1]     = val1;
		   eta[integrationOrder][1]    = val1;
		   weight[integrationOrder][1] = weightVal;

		   xi[integrationOrder][2]     = val1;
		   eta[integrationOrder][2]    = val2;
		   weight[integrationOrder][2] = weightVal;

		   xi[integrationOrder][3]     = val2;
		   eta[integrationOrder][3]    = val1;
		   weight[integrationOrder][3] = weightVal;
		   
		   
		   val1 = 0.470142064105;
		   val2 = 0.059715871789;
		   weightVal = 0.132394152788 * 0.5;

		   xi[integrationOrder][4]     = val1;
		   eta[integrationOrder][4]    = val1;
		   weight[integrationOrder][4] = weightVal;

		   xi[integrationOrder][5]     = val1;
		   eta[integrationOrder][5]    = val2;
		   weight[integrationOrder][5] = weightVal;

		   xi[integrationOrder][6]     = val2;
		   eta[integrationOrder][6]    = val1;
		   weight[integrationOrder][6] = weightVal;
		   
		   
		   break;
	}

 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_eta)    =    eta[integrationOrder];
	 (*p_weight) = weight[integrationOrder];

	 return;
 }//end of  getQuadraturePointsForTri
//..............................................................
void TriQuadPoints::deleteAll()
{
	int i;
	for(i=0;i<5;i++) {
		if(xi[i]){
			delete [] xi[i];xi[i]=0;
			delete [] eta[i];eta[i]=0;
			delete [] weight[i];weight[i]=0;
		}
	}
}
//..............................................................
// Static data
double * TriQuadPoints::xi[5];
double * TriQuadPoints::eta[5];
double * TriQuadPoints::weight[5];

//===============================================================
#include "hexQuadPoint.hpp"

void HexQuadPoints::getQuadraturePointsForHex( int integrationOrder,
									  double **p_xi, double **p_eta, 
									  double **p_zeta,
									  double **p_weight )
{
	GAUSS gauss(integrationOrder);
	double xiPoint, etaPoint;
	double weightXi, weightEta;
	int ip = -1;

   if(xi[integrationOrder]==NULL){
	   int length =integrationOrder*integrationOrder*integrationOrder; 
	   xi[integrationOrder] = new double [ length ];
	   eta[integrationOrder] = new double [ length ];
	   zeta[integrationOrder] = new double [ length ];
	   weight[integrationOrder] = new double [ length ];
	   cout<<"Initialized order:HexQuadPoints::getQuadraturePointsForHex "
		     <<integrationOrder<<endl;
   }
	else
	{
 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_eta)    =    eta[integrationOrder];
	 (*p_zeta)   =   zeta[integrationOrder];
	 (*p_weight) = weight[integrationOrder];
	 return;
	}

	for(int i=0;i<integrationOrder;i++)
		{
		  xiPoint = gauss.getCoord(i);
		  weightXi= gauss.getWeight(i);
		for(int j=0;j<integrationOrder;j++)
		  {
			etaPoint = gauss.getCoord(j);
			weightEta= gauss.getWeight(j);

			for(int k=0;k<integrationOrder;k++)
				{
				ip++;
					 xi[integrationOrder][ip] = xiPoint;
					eta[integrationOrder][ip] = etaPoint;
				  zeta[integrationOrder][ip] = gauss.getCoord(k);
				weight[integrationOrder][ip] = weightXi * weightEta * gauss.getWeight(k);
				}
			 }
		}
 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_eta)    =    eta[integrationOrder];
	 (*p_zeta)   =   zeta[integrationOrder];
	 (*p_weight) = weight[integrationOrder];

/*
	 cout<<"Integration points and weights:  "
		  <<"getQuadraturePointsForHex"<<endl;
 for(ip= 0; ip<integrationOrder*integrationOrder*integrationOrder; ip++)
 {cout<<(*p_xi)[ip]<<"  "<<(*p_eta)[ip]<<"  "
      <<(*p_zeta)[ip]<<"  "<<(*p_weight)[ip]<<endl;  }
 //exit(1);
 */

	 return;
 }//end of  getQuadraturePointsForHex
//..............................................................
void HexQuadPoints::deleteAll()
{
	int i;
	for(i=0;i<5;i++) {
		if(xi[i]){
			delete [] xi[i];xi[i]=0;
			delete [] eta[i];eta[i]=0;
			delete [] zeta[i];zeta[i]=0;
			delete [] weight[i];weight[i]=0;
		}
	}
}
//..............................................................
// Static data
double * HexQuadPoints::xi[5];
double * HexQuadPoints::eta[5];
double * HexQuadPoints::zeta[5];
double * HexQuadPoints::weight[5];

//===============================================================
//==============================================================

#include "tetQuadPoint.hpp"

void TetQuadPoints::getQuadraturePointsForTet( int integrationOrder,
									  double **p_xi, double **p_eta, 
									  double **p_zeta,
									  double **p_weight )
{

   if(xi[integrationOrder]==NULL){
	   int length = 0;
	   switch (integrationOrder){
		   case 1:
			   length = 1;
			   break;
		   case 2:
			   length = 4;
			   break;
		   case 3:
			   length = 15;
			   break;
		   default:
			   cout << "Integration Order " << integrationOrder << " is not supported for tets." << endl;
			   exit(1);
	   } 
	   xi[integrationOrder] = new double [ length ];
	   eta[integrationOrder] = new double [ length ];
	   zeta[integrationOrder] = new double [ length ];
	   weight[integrationOrder] = new double [ length ];
	   cout<<"Initialized order:TetQuadPoints::getQuadraturePointsForTet "
		     <<integrationOrder<<endl;
   }
	else
	{
 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_eta)    =    eta[integrationOrder];
	 (*p_zeta)   =   zeta[integrationOrder];
	 (*p_weight) = weight[integrationOrder];
	 return;
	}

   switch(integrationOrder){
	   double val1, val2, weightVal;
	   case 1:
		   xi[integrationOrder][0]     = 0.25;
		   eta[integrationOrder][0]    = 0.25;
		   zeta[integrationOrder][0]   = 0.25;
		   weight[integrationOrder][0] = 1.0/6;
		   break;
	   case 2:
		   val1 = (5.0-sqrt(5.0))/20;
		   val2 = (5.0+3*sqrt(5.0))/20;
		   xi[integrationOrder][0]     = val1;
		   eta[integrationOrder][0]    = val1;
		   zeta[integrationOrder][0]   = val1;
		   weight[integrationOrder][0] = 1.0/24;

		   xi[integrationOrder][1]     = val2;
		   eta[integrationOrder][1]    = val1;
		   zeta[integrationOrder][1]   = val1;
		   weight[integrationOrder][1] = 1.0/24;

		   xi[integrationOrder][2]     = val1;
		   eta[integrationOrder][2]    = val2;
		   zeta[integrationOrder][2]   = val1;
		   weight[integrationOrder][2] = 1.0/24;

		   xi[integrationOrder][3]     = val1;
		   eta[integrationOrder][3]    = val1;
		   zeta[integrationOrder][3]   = val2;
		   weight[integrationOrder][3] = 1.0/24;
		   break;
	   case 3:
		   xi[integrationOrder][0]     = 0.25;
		   eta[integrationOrder][0]    = 0.25;
		   zeta[integrationOrder][0]   = 0.25;
		   weight[integrationOrder][0] = 16.0/810;

		   val1 = (7.0 - sqrt(15.0))/34;
		   val2 = (13.0 + 3.0*sqrt(15.0))/34;
		   weightVal = (2665.0 + 14.0*sqrt(15.0))/226800;

		   xi[integrationOrder][1]     = val1;
		   eta[integrationOrder][1]    = val1;
		   zeta[integrationOrder][1]   = val1;
		   weight[integrationOrder][1] = weightVal;

		   xi[integrationOrder][2]     = val2;
		   eta[integrationOrder][2]    = val1;
		   zeta[integrationOrder][2]   = val1;
		   weight[integrationOrder][2] = weightVal;

		   xi[integrationOrder][3]     = val1;
		   eta[integrationOrder][3]    = val2;
		   zeta[integrationOrder][3]   = val1;
		   weight[integrationOrder][3] = weightVal;

		   xi[integrationOrder][4]     = val1;
		   eta[integrationOrder][4]    = val1;
		   zeta[integrationOrder][4]   = val2;
		   weight[integrationOrder][4] = weightVal;

		   val1 = (7.0 + sqrt(15.0))/34;
		   val2 = (13.0 - 3.0*sqrt(15.0))/34;
		   weightVal = (2665.0 - 14.0*sqrt(15.0))/226800;

		   xi[integrationOrder][5]     = val1;
		   eta[integrationOrder][5]    = val1;
		   zeta[integrationOrder][5]   = val1;
		   weight[integrationOrder][5] = weightVal;

		   xi[integrationOrder][6]     = val2;
		   eta[integrationOrder][6]    = val1;
		   zeta[integrationOrder][6]   = val1;
		   weight[integrationOrder][6] = weightVal;

		   xi[integrationOrder][7]     = val1;
		   eta[integrationOrder][7]    = val2;
		   zeta[integrationOrder][7]   = val1;
		   weight[integrationOrder][7] = weightVal;

		   xi[integrationOrder][8]     = val1;
		   eta[integrationOrder][8]    = val1;
		   zeta[integrationOrder][8]   = val2;
		   weight[integrationOrder][8] = weightVal;

		   val1 = (10.0 - 2.0*sqrt(15.0))/40;
		   val2 = (10.0 + 2.0*sqrt(15.0))/40;
		   weightVal = 20.0/2268;

		   xi[integrationOrder][9]     = val2;
		   eta[integrationOrder][9]    = val1;
		   zeta[integrationOrder][9]   = val1;
		   weight[integrationOrder][9] = weightVal;

		   xi[integrationOrder][10]     = val1;
		   eta[integrationOrder][10]    = val2;
		   zeta[integrationOrder][10]   = val1;
		   weight[integrationOrder][10] = weightVal;

		   xi[integrationOrder][11]     = val1;
		   eta[integrationOrder][11]    = val1;
		   zeta[integrationOrder][11]   = val2;
		   weight[integrationOrder][11] = weightVal;

		   xi[integrationOrder][12]     = val1;
		   eta[integrationOrder][12]    = val2;
		   zeta[integrationOrder][12]   = val2;
		   weight[integrationOrder][12] = weightVal;

		   xi[integrationOrder][13]     = val2;
		   eta[integrationOrder][13]    = val1;
		   zeta[integrationOrder][13]   = val2;
		   weight[integrationOrder][13] = weightVal;

		   xi[integrationOrder][14]     = val2;
		   eta[integrationOrder][14]    = val2;
		   zeta[integrationOrder][14]   = val1;
		   weight[integrationOrder][14] = weightVal;
		   break;
	}

 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_eta)    =    eta[integrationOrder];
	 (*p_zeta)   =   zeta[integrationOrder];
	 (*p_weight) = weight[integrationOrder];

	 return;
 }//end of  getQuadraturePointsForTet
//..............................................................
void TetQuadPoints::deleteAll()
{
	int i;
	for(i=0;i<5;i++) {
		if(xi[i]){
			delete [] xi[i];xi[i]=0;
			delete [] eta[i];eta[i]=0;
			delete [] zeta[i];zeta[i]=0;
			delete [] weight[i];weight[i]=0;
		}
	}
}
//..............................................................
// Static data
double * TetQuadPoints::xi[5];
double * TetQuadPoints::eta[5];
double * TetQuadPoints::zeta[5];
double * TetQuadPoints::weight[5];

//===============================================================
#include "wedgeQuadPoint.hpp"

void WedgeQuadPoints::getQuadraturePointsForWedge( int integrationOrder,
									  double **p_xi, double **p_eta, 
									  double **p_zeta,
									  double **p_weight )
{

   if(xi[integrationOrder]==NULL){
	   int length = 0;
	   switch (integrationOrder){
		   case 1:
			   length = 2;
			   break;
		   case 2:
			   length = 9;
			   break;
		   case 3:
			   length = 18;
			   break;
		   default:
			   cout << "Integration Order " << integrationOrder << " is not supported for wedges." << endl;
			   exit(1);
	   } 
	   xi[integrationOrder] = new double [ length ];
	   eta[integrationOrder] = new double [ length ];
	   zeta[integrationOrder] = new double [ length ];
	   weight[integrationOrder] = new double [ length ];
	   cout<<"Initialized order:WedgeQuadPoints::getQuadraturePointsForWedge "
		     <<integrationOrder<<endl;
   }
	else
	{
 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_eta)    =    eta[integrationOrder];
	 (*p_zeta)   =   zeta[integrationOrder];
	 (*p_weight) = weight[integrationOrder];
	 return;
	}

   switch(integrationOrder){
	   double val1, val2, val3, weightVal;
	   case 1:
		   xi[integrationOrder][0]     = 1.0/3;
		   eta[integrationOrder][0]    = 1.0/3;
		   zeta[integrationOrder][0]   = 1.0/sqrt(3.0);
		   weight[integrationOrder][0] = 1.0/2;

		   xi[integrationOrder][1]     = 1.0/3;
		   eta[integrationOrder][1]    = 1.0/3;
		   zeta[integrationOrder][1]   = -1.0/sqrt(3.0);
		   weight[integrationOrder][1] = 1.0/2;
		   break;
	   case 2:
		   val1 = 1.0/6;
		   val2 = 4.0/6;
		   val3 = sqrt(3.0/5);
		   weightVal = 5.0/54;

		   xi[integrationOrder][0]     = val1;
		   eta[integrationOrder][0]    = val1;
		   zeta[integrationOrder][0]   = val3;
		   weight[integrationOrder][0] = weightVal;

		   xi[integrationOrder][1]     = val1;
		   eta[integrationOrder][1]    = val2;
		   zeta[integrationOrder][1]   = val3;
		   weight[integrationOrder][1] = weightVal;

		   xi[integrationOrder][2]     = val2;
		   eta[integrationOrder][2]    = val1;
		   zeta[integrationOrder][2]   = val3;
		   weight[integrationOrder][2] = weightVal;

		   xi[integrationOrder][3]     = val1;
		   eta[integrationOrder][3]    = val1;
		   zeta[integrationOrder][3]   = -val3;
		   weight[integrationOrder][3] = weightVal;

		   xi[integrationOrder][4]     = val1;
		   eta[integrationOrder][4]    = val2;
		   zeta[integrationOrder][4]   = -val3;
		   weight[integrationOrder][4] = weightVal;

		   xi[integrationOrder][5]     = val2;
		   eta[integrationOrder][5]    = val1;
		   zeta[integrationOrder][5]   = -val3;
		   weight[integrationOrder][5] = weightVal;
		   
		   val3 = 0.0;
		   weightVal = 8.0/54;

		   xi[integrationOrder][6]     = val1;
		   eta[integrationOrder][6]    = val1;
		   zeta[integrationOrder][6]   = val3;
		   weight[integrationOrder][6] = weightVal;

		   xi[integrationOrder][7]     = val1;
		   eta[integrationOrder][7]    = val2;
		   zeta[integrationOrder][7]   = val3;
		   weight[integrationOrder][7] = weightVal;

		   xi[integrationOrder][8]     = val2;
		   eta[integrationOrder][8]    = val1;
		   zeta[integrationOrder][8]   = val3;
		   weight[integrationOrder][8] = weightVal;
		   break;
	   case 3:
		   val1 = 1.0/6;
		   val2 = 4.0/6;
		   val3 = sqrt(3.0/5);
		   weightVal = 1.0/12;

		   xi[integrationOrder][0]     = val1;
		   eta[integrationOrder][0]    = val1;
		   zeta[integrationOrder][0]   = val3;
		   weight[integrationOrder][0] = weightVal;

		   xi[integrationOrder][1]     = val1;
		   eta[integrationOrder][1]    = val2;
		   zeta[integrationOrder][1]   = val3;
		   weight[integrationOrder][1] = weightVal;

		   xi[integrationOrder][2]     = val2;
		   eta[integrationOrder][2]    = val1;
		   zeta[integrationOrder][2]   = val3;
		   weight[integrationOrder][2] = weightVal;

		   xi[integrationOrder][3]     = val1;
		   eta[integrationOrder][3]    = val1;
		   zeta[integrationOrder][3]   = -val3;
		   weight[integrationOrder][3] = weightVal;

		   xi[integrationOrder][4]     = val1;
		   eta[integrationOrder][4]    = val2;
		   zeta[integrationOrder][4]   = -val3;
		   weight[integrationOrder][4] = weightVal;

		   xi[integrationOrder][5]     = val2;
		   eta[integrationOrder][5]    = val1;
		   zeta[integrationOrder][5]   = -val3;
		   weight[integrationOrder][5] = weightVal;
		   
		   val3 = 0.0;
		   weightVal = 2.0/15;

		   xi[integrationOrder][6]     = val1;
		   eta[integrationOrder][6]    = val1;
		   zeta[integrationOrder][6]   = val3;
		   weight[integrationOrder][6] = weightVal;

		   xi[integrationOrder][7]     = val1;
		   eta[integrationOrder][7]    = val2;
		   zeta[integrationOrder][7]   = val3;
		   weight[integrationOrder][7] = weightVal;

		   xi[integrationOrder][8]     = val2;
		   eta[integrationOrder][8]    = val1;
		   zeta[integrationOrder][8]   = val3;
		   weight[integrationOrder][8] = weightVal;

		   val1 = 1.0/2;
		   val2 = 0.0;
		   val3 = sqrt(3.0/5);
		   weightVal = 1.0/108;

		   xi[integrationOrder][9]     = val1;
		   eta[integrationOrder][9]    = val1;
		   zeta[integrationOrder][9]   = val3;
		   weight[integrationOrder][9] = weightVal;

		   xi[integrationOrder][10]     = val1;
		   eta[integrationOrder][10]    = val2;
		   zeta[integrationOrder][10]   = val3;
		   weight[integrationOrder][10] = weightVal;

		   xi[integrationOrder][11]     = val2;
		   eta[integrationOrder][11]    = val1;
		   zeta[integrationOrder][11]   = val3;
		   weight[integrationOrder][11] = weightVal;

		   xi[integrationOrder][12]     = val1;
		   eta[integrationOrder][12]    = val1;
		   zeta[integrationOrder][12]   = -val3;
		   weight[integrationOrder][12] = weightVal;

		   xi[integrationOrder][13]     = val1;
		   eta[integrationOrder][13]    = val2;
		   zeta[integrationOrder][13]   = -val3;
		   weight[integrationOrder][13] = weightVal;

		   xi[integrationOrder][14]     = val2;
		   eta[integrationOrder][14]    = val1;
		   zeta[integrationOrder][14]   = -val3;
		   weight[integrationOrder][14] = weightVal;
		   
		   val3 = 0.0;
		   weightVal = 2.0/135;

		   xi[integrationOrder][15]     = val1;
		   eta[integrationOrder][15]    = val1;
		   zeta[integrationOrder][15]   = val3;
		   weight[integrationOrder][15] = weightVal;

		   xi[integrationOrder][16]     = val1;
		   eta[integrationOrder][16]    = val2;
		   zeta[integrationOrder][16]   = val3;
		   weight[integrationOrder][16] = weightVal;

		   xi[integrationOrder][17]     = val2;
		   eta[integrationOrder][17]    = val1;
		   zeta[integrationOrder][17]   = val3;
		   weight[integrationOrder][17] = weightVal;
		   break;
	}

 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_eta)    =    eta[integrationOrder];
	 (*p_zeta)   =   zeta[integrationOrder];
	 (*p_weight) = weight[integrationOrder];

	 return;
 }//end of  getQuadraturePointsForWedge
//..............................................................
void WedgeQuadPoints::deleteAll()
{
	int i;
	for(i=0;i<5;i++) {
		if(xi[i]){
			delete [] xi[i];xi[i]=0;
			delete [] eta[i];eta[i]=0;
			delete [] zeta[i];zeta[i]=0;
			delete [] weight[i];weight[i]=0;
		}
	}
}
//..............................................................
// Static data
double * WedgeQuadPoints::xi[5];
double * WedgeQuadPoints::eta[5];
double * WedgeQuadPoints::zeta[5];
double * WedgeQuadPoints::weight[5];

//===============================================================
//===============================================================
#include "quadQuadPoint.hpp"

void QuadQuadPoints::getQuadraturePointsForQuad( int integrationOrder,
									  double **p_xi, double **p_eta, 
									  double **p_weight )
{
	//cout<<"entering getQuadraturePointsForQuad"<<endl;
	GAUSS gauss(integrationOrder);
	double xiPoint;
	double weightXi;
	int ip = -1;

   if(xi[integrationOrder]==NULL)
	{int length =integrationOrder*integrationOrder; 
   	  xi[integrationOrder] = new double [ length ];
		 eta[integrationOrder] = new double [ length ];
	 weight[integrationOrder] = new double [ length ];
		 cout<<"Initialized order:QuadQuadPoints::getQuadraturePointsForQuad "
		     <<integrationOrder<<endl;
   }
	else
	{
 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_eta)    =    eta[integrationOrder];
	 //(*p_zeta)   =   zeta[integrationOrder];
	 (*p_weight) = weight[integrationOrder];
	 return;
	}

	for(int i=0;i<integrationOrder;i++)
		{
		  xiPoint = gauss.getCoord(i);
		  weightXi= gauss.getWeight(i);
		
			for(int k=0;k<integrationOrder;k++)
				{
				ip++;
					 xi[integrationOrder][ip] = xiPoint;
				   eta[integrationOrder][ip] = gauss.getCoord(k);
				weight[integrationOrder][ip] = weightXi * gauss.getWeight(k);
				}// k
		}// i
 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_eta)    =    eta[integrationOrder];
	 (*p_weight) = weight[integrationOrder];

/*
	 cout<<"Integration points and weights:  "
		  <<"getQuadraturePointsForQuad"<<endl;
 for(ip= 0; ip<integrationOrder*integrationOrder; ip++)
 {cout<<(*p_xi)[ip]<<"  "<<(*p_eta)[ip]<<"  "
      <<(*p_weight)[ip]<<endl;  }
 */
	 return;
 }//end of  getQuadraturePointsForQuad
//..............................................................
void QuadQuadPoints::deleteAll()
{
	int i;
	for(i=0;i<5;i++) {
		if(xi[i]){
			delete [] xi[i];xi[i]=0;
			delete [] eta[i];eta[i]=0;
			delete [] weight[i];weight[i]=0;
		}
	}
}
//..............................................................
// Static data

double * QuadQuadPoints::xi[5];
double * QuadQuadPoints::eta[5];
//double * QuadQuadPoints::zeta[5];
double * QuadQuadPoints::weight[5];

//==============================================================
//==============================================================
//==============================================================
#include "lineQuadPoint.hpp"
void LineQuadPoints::getQuadraturePointsForLine( int integrationOrder,
									  double **p_xi, 
									  double **p_weight )
{
	//cout<<"entering getQuadraturePointsForQuad"<<endl;
	GAUSS gauss(integrationOrder);

   if(xi[integrationOrder]==NULL)
	{int length =integrationOrder; 
   	  xi[integrationOrder] = new double [ length ];
	 weight[integrationOrder] = new double [ length ];
		 cout<<"Initialized order:LineQuadPoints::getQuadraturePointsForLine "
		     <<integrationOrder<<endl;
   }
	else
	{
 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_weight) = weight[integrationOrder];
	 return;
	}

	for(int i=0;i<integrationOrder;i++)
		{
        xi[integrationOrder][i] = gauss.getCoord(i);;
		weight[integrationOrder][i] = gauss.getWeight(i) ;
		}// i
 	 (*p_xi)     =     xi[integrationOrder];
	 (*p_weight) = weight[integrationOrder];

	 return;
 }//end of  getQuadraturePointsForLine

//..............................................................
void LineQuadPoints::deleteAll()
{
	int i;
	for(i=0;i<5;i++) {
		if(xi[i]){
			delete [] xi[i];xi[i]=0;
			delete [] weight[i];weight[i]=0;
		}
	}
}
//..............................................................

// Static data

double * LineQuadPoints::xi[5];
double * LineQuadPoints::weight[5]; 

//==============================================================
void releaseGaussQuadratureStaticMemory()
{
	GAUSS::deleteAll();
	LineQuadPoints::deleteAll();
	QuadQuadPoints::deleteAll();
	HexQuadPoints::deleteAll();

}