#include "stdafx.h"

#include "PointLoad.hpp"
#include "BCs/SurfaceLoad.hpp"

#include <cmath>

void shapeQuad8( double * S, double xi, double eta);
void derivQuad8(double * Szeta,	double * Seta, double zeta, double eta);
void shapeQuad4( double * S, double xi, double eta);
void derivQuad4(double * Szeta,	double * Seta, double zeta, double eta);
void shapeTri6( double * S, double xi, double eta);
void derivTri6(double * Sxi, double * Seta, double xi, double eta);

//==============================================================================================+
void SurfaceLoad::applyTo(Equations &equations, NodeGroup *node, int numNodes, BasicModel *model)
{                         
	if (!haveReadMesh) 
	{
		cout << "\nWill read surface mesh soon for the first time." << endl;
		surfaceMesh.setFileManager(filemanager);
		surfaceMesh.ReadMesh(filename);
		haveReadMesh = true;
	}
	
	int maxSurfaceNodeNum = surfaceMesh.numNodes;
	int maxSurfaceElemNum = surfaceMesh.numElements;
	int dimensionOfMesh = 2;  // Surface mesh, do not change from two
	int numIntePt;
	
	double xi[9]; // Gauss point coordinates
	double eta[9];
	double wt[9]; // Gauss point weights
	
	double elementalNodalFlux[8];  // Value of the applied flux on the surface
	double elementalNodalForce[8]; // Computed nodal force after integration of the flux
	
	double gaussFlux; // value of the flux at the gauss point
	double modG, g1, g2, g3;      // surface scalling for a element
	
	int nodalConectivity[8];
	
	double shape[8];
	double dSdxi[8];
	double dSdeta[8];
	
	// x, y and z coordinates of nodes of a element
	double xCordElemNode[8]; 
	double yCordElemNode[8];
	double zCordElemNode[8];
	
	// Loop over surface mesh starts
	for (int ie = 0; ie < maxSurfaceElemNum; ie++)  
	{
		//cout <<"element "<<ie<<"Num node"<<(surfaceMesh.element[ie]).node.getCapacity()<<endl<<endl;
		const int numNodesPerElem = (surfaceMesh.element[ie]).node.getCapacity();
		
		if ( (numNodesPerElem == 8) || (numNodesPerElem == 4) )
		{
			numIntePt = 9;

			xi[0]= 0.774596669241483;
			xi[1]= 0.;
			xi[2]=-0.774596669241483;
			xi[3]= 0.774596669241483;
			xi[4]= 0.;
			xi[5]=-0.774596669241483;
			xi[6]= 0.774596669241483;
			xi[7]= 0.;
			xi[8]=-0.774596669241483;

			eta[0]= 0.774596669241483;
			eta[1]= 0.774596669241483;
			eta[2]= 0.774596669241483;
			eta[3]= 0.;
			eta[4]= 0.;
			eta[5]= 0.;
			eta[6]=-0.774596669241483;
			eta[7]=-0.774596669241483;
			eta[8]=-0.774596669241483;
    
			wt[0]=25./81.;
			wt[1]=40./81.;
			wt[2]=25./81.;
			wt[3]=40./81.;
			wt[4]=64./81.;
			wt[5]=40./81.;
			wt[6]=25./81.;
			wt[7]=40./81.;
			wt[8]=25./81.;
		}
		else if ( numNodesPerElem == 6 )
		{
			numIntePt = 7;

		    xi[0]     = 1.0/3;
		    eta[0]    = 1.0/3;
		    wt[0] = 0.225*0.5;


		    double val1 = 0.101286507323;
		    double val2 = 0.797426985353;
		    double weightVal = 0.125939180544 * 0.5;

		   xi[1]     = val1;
		   eta[1]    = val1;
		   wt[1] = weightVal;

		   xi[2]     = val1;
		   eta[2]    = val2;
		   wt[2] = weightVal;

		   xi[3]     = val2;
		   eta[3]    = val1;
		   wt[3] = weightVal;
		   
		   
		   val1 = 0.470142064105;
		   val2 = 0.059715871789;
		   weightVal = 0.132394152788 * 0.5;

		   xi[4]     = val1;
		   eta[4]    = val1;
		   wt[4] = weightVal;

		   xi[5]     = val1;
		   eta[5]    = val2;
		   wt[5] = weightVal;

		   xi[6]     = val2;
		   eta[6]    = val1;
		   wt[6] = weightVal;
		}
		else
		{
			cout << "\n\n Number of elements per node is wrong. Surface ele #"<<ie<<endl;
		}

		// Extract nodal connectivity and nodal coordinates 
		// of the nodes of this element
		for (int i = 0; i < numNodesPerElem; i++)
		{
			nodalConectivity[i] = (surfaceMesh.element[ie]).node[i].nodeNum;
			xCordElemNode[i] = (surfaceMesh.element[ie]).node[i].x;
			yCordElemNode[i] = (surfaceMesh.element[ie]).node[i].y;
			zCordElemNode[i] = (surfaceMesh.element[ie]).node[i].z;
			
			// read in flux value at nodes from the shared memory
			elementalNodalFlux[i] = 4*10000*0.01/(4*500*0.2); // (d3dot*w/4) = volumatric heat generation, ro*c=500*0.2
			//elementalNodalFlux[i] = 1; // for area cal
			
			//Initialize elemental nodal force to zero
			elementalNodalForce[i] = 0;
		}		
		
		//Loop over gauss point starts
		for (int m = 0; m < numIntePt; m++)
		{
			if ( numNodesPerElem == 8 )
			{ 
				shapeQuad8(shape,xi[m],eta[m]);
				derivQuad8(dSdxi,dSdeta,xi[m],eta[m]);
			}
			if ( numNodesPerElem == 4 )
			{ 
				shapeQuad4(shape,xi[m],eta[m]);
				derivQuad4(dSdxi,dSdeta,xi[m],eta[m]);
			}
			else if ( numNodesPerElem == 6 )
			{
				shapeTri6(shape,xi[m],eta[m]);
				derivTri6(dSdxi,dSdeta,xi[m],eta[m]);
			}
			else
			{
				cout << "\n\n Number of elements per node is wrong. Surface ele #"<<ie<<endl;
			}
			
			// Calculate modG for scalling the element area
			double dxdxi = 0;
			double dydxi = 0;
			double dzdxi = 0;
			
			double dxdeta = 0;
			double dydeta = 0;
			double dzdeta = 0;
			
			for (int i = 0; i < numNodesPerElem; i ++)
			{
				dxdxi = dxdxi + dSdxi[i]*xCordElemNode[i];
				dydxi = dydxi + dSdxi[i]*yCordElemNode[i];
				dzdxi = dzdxi + dSdxi[i]*zCordElemNode[i];
			
				dxdeta = dxdeta + dSdeta[i]*xCordElemNode[i];
				dydeta = dydeta + dSdeta[i]*yCordElemNode[i];
				dzdeta = dzdeta + dSdeta[i]*zCordElemNode[i];			
			}
				
			g1 = dydxi*dzdeta - dydeta*dzdxi;	
			g2 = dzdxi*dxdeta - dzdeta*dxdxi;
			g3 = dxdxi*dydeta - dxdeta*dydxi;
			
			//ref: http://www.softeng.rl.ac.uk/st/projects/felib3/Docs/PDF/Intro.pdf
			// Page 70 of above reference has the expression for surface mapping
			modG = sqrt(g1*g1 + g2*g2 + g3*g3);
			
			// Calculate the flux at the gauss point
			gaussFlux = 0.0;
            for(int i=0; i < numNodesPerElem; i++)
			{
               gaussFlux=gaussFlux+shape[i]*elementalNodalFlux[i];
			}
			
			// Calculate the elemental nodel force
            for(int i=0; i < numNodesPerElem; i++)
			{
                elementalNodalForce[i]=elementalNodalForce[i]+wt[m]*gaussFlux*shape[i]*modG;
			}
		} // Loop over gauss point ends
			
		// Apply the elemental nodal force of this element to the equation 
		for (int i = 0; i < numNodesPerElem; i ++)
		{
			int direction = 1;
			int globalDof = (*node)[nodalConectivity[i]].getFirstDof()+direction-1;
			equations.addLoadToEqn(elementalNodalForce[i],globalDof);
		}

		//double area = 0;
		//for (int i = 0; i < numNodesPerElem; i++)
		//{
		//	area = area + elementalNodalForce[i];
		//}
		//cout << "\n\n\n******************\n area of element: "<<ie<<" :" << area<<endl<<endl;
		
		
	} // Loop over surface element ends

}
//==========================================================
bool SurfaceLoad::read(istream * inStream)
{
	// did not work BETA_OUT << "\n\n We are in read load" <<endl;
	haveReadMesh = false;

	string token;
	*inStream >> token;
	if(COMPARE(token.c_str(),"end") == 0){
		return false;
	} else {
		strcpy(filename,token.c_str());
	}

	print(BETA_OUT);
	return true;
}
//==============================================================================================
void SurfaceLoad::print(ostream &ostrm)
{
	ostrm << "Surface Load read: ";
	return;
}