#include "stdafx.h"

#include "models/ElasticityModel.hpp"

struct AdditionalConstraint
{
        char s[40];
        int node;
        int dof;
        double factor;
};

//const double EPSILON = 1e-4;
//const double EPSILON = 1e-6;
double MPCTransformationTolerance = 1e-6;

enum EntityType {DUMMY_NODE, POINT, LINE, PLANE}; 
//==============================================================================
void PointsInEntity(BasicMesh *mesh, double T[3][3], double dx, double dy, double dz,
							  int Dim, double Ksi, double Eta, 
							  int &numPE, int *pl) 
// *pl should be allocated before calling this function
{
	int i;
	double a,b,c,  a1, b1, c1;

	Node *n;
	numPE=0;
	for (i=0; i<mesh->numNodes; i++) {
		n=&mesh->node[i];
		a1=n->x-dx;
		b1=n->y-dy;
		c1=n->z-dz;

		a=T[0][0]*a1+T[0][1]*b1+T[0][2]*c1;  // ksi
		b=T[1][0]*a1+T[1][1]*b1+T[1][2]*c1;  // eta
		c=T[2][0]*a1+T[2][1]*b1+T[2][2]*c1;  // z

//		if ((fabs(c) < EPSILON) && (((a >=0.0) && ( a <= (Ksi+ EPSILON)))||(Dim < 0))) {  
		if ((fabs(c) < MPCTransformationTolerance) && (((a >= -MPCTransformationTolerance) && ( a <= (Ksi+ MPCTransformationTolerance)))||(Dim < 0))) {  
			// on the plane and in the region
		    if ( ((Dim==2) && (fabs(b) < MPCTransformationTolerance)) || (Dim==-2) )  { pl[numPE]=i; numPE++; } // line
//		    else if( ((Dim==3) && (b >= 0) && ( b <= (Eta+EPSILON) )) || (Dim==-3) ) { pl[numPE]=i; numPE++; } // plane
		    else if( ((Dim==3) && (b >= -MPCTransformationTolerance) && ( b <= (Eta+MPCTransformationTolerance) )) || (Dim==-3) ) { pl[numPE]=i; numPE++; } // plane
		}
	}
}
//==============================================================================
void GetMatchedList(BasicMesh *mesh,
		double Ts[3][3], double sdx, double sdy, double sdz,
		double Tm[3][3], double mdx, double mdy, double mdz,
		int Dim, int numPE, int *pS, int *pM)
// *pS and *pM should be allocated before calling this function
{
	int i, id;
	double a, b, c=0.0, a1, b1, c1;

	for(i=0; i< numPE; i++) {
		a1=mesh->node[pS[i]].x-sdx;
		b1=mesh->node[pS[i]].y-sdy;
		c1=mesh->node[pS[i]].z-sdz;
		a=Ts[0][0]*a1+Ts[0][1]*b1+Ts[0][2]*c1;

		if (abs(Dim)==2) {b=0.0; }  // line 
		else {b=Ts[1][0]*a1+Ts[1][1]*b1+Ts[1][2]*c1;} // plane

		// find x, y, z on Master entity; 
		// note that Tm is orthogonal matrix

		a1=Tm[0][0]*a+Tm[1][0]*b + mdx;  // ksi
		b1=Tm[0][1]*a+Tm[1][1]*b + mdy;  // eta
		c1=Tm[0][2]*a+Tm[1][2]*b + mdz;  // z

		id=mesh->getNodeIndex(a1,b1,c1);
		if (id >=0 ) {pM[i]=id;} 
		else {
			cerr << "No corresponding node on Master plane for node "<< pS[i];
			exit(1);
//			FatalError("No cooresponding node");
		}
	}
}
//==============================================================================
void getTransT(double T[3][3], double c[3][3], double & ksi, double & eta, int Dim)
{
	double a1,b1,c1,d1, a2,b2,c2, d2, a3, b3, c3, d3;

	if (abs(Dim) < 2) { ksi=0.0; eta=0.0; return; }
	
	switch (abs(Dim)) {
	case 2:
		eta=0.0;
		a1=c[1][0]-c[0][0];  // p1 origin, p2 end
		b1=c[1][1]-c[0][1];
		c1=c[1][2]-c[0][2];
		ksi=sqrt(a1*a1+b1*b1+c1*c1);
		if (fabs(a1) < MPCTransformationTolerance) {a2=1.0;} else { a2=0.0;}
		if (fabs(b1) < MPCTransformationTolerance) {b2=1.0;} else { b2=0.0;}
		if (fabs(c1) < MPCTransformationTolerance) {c2=1.0;} else { c2=0.0;}
		if (fabs(a2)+fabs(b2)+fabs(c2) < MPCTransformationTolerance) {a2=1.0;}
		a3=b1*c2-b2*c1;
		b3=c1*a2-a1*c2;
		c3=a1*b2-a2*b1;
		a2=b1*c3-b3*c1;
		b2=c1*a3-a1*c3;
		c2=a1*b3-a3*b1;
		d1=ksi;
		d2=sqrt(a2*a2+b2*b2+c2*c2);
		d3=sqrt(a3*a3+b3*b3+c3*c3);
        T[0][0]=a1/d1; T[1][0]=a2/d2; T[2][0]=a3/d3;
		T[0][1]=b1/d1; T[1][1]=b2/d2; T[2][1]=b3/d3;
		T[0][2]=c1/d1; T[1][2]=c2/d2; T[2][2]=c3/d3;
				break;
	case 3:
		a1=c[0][0]-c[1][0];  // p2 origin, p2p1 -> x; p2p3 -> y
		b1=c[0][1]-c[1][1];
		c1=c[0][2]-c[1][2];
		ksi=sqrt(a1*a1+b1*b1+c1*c1);
		
		a2=c[2][0]-c[1][0];
		b2=c[2][1]-c[1][1];
		c2=c[2][2]-c[1][2];
		eta=sqrt(a2*a2+b2*b2+c2*c2);
		
		a3=b1*c2-b2*c1;
		b3=c1*a2-a1*c2;
		c3=a1*b2-a2*b1;

		d1=ksi;
		d2=eta;
		d3=sqrt(a3*a3+b3*b3+c3*c3);
		T[0][0]=a1/d1; T[1][0]=a2/d2; T[2][0]=a3/d3;
		T[0][1]=b1/d1; T[1][1]=b2/d2; T[2][1]=b3/d3;
		T[0][2]=c1/d1; T[1][2]=c2/d2; T[2][2]=c3/d3;
		break;
	}

}
//==============================================================================
void printMPC(int sNode, int sdof, double diff,AdditionalConstraint *aC, int numACs, ostream* outStream) 
{
	int i;

	(*outStream) << sNode << ' ' << sdof;
	for(i=0;i<numACs;i++)
	   (*outStream) <<  ' ' << aC[i].node << ' ' << aC[i].dof << ' ' << aC[i].factor; 
	(*outStream) <<' '<< diff << '\n';
}
//==============================================================================
void printMPC_2_File(ofstream *os, int sNode, int sdof, double diff,AdditionalConstraint *aC, int numACs) 
{
	int i;

	(*os) << sNode << ' ' << sdof;
	if(numACs<1){
		(*os) <<"No master dofs in the additional constraint!" << endl;
		cout <<"No master dofs in the additional constraint!" << endl;
	}
	for(i=0;i<numACs;i++)
	   (*os) <<  ' ' << aC[i].node << ' ' << aC[i].dof << ' ' << aC[i].factor; 
	(*os) <<' '<< diff << endl;
}
//==============================================================================
void printMPCNeutralFile(ofstream *os, int sNode, int sdof, double diff,AdditionalConstraint *aC, int numACs, NodeGroup *node)
{
	(*os) << (*node)[sNode].getFirstDof() + sdof  ;

	for (int i=0; i<numACs; i++) {
	   (*os) <<  ' ' << (*node)[aC[i].node].getFirstDof() +  aC[i].dof << ' ' << aC[i].factor; 
	}
	(*os) <<' '<< diff << endl;
}
//*********************************************************************************
void getMPC(char *S, int sNode, int sdof, double diff,AdditionalConstraint *aC, int numACs) 
{
	int i,j;
	j  = sprintf(S, "%d  %d ",sNode,sdof);
	for(i=0;i<numACs;i++)
		j += sprintf(S+j, "%d  %d %Lf ", aC[i].node , aC[i].dof, aC[i].factor); 
	j += sprintf(S+j, "%Lf \n", diff);
	// BETA_OUT <<S; // 990527
}
//==============================================================================
//+xtang 990116's version
void getAnEqn(int sNode,int sDof, double Diff, 
			  AdditionalConstraint *ac, int numACs, NodeGroup *node, Equations *equations)
{
	//char S[300];     
	//getMPC(S,sNode, sDof, Diff, ac, numACs); // for output 
	AdditionalEquation *eqn;
	eqn = new AdditionalEquation;
	eqn->setSlave( (*node)[sNode].getFirstDof() + sDof);

	for (int i=0; i<numACs; i++) {
		eqn->addAdditionalTerm( (*node)[ac[i].node].getFirstDof() + ac[i].dof,ac[i].factor);
	}

	eqn->setConstant(Diff);
//	cout <<"-";
	equations->addAdditionalEquation(eqn);
//	cout <<"+";
}
//-xtang 990116's version

//*********************************************************************************
//  readGeneralMPC
//
//  Xiaodong Tang 980908 
//  
//*********************************************************************************
//  Format of data:
//
//      sDof = mDof_1*factor_1 + mDof_2*factor_2 +...+ difference;
//
//      readGeneralMPC
//      SlaveEntityID  SlaveEntityData_1 ...  MasterEntityID  MasterEntityData_1 ...
//      SlaveDof  [MasterDof  Factor]  Difference
//      ...
//      -1
//      ...        
//      exitReadGeneralMPC 
//
//  where:
//      (1) SlaveEntityID,  3 = plane (rectangular region)
//                          2 = line  (line segment)
//                          1 = point
//
//      (2) MasterEntityID,  3 = plane (rectangular region)
//                           2 = line (line segment)
//                           1 = point
//                           0 = dummy node
//      (3) EntityData,  there are 3 node numbers for plane, 
//          2 for line, 1 for point and 0 for dummy node.  User
//          need to get these node numbers from a mesh.  Use Plot2000 for
//          this purpose.
//
//          The convention for specifying an entity:
//          (a) plane:  first 3 vertex node numbers, p1, p2, p3,  ccw or cw
//                    p3----------p4
//                     |           | 
//                     |           | 
//                     |           |
//                    p2----------p1
//          (b) line:
//                    p1-----------p2
//          (c) point:
//                    o p1
//          (d) dummy node:
//                    no data.
//
//      (4) Dof,  0=Ux, 1=Uy, 2=Uz, ...  at user's discretion
//*********************************************************************************
/*     
void Model::readGeneralMPC(istream * inStream) 
{
//code taken off from file 
}
*/

//+xtang 981214
//
//     Usage of this new format for specifying Multipoint Constraints:
//
//     Author:  Xiaodong Tang
//       Date:  1998-Dec.-15
//     Origin:  Texas A&M University
//
//     Entity types:
//     o  Plane			3
//     o  Line			2
//     o  Point			1
//	   o  Dummy Node	0
//
//     Constraint types;
//     -  PlaneToPlane
//     -  PlaneToPoint
//     -  LineToLine
//     -  LineToLine
//     -  PointToPoint


string createNumericString(string prefix, int num, string suffix)
{
	char buffer [10];
    sprintf(buffer,"%d",num);
	string snum(buffer);
	string numericString = prefix + snum + suffix;

	return numericString;

}

void ElasticityModel::ReadMultiPointConstraints(istream * inStream) 
{
	int					numberOfTokens;
	char				*tokenList[20];

	int numNodes=mesh->numNodes;
	int numDims=mesh->numDims;
	NodeGroup *node=&(mesh->node);

	int sNode, i;
	double difference;

	BETA_OUT	<< "Current MPCTransformationTolerance = " << DOUBLE_FORMAT << MPCTransformationTolerance << endl;
	cout		<< "Current MPCTransformationTolerance = " << DOUBLE_FORMAT << MPCTransformationTolerance << endl;


	ofstream* oMPC=filemanager->OpenOutputStream("MPC.txt");
	
	// these data defined by xtang
	int sEnt, mEnt;
	double cs[3][3], cm[3][3];
	double Ts[3][3], Tm[3][3];
	int j;
	int sdof, mdof, ddof, count, numPE;
	double mfact, dfact;
	double sKsi, sEta, mKsi, mEta;

	//JV051304
	for (i=0;i<3;i++){
		for (j=0;j<3;j++){
			cs[i][j]=0;//JV051304 - initializing coordinates to 0 
			cm[i][j]=0;//or it does not pick node properly for 2d case
		}
	}

	int *pS, *pM;
	pS = new int [numNodes];  // open full length of array
	pM = new int [numNodes];  // open full length of array

   	AdditionalConstraint ac[6];  // sufficient number

    while (true) {

		_NextEntity:   // loop for entity

	    BETA_OUT << endl;

		// get specifications of slave and master entities
	
		exitFlag=getLineAndTokenize(inStream,"exitReadMultiPointConstraints",
			                         tokenList, numberOfTokens);
		
		if(exitFlag==1) goto _EndReadMultiPointConstraints_New;
	
		//Get Slave and Master entity modes

		BETA_OUT << endl;
		// Output slave information
		BETA_OUT<<"========================== "<<endl; 
		if(COMPARE(tokenList[0],"PlaneToPlane")==0 ) {
			sEnt=3; mEnt=3;
			BETA_OUT << " PLANE slaved to PLANE\n";
		}
		else if(COMPARE(tokenList[0],"PlaneToPoint")==0 ) {
			sEnt=3; mEnt=1;
			BETA_OUT << " PLANE slaved to POINT\n";
		}
		else if(COMPARE(tokenList[0],"LineToLine")==0 )   {
			sEnt=2; mEnt=2;
			BETA_OUT << " LINE slaved to LINE\n";
		}
		else if(COMPARE(tokenList[0],"LineToPoint")==0 )  {
			sEnt=2; mEnt=1;
			BETA_OUT << " LINE slaved to POINT (LineToPoint)\n";
		}
		else if(COMPARE(tokenList[0],"PointToPoint")==0 ) {
			sEnt=1; mEnt=1;
			BETA_OUT << " POINT slaved to POINT\n";
		}
		else 
		{
			sEnt=-1; mEnt=-1;
			BETA_OUT << " SLAVE and/or MASTER are not defined ! \n";
		}  // Entity not defined
		BETA_OUT<<"========================== "<<endl; 


		// get Slave entity definition
		exitFlag=getLineAndTokenize(inStream,"exitReadMultiPointConstraints",
			tokenList, numberOfTokens);

		if(exitFlag==1) goto _EndReadMultiPointConstraints_New;


		if (COMPARE(tokenList[0],"Node")==0 ) {
            j=0;
			if ( (numberOfTokens-1) < sEnt) {
				BETA_OUT << " Wrong number of nodes in defining an entity...slave entity not defined! " << endl;
			}
			else 
			for (i=1; i <= sEnt; i++) {
				cs[j][0]=(*node)[atoi(tokenList[i])].x;
		        cs[j][1]=(*node)[atoi(tokenList[i])].y;
				cs[j][2]=(*node)[atoi(tokenList[i])].z;
                				
				j++;
			}
		}
		else if (COMPARE(tokenList[0],"Coord")==0) {
			if ( (numberOfTokens-1) < sEnt*numDims) {
				EXIT_BETA("Wrong number of coords in defining an entity...slave entity not defined!");
			}
			if ( (numberOfTokens-1) > sEnt*numDims) {
				BETA_OUT << "This model has a " <<numDims<<"D mesh. Specified entity can be defined using " <<sEnt*numDims<<" coordinates." << endl;
				BETA_OUT << "You have used " <<numberOfTokens-1<<" coordinates to define this entitiy." << endl;
				BETA_OUT << "Excess number of coords used to define an entity. Check your input script to make sure this was intended! " << endl;
			}
			//the following assumes that the mesh is ATLEAST 2D. that might not be a good idea.
			for (j=0; j < sEnt; j++) {
				cs[j][0]=atof(tokenList[j*numDims+1]);
		        cs[j][1]=atof(tokenList[j*numDims+2]);
				if (numDims==3) {cs[j][2]=atof(tokenList[j*numDims+3]);}
			}
		} 

		if (COMPARE(tokenList[0],"exit")==0) goto _NextEntity;

		// Get Master entity definition
		exitFlag=getLineAndTokenize(inStream,"exitReadMultiPointConstraints",
			tokenList, numberOfTokens);
		if(exitFlag==1) goto _EndReadMultiPointConstraints_New;

		if (COMPARE(tokenList[0],"Node")==0 ) {
            j=0;
			if ( (numberOfTokens-1) < mEnt) {
				BETA_OUT << " Wrong number of nodes in defining an entity...master entity not defined! " << endl;
			}
			else 
			for (i=1; i <= mEnt; i++) {
				cm[j][0]=(*node)[atoi(tokenList[i])].x;
		        cm[j][1]=(*node)[atoi(tokenList[i])].y;
				cm[j][2]=(*node)[atoi(tokenList[i])].z;
				j++;
			}
		}
		else if (COMPARE(tokenList[0],"Coord")==0) {
			if ( (numberOfTokens-1) < mEnt*numDims) {
				BETA_OUT << " Wrong number of coords in defining an entity...master entity not defined! " << endl;
			}
			else 
			for (j=0; j < mEnt; j++) {
				cm[j][0]=atof(tokenList[j*numDims+1]);
		        cm[j][1]=atof(tokenList[j*numDims+2]);
				if (numDims==3) {cm[j][2]=atof(tokenList[j*numDims+3]);}
			}
		} 
		else if (COMPARE(tokenList[0],"Dummy")==0) {
			mEnt=DUMMY_NODE; 
		} // For Master only 

		if (COMPARE(tokenList[0],"exit")==0) goto _NextEntity;


		switch (sEnt) {
		case 1:  // master must be a point no matter of what master entity type specified
			sNode=mesh->getNodeIndex(cs[0][0],cs[0][1],cs[0][2]);
			pS[0] = sNode;
			if (mEnt==DUMMY_NODE) {ac[0].node=numNodes-1;} // 981202 xt number of the last node in the mesh
			else ac[0].node=mesh->getNodeIndex(cm[0][0],cm[0][1],cm[0][2]);
			break;
		case 2:
		case 3:  // master must be dummy, a point or a line;
			getTransT(Ts, cs, sKsi, sEta, sEnt);
            if (abs(sEnt)==2) { 
				PointsInEntity(mesh, Ts, cs[0][0], cs[0][1],cs[0][2], sEnt, sKsi, sEta, numPE, pS);
			} 
			else if (abs(sEnt)==3) {
				PointsInEntity(mesh, Ts, cs[1][0], cs[1][1],cs[1][2], sEnt, sKsi, sEta, numPE, pS);
			}
			
			break;
		} // switch


//// *************  
// xtang 19991016
			BETA_OUT<<"SLAVE: "; 
			if (sEnt==3) {BETA_OUT<<"Plane"; }
			else if (sEnt==2) {BETA_OUT<<"Line"; }
			else if (sEnt==1) {BETA_OUT<<"Point"; }

            BETA_OUT << endl;

            // output vertex info
            
			int sDump;

			BETA_OUT << "Points used to define slaved entity: "<<sEnt<<endl;
			for (i=0; i < sEnt; i++) {
			  sDump=mesh->getNodeIndex(cs[i][0],cs[i][1],cs[i][2]);
			  BETA_OUT <<"Point "<<i+1<<":   "<< DOUBLE_FORMAT << cs[i][0]<<"   "<<DOUBLE_FORMAT << cs[i][1]<<"   "<<DOUBLE_FORMAT << cs[i][2]<<"  "<<sDump<<endl;
			}
            
			int irow;

			if (sEnt > 1) {
			BETA_OUT <<numPE << " node(s) selected."<<endl;

			irow=0;
			for (i=0; i< numPE; i++) {
				BETA_OUT << pS[i] << ' ';
				irow++;
				if (irow==20) {
					BETA_OUT << endl;
					irow=0;
				}
			}
			BETA_OUT<<endl;
			}

			BETA_OUT<<endl<<"MASTER: "; 
			if (mEnt==3) {BETA_OUT<<"Plane"; }
			else if (mEnt==2) {BETA_OUT<<"Line"; }
			else if (mEnt==1) {BETA_OUT<<"Point"; }
			else if (mEnt==0) {BETA_OUT<<"Dummy"; }
            BETA_OUT << endl;			
            
			BETA_OUT << "Points used to define master entity: "<<mEnt<<endl;
			for (i=0; i < mEnt; i++) {
			  sDump=mesh->getNodeIndex(cm[i][0],cm[i][1],cm[i][2]);
			  BETA_OUT <<"Point "<<i+1<<":   " <<DOUBLE_FORMAT << cm[i][0]<<"   "<<DOUBLE_FORMAT << cm[i][1]<<"   "<<DOUBLE_FORMAT << cm[i][2]<<"  "<<sDump<<endl;
			}

			if (abs(mEnt) > 1) {
				getTransT(Tm, cm, mKsi, mEta, mEnt);
				if (fabs(mKsi-sKsi) > MPCTransformationTolerance) {
			        cerr << "Regions not identical";
					FatalError("Regions not identical");
				}
				if (abs(mEnt)==2) {
					GetMatchedList(mesh, Ts, cs[0][0], cs[0][1],cs[0][2],	Tm, cm[0][0], cm[0][1],cm[0][2], mEnt, numPE, pS, pM);
				}
				else if (abs(mEnt)==3) {
					GetMatchedList(mesh, Ts, cs[1][0], cs[1][1],cs[1][2],	Tm, cm[1][0], cm[1][1],cm[1][2], mEnt, numPE, pS, pM);
				}

			
				if(mEnt>1) {
				BETA_OUT <<  numPE << " node(s) selected."<<endl;

				irow=0;
				for (i=0; i< numPE; i++) {
					BETA_OUT << pM[i] << ' ';
					irow++;
					if (irow==20) {
						BETA_OUT << endl;
						irow=0;
					}
				}
				BETA_OUT<<endl;
				}

			}
			else if (abs(mEnt)==0) {
				ac[0].node=numNodes-1;
			    BETA_OUT << endl;
				BETA_OUT << ac[0].node << ' ';
			    BETA_OUT<<endl;
			} 
			else if (abs(mEnt)==1) {
				ac[0].node=mesh->getNodeIndex(cm[0][0],cm[0][1],cm[0][2]);
				BETA_OUT << endl;
			    BETA_OUT << ac[0].node << ' ';
			    BETA_OUT<<endl;
			}

//// *************

		BETA_OUT << endl;
		

		_MPCLoopForAnEntity:   // loop for MPC's of an entity

		// get MPC specifications
		exitFlag=getLineAndTokenize(inStream,"exitReadMultiPointConstraints",
			                         tokenList, numberOfTokens);
		if(exitFlag==1) {
			cerr << "Premature end of data block...check !";
			FatalError("Premature end of data block");
		}
		
		count=0;
		
		if (COMPARE(tokenList[0],"exit")==0) goto _NextEntity;

		sdof= abs(atoi(tokenList[0]));

  
		mdof= atoi(tokenList[1]);
		mfact= atof(tokenList[2]);
		
		ac[count].dof=mdof;
		ac[count].factor=mfact;   // ac[count].node to be assigned later !


		// 990527
		string s1 = createNumericString("S",sdof+1,"");
		string s2;
		if (mEnt == 0) {s2=createNumericString("<D",mdof+1,">");} else {s2=createNumericString("M",mdof+1,"");}
		
		BETA_OUT <<s1<<" = "<<s2<<" * "<<mfact<<" + ";  // BCO 032111
		
		i=3;
		count++;  // here we consider the dummy node masters
		while(i+2 < numberOfTokens) {
			ddof = atoi(tokenList[i++]);
			dfact= atof(tokenList[i++]);
			if(count > 5) FatalError("Trying to write more AdditionalConstraints that memory allocated for. Fix!");
			ac[count].node=numNodes-1; 
			ac[count].dof = ddof;
			ac[count].factor = dfact;
			s2=createNumericString("<D",ddof+1,">");

			BETA_OUT << s2<<" * "<<dfact<<" + ";   // BCO 032111

			count++;
		}
		difference = atof(tokenList[i++]);
		
		BETA_OUT << difference<<endl;  //990527

//		cout << "imposing eqn ";
	
		switch (abs(sEnt)) {
		case 1:  // master must be a point no matter of what master entity type specified
			if ((sNode==ac[0].node) && (sdof==ac[0].dof)) break; // skip coincide dof's
			getAnEqn(sNode,sdof, difference, ac, count, &(mesh->node), &equations);
//			printMPC_2_File(oMPC,pS[i],sdof, difference, ac, count);
			printMPCNeutralFile(oMPC, pS[0], sdof, difference, ac, count, &(mesh->node));
			break;
		case 2:
		case 3:  // master must be dummy, a point or a line;
			for (i=0; i<numPE; i++) {
				//+xtang 990114 fixed to skip coincide node and dof's
				//JV030509 why are we skipping coincident nodes? it is handled by the eqn class anyway
				//JV030509 does not work correctly when trying to constrain dofs to 0 using MPC. so not skipping anymore.
				if (abs(mEnt) > 1) {ac[0].node=pM[i]; }
				//if ( (pS[i]==ac[0].node) && (sdof==ac[0].dof) ) continue; // skip coincide dof's //JV030509 not skipping anymore
//				cout << i << " ";
//				if(i%10==0) cout << endl;
				getAnEqn(pS[i],sdof, difference, ac, count, &(mesh->node), &equations);
//				cout << ".";
//				printMPC_2_File(oMPC,pS[i],sdof, difference, ac, count);
				printMPCNeutralFile(oMPC, pS[i], sdof, difference, ac, count, &(mesh->node));
//				cout << "x";
			} // for
/*
				char aelistfilename[200];
				sprintf(aelistfilename,"aelist_%d_%d.txt",sDump,i);
				ofstream* os=filemanager->OpenOutputStream(aelistfilename);
				equations.printAdditionalEquations(os);
				filemanager->CloseOutputStream(os);
*/
		} // switch
//		cout << "... finished imposing eqn" << endl;
		
		goto _MPCLoopForAnEntity;  // highly efficient

	} // while
	_EndReadMultiPointConstraints_New:

	// xtang 000422 output finalized additional equations 
//	equations.printAdditionalEquations();

	BETA_OUT << endl;
//    cout << "exitReadMultiPointConstraints - about to release memory" << endl;
	delete [] pS;
	delete [] pM;
//    cout << "exitReadMultiPointConstraints - released memory" << endl;
	*oMPC<<"exitReadMultiPointConstraints \n";
	filemanager->CloseOutputStream(oMPC);
//    cout << "exitReadMultiPointConstraints - finished" << endl;
	
	return;
}
