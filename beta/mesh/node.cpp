#include "stdafx.h"
#include "utility/utility.h"

#include "node.hpp" 

//==========================================================
Node::Node()
{
	Initialize();
}

//==========================================================
Node::Node(const Point &n)
{
	Initialize();
	x=n.x; y=n.y; z=n.z;
}
//==========================================================
Node::Node(const int NodeNum, const double X, const double Y, const double Z)
{
	Initialize();
	x=X; y=Y; z=Z;
	nodeNum = NodeNum;
	newNodeNum = NodeNum;
}
//==========================================================
void Node::Initialize()
{
	x=y=z=0.0;
	nodeNum=newNodeNum=-1;
	firstDof=-1;
	numDof=0;
	numElementsAtNode=0;
	eList=0;
	isActive=true;
	isSelected=false;
	isMidNode=false;
	isFiber=false;
}
//==========================================================
Node::~Node()
{
	if (numElementsAtNode > 0) delete [] eList;
	numElementsAtNode=0;
}

//==========================================================
void Node::readAllNodalCoordinates( Node * node,
												int numNodes,
												int numDims,
												istream *inputStream)
{
cout<<"This is not up to date!... don't use now\n";  exit(1);
/*	for(int i=0; i< numNodes; i++){
		int number;
		(*inputStream)>>number;
		node[number].read(*inputStream,numDims);
		node[number].nodeNum=number;
	}
	*/
}

//==========================================================
bool Node::read(istream &istrm,CoordinateReadOption readOption)
{
 switch(readOption) {
case X:
		istrm >> x;
		y=z=0;
		break;
case XY:
		istrm >> x >> y;
		z=0;
		break;
case YZ:    // Quasi-3D mesh  
		//pre-existing capability to plot '2d-element' in 3D space (e.g. shell/plate elements) have been REMOVED
		istrm >> y >> z;
		x=0;
		break;
case XYZ:
		istrm >> x >> y >> z;
		break;
	default:
		return false;

/*	case 1:
		istrm >> x;
		y=z=0;
		break;
	case 2:
		istrm >> x >> y;
		z=0;
		break;
	case -2:    // Quasi-3D mesh  
		//pre-existing capability to plot '2d-element' in 3D space (e.g. shell/plate elements) have been REMOVED
		istrm >> y >> z;
		x=0;
		break;
	case 3:
		istrm >> x >> y >> z;
		break;
	default:
		return false;
		*/

//(*defaultOutStream)<<"x y z="<<x<<"  "<<y<<"  "<<z<<endl;

	}//end of switch

	return true;
}
//==========================================================
void Node::readBinary(FILE *fp)
{
	const int SIZEOF_INT	=sizeof(int);
	const int SIZEOF_DOUBLE	=sizeof(Coordinate);

	fread(&nodeNum, SIZEOF_INT, 1, fp);
	fread(&x, SIZEOF_DOUBLE, 1, fp);
	fread(&y, SIZEOF_DOUBLE, 1, fp);
	fread(&z, SIZEOF_DOUBLE, 1, fp);
}
//========================================================== 
void Node::writeBinary(FILE *fp)
{
	const int SIZEOF_INT	=sizeof(int);
	const int SIZEOF_DOUBLE	=sizeof(Coordinate);

	fwrite(&nodeNum, SIZEOF_INT, 1, fp);
	fwrite(&x, SIZEOF_DOUBLE, 1, fp);
	fwrite(&y, SIZEOF_DOUBLE, 1, fp);
	fwrite(&z, SIZEOF_DOUBLE, 1, fp);
}
//==========================================================
ostream& operator<<(ostream &o, const Node& node)
{
	o << node.nodeNum << " : (" << node.x << ',' << node.y << ',' << node.z << 
		")\n\tfirstDof: " << node.firstDof << "\tnumDof: " << node.numDof << endl;
	return o;
}
//==========================================================
void Node::print(ostream &ostrm)
{
	ostrm << (*this);
}
//==========================================================
void Node::write(ostream &ostrm, int numDims)
{
	ostrm << nodeNum << ' ' << x << ' ' << y;
	if(numDims==3) ostrm << ' ' << z << '\n';
}

void Node::writeNodeInfo(ostream &ostrm, int numDims) // xtang 031902
{
	ostrm <<nodeNum<<"\t"<<newNodeNum<<"\t"<<numDof<<"\t"<<firstDof<<"\t"<<isActive<<endl;
}
//==========================================================
void Node::setNumDofPerNode(istream * inStream, ostream *out,
							int numNodes, NodeGroup *node, int &totalNumDof)
{
	int numDofAtNode,first,last,increment,i,dof;
	bool *dofSet,dofCheck;	
	dofSet = new bool[numNodes];
	for(i=0;i<numNodes;i++) dofSet[i] = false;

	char *localTokenList[20];
	int  numberOfTokens;

	while (true) {
		if( getLineAndTokenize(inStream,"end",localTokenList, numberOfTokens)==1) {
			goto ExitSetNumDofPerNode;
		}
		if(numberOfTokens==2 && COMPARE(localTokenList[1],"all")==0 ){
			numDofAtNode=atoi(localTokenList[0]);
			first=0;
			last=numNodes-1;
			increment=1;
		}else if(numberOfTokens==4){
			numDofAtNode=atoi(localTokenList[0]);
			first=atoi(localTokenList[1]);
			last=atoi(localTokenList[2]);
			increment=atoi(localTokenList[3]);
		}else{
			numDofAtNode=atoi(localTokenList[0]);
			if(numDofAtNode<0) goto ExitSetNumDofPerNode;
		}

		(*out)<<first<<" "<<last<<" "<<increment<<endl;
		(*out)<< "Number of DOF = " << numDofAtNode << " for nodes "
				<< first <<" to "<<last << " by " << increment << '\n';

		if ( (first < 0) || (last >= numNodes) ) {
			cerr << "Node::setNumDofPerNode : specified first and last node numbers - " << first << " and " << last << endl;
			cerr << "numNodes =" << numNodes << endl;
			cerr << "you are trying to set the dofs for a non-existant node!" << endl;
			exit(1);
		}

		//_Set:
		for(i=first;i<=last && i<numNodes;i+=increment) {
			dofSet[i] = true;
			(*node)[i].setNumDof(numDofAtNode);
		}
	}
	ExitSetNumDofPerNode:
	totalNumDof = 0;
	for(i=0;i<numNodes;i++) totalNumDof += (*node)[i].getNumDof();
	(*out) << "TotalNumDof = " << totalNumDof << endl;
	
	(*out)<<"\t**Setting first dof for each node"<<endl;
	for(i=0,dof=0;i<numNodes;i++) {
		(*node)[i].setFirstDof(dof);
		dof += (*node)[i].getNumDof();
	}
	
	dofCheck = true;
	for(i=0;i<numNodes;i++)
	if(dofSet[i]==false) dofCheck = false;
	if(!dofCheck) {
		(*out) << "\n\n";
		for(i=0;i<numNodes;i++)
		if(dofSet[i]==false)
				(*out) << "Node " << (*node)[i].nodeNum << ": numDof not set...\n";
		(*out)<<"Some dof not set..."<<endl;
	}
	delete [] dofSet;
}

//==========================================================
