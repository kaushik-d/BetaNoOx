#include "stdafx.h"
#include "NodeSet.hpp"
#include "utility/utility.h"

/////////////////////////////////////////////////////////////////////////////
//Helper function prototypes
bool nodeOnPlane(const Node* nodeToCheck, 
                 const double normal[], 
                 const double pointOnPlane[], 
                 const double tolerance,
                 double& distance);

void crossProduct3D(const double vector1[], 
                    const double vector2[], 
                    double result[]);

double dotProduct3D(const double vector1[], 
                    const double vector2[]);

void normalizeVector3D(double vector[]);

void scaleVector3D(double vector[], 
                   const double& scaleFactor);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Default Constructor
NodeSet::NodeSet():
    sortedNodes(0),
    globalSortedSet(),
    nodeSearchTolerance(1.0e-8)
{   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

NodeSet::NodeSet(const double &tolerance):
    sortedNodes(0),
    globalSortedSet(),
    nodeSearchTolerance(tolerance)
{   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Copy Constructor
NodeSet::NodeSet(const NodeSet &arg):
    sortedNodes(arg.sortedNodes),
    globalSortedSet(arg.globalSortedSet),
    nodeSearchTolerance(arg.nodeSearchTolerance)
{   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NodeSet::addNode(Node* nodeToAdd)
{
    if(globalSortedSet.count(nodeToAdd)==0)
    {
        globalSortedSet.insert(nodeToAdd);
        sortedNodes.push_back(nodeToAdd);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool NodeSet::read(istream * fin, BasicMesh* mesh)
{
    string command;
    *fin >> command; //First read in how to process the rest of the line

    //================   Node   ================//
    if(COMPARE(command.c_str(),"Node") == 0){
        int nodeNum;
        *fin >> nodeNum;
        if(nodeNum < 0 || nodeNum >= mesh->numNodes){
            cout << "Invalid Node number specified for NodeSet.\n";
            exit(1);
        }
        addNode(&(mesh->node[nodeNum]));
        return true;
    }
    //=============  MultiNode  ================//
    else if( COMPARE(command.c_str(),"MultiNode")==0 ){
        int first, last, increment;
        *fin >> first >> last >> increment;
        if(first < 0 || 
           last >= mesh->numNodes ||
           first > last){
            cout << "Invalid Node number specified for NodeSet.\n";
            exit(1);
        }
        for(int i = first;i<=last;i+=increment){
            addNode(&(mesh->node[i]));
        }
        return true;
    }
    //=============  NodeNumberList  ===========//
    else if(COMPARE(command.c_str(),"NodeNumberList") == 0){
        int numNodes, thisNodeNum;
        *fin >> numNodes;
        for(int i=0;i<numNodes;++i){
            *fin >> thisNodeNum;
            if(thisNodeNum < 0 || thisNodeNum >= mesh->numNodes){
                cout << "Invalid Node number specified for NodeSet.\n";
                exit(1);
            }
            addNode(&(mesh->node[thisNodeNum]));
        }
        return true;
    }
    //==============  Point  ===================//
    else if(COMPARE(command.c_str(),"Point") == 0){
        string restOfLine;
        getline(*fin, restOfLine);
        vector<double> CoordVals = TokenizeStringToVectorOfDoubles(restOfLine);
        double coords[] = {0.0,0.0,0.0};
		switch(CoordVals.size()){
			case 3:
				coords[2]=CoordVals[2];
			case 2:
				coords[1]=CoordVals[1];
			case 1:
				coords[0]=CoordVals[0];
		}
        addNodesOnPoint(mesh,coords);
        return true;
    }
    //=============  NodesOnPlane  =============//
    else if(COMPARE(command.c_str(),"NodesOnPlane") == 0){
        int dir;
        double coord;
        *fin >> dir >> coord;
        addNodesOnPlane(mesh,dir,coord);
        return true;
    }
    //===========  OrderedNodesOnLine  ==========//
    // This command needs to be generalized so that you can input a node numbers or a coordinate to
    // define a given point...
    else if(COMPARE(command.c_str(),"OrderedNodesOnLine") == 0){
        double coord1[] = {0.0,0.0,0.0};
        double coord2[] = {0.0,0.0,0.0};
        string restOfLine;
        getline(*fin, restOfLine);
        vector<double> CoordVals = TokenizeStringToVectorOfDoubles(restOfLine);
        switch(CoordVals.size()){
			case 2:
				coord1[0] = CoordVals[0];
                coord2[0] = CoordVals[1];
                break;
			case 4:
				coord1[0] = CoordVals[0];
                coord1[1] = CoordVals[1];
                coord2[0] = CoordVals[2];
                coord2[1] = CoordVals[3];
                break;
			case 6:
				coord1[0] = CoordVals[0];
                coord1[1] = CoordVals[1];
                coord1[2] = CoordVals[2];
                coord2[0] = CoordVals[3];
                coord2[1] = CoordVals[4];
                coord2[2] = CoordVals[5];
                break;
            default:
                cout << "Invalid number of coordinates given in defining NodeSet with OrderedNodesOnLine.\n";
                exit(1);
		}
        addOrderedNodesOnLineSegment(mesh,coord1,coord2);
        return true;
    }
    //==========  OrderedNodesOnPlane  ==========//
    // This command needs to be generalized so that you can input a node numbers or a coordinate to
    // define a given point...
    else if(COMPARE(command.c_str(),"OrderedNodesOnPlane") == 0){
        double coord1[] = {0.0,0.0,0.0};
        double coord2[] = {0.0,0.0,0.0};
        double coord3[] = {0.0,0.0,0.0};
        string restOfLine;
        getline(*fin, restOfLine);
        vector<double> CoordVals = TokenizeStringToVectorOfDoubles(restOfLine);
        switch(CoordVals.size()){
			case 3:
				coord1[0] = CoordVals[0];
                coord2[0] = CoordVals[1];
                coord3[0] = CoordVals[2];
                break;
			case 6:
				coord1[0] = CoordVals[0];
                coord1[1] = CoordVals[1];
                coord2[0] = CoordVals[2];
                coord2[1] = CoordVals[3];
                coord3[0] = CoordVals[4];
                coord3[1] = CoordVals[5];
                break;
			case 9:
				coord1[0] = CoordVals[0];
                coord1[1] = CoordVals[1];
                coord1[2] = CoordVals[2];
                coord2[0] = CoordVals[3];
                coord2[1] = CoordVals[4];
                coord2[2] = CoordVals[5];
                coord3[0] = CoordVals[6];
                coord3[1] = CoordVals[7];
                coord3[2] = CoordVals[8];
                break;
            default:
                cout << "Invalid number of coordinates given in defining NodeSet with OrderedNodesOnPlane.\n";
                exit(1);
		}
        addOrderedNodesOnFiniteRectangularPlane(mesh,coord1,coord2,coord3);
        return true;
    }
    return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NodeSet::addNodesOnPoint(BasicMesh *mesh, 
                              double coords[])
{
    //Note - implicitly assumes that nodes in 2D and 1D models have extra coordinates (z or y,z)
    //set to zero, which they should based on Node::initialize()

    //Stores the distance from the node to the input coordinates
    double distance=0.0;

    //Loop through the mesh and add any nodes within the tolerance distance to sortedNodes if it isn't already in sortedNodes
    for(int i=0;i<mesh->numNodes;++i){
        distance = sqrt((mesh->node[i].x-coords[0])*(mesh->node[i].x-coords[0]) + 
                        (mesh->node[i].y-coords[1])*(mesh->node[i].y-coords[1]) +
                        (mesh->node[i].z-coords[2])*(mesh->node[i].z-coords[2]));
        if(distance<nodeSearchTolerance){
            addNode(&(mesh->node[i]));
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NodeSet::addNodesOnPlane(BasicMesh* mesh,int dir,double coord)
{
    double NormVec[] = {0.0,0.0,0.0};
    double Coord[] = {0.0,0.0,0.0};
    double distance;

    NormVec[dir-1] = 1.0;
    Coord[dir-1] = coord;
    for(int i=0; i<mesh->numNodes; ++i){
        if( nodeOnPlane(&(mesh->node[i]),NormVec,Coord,nodeSearchTolerance,distance) ){
            addNode(&(mesh->node[i]));
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NodeSet::addOrderedNodesOnLineSegment(BasicMesh* mesh,
                                           double startCoords[],
                                           double endCoords[])
{
    //Note - implicitly assumes that nodes in 2D and 1D models have extra coordinates (z or y,z)
    //set to zero, which they should based on Node::initialize()

    //Store the geometry data into geometryCoords
    geometryCoords.push_back(vector<Point>(0));
    geometryCoords.back().push_back(Point(startCoords));
    geometryCoords.back().push_back(Point(endCoords));

    //Indices of nodes sorted by local coordinates
    set<NodeLocalCoordData,NodeLocalCoordComparator> nodesToAdd;

    //More for readablitiy, define the coordinates of the two endpoints of the line
    //segment to be (x0,y0,z0) and (x1,y1,z1)
    double x0 = startCoords[0], y0 = startCoords[1], z0 = startCoords[2];
    double x1 = endCoords[0], y1 = endCoords[1], z1 = endCoords[2];

    //Vx, Vy, Vz are the components of a vector going from one end of the segment to another
    double Vx = x1 - x0;
    double Vy = y1 - y0;
    double Vz = z1 - z0;

    //Define the length of the line segment.  Storing the square will save some calculations
    //in the loop over the nodes in the mesh
    double segmentLengthSquared = Vx*Vx + Vy*Vy + Vz*Vz;
    double segmentLength = sqrt(segmentLengthSquared);


    //The x, y, and z coordinates of the node being checked
    double xNode,yNode,zNode;

    //The parametric variable t for the closest point on the line to the node
    double t;

    //The components of a vector going from the first
    //point on the line segement to the node
    double Vx0Node,Vy0Node,Vz0Node;

    //The components of a vector going from the second
    //point on the line segement to the node
    double Vx1Node,Vy1Node,Vz1Node;

    //The minimum distance from the node to the line passing through the endpoints
    double distance;

    //Now that these values are defined, loop through the mesh to find the nodes which are on the line segment
    //Each node is checked in ta manner that is outlined by the Wolfram math website at 
    //http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    
    for(int i=0;i<mesh->numNodes;++i){
        //Get the coordinates of the node being examined:
        xNode=mesh->node[i].x;
        yNode=mesh->node[i].y;
        zNode=mesh->node[i].z;

        //Define a vector going from the first point on the segment to the node
        Vx0Node = xNode - x0;
        Vy0Node = yNode - y0;
        Vz0Node = zNode - z0;

        //Define a vector going from the second point on the segment to the node
        //Vx1Node = xNode - x1;
        //Vy1Node = yNode - y1;
        //Vz1Node = zNode - z1;
        
        //Determine the distance along the segment where the intersection occurs
        t = (Vx*Vx0Node + Vy*Vy0Node + Vz*Vz0Node)/segmentLengthSquared;

        //If the node is beyond either endpoint of the line segment, then determine
        //if it is within nodeSearchTolerance of the endpoint.  If so, it is on
        //the line segment.  Otherwise, it is not.

        if( t < 0.0){
            if(sqrt(Vx0Node*Vx0Node + Vy0Node*Vy0Node + Vz0Node*Vz0Node) < nodeSearchTolerance) 
            {
                NodeLocalCoordData thisNode;
                thisNode.nodeNum = i;
                thisNode.LocalCoords.push_back(t*segmentLength);
                nodesToAdd.insert(thisNode);
                continue;
            }else{continue;}
        }
        if( t > 1.0){
            //calculate vector going from the second point on the segment to the node
            Vx1Node = xNode - x1;
            Vy1Node = yNode - y1;
            Vz1Node = zNode - z1;

            if(sqrt(Vx1Node*Vx1Node + Vy1Node*Vy1Node + Vz1Node*Vz1Node) < nodeSearchTolerance) 
            {
                NodeLocalCoordData thisNode;
                thisNode.nodeNum = i;
                thisNode.LocalCoords.push_back(t*segmentLength);
                nodesToAdd.insert(thisNode);
                continue;
            }else{continue;}
        }

        //Now calculate the minimum distance from the node to the line 
        //that passes between the endpoints.

        //Cross product between vector from segment start to end and 
        //vector from segment start to node

        double CrossX = Vy*Vz0Node - Vz*Vy0Node;
        double CrossY = Vz*Vx0Node - Vx*Vz0Node;
        double CrossZ = Vx*Vy0Node - Vy*Vx0Node;

        distance = sqrt(CrossX*CrossX + CrossY*CrossY + CrossZ*CrossZ)/segmentLength;

        //Determine if the distance is within the tolerance and add the node to 
        //sortedNodes if it is (and isn't already in sortedNodes)
        if(distance<nodeSearchTolerance){
            NodeLocalCoordData thisNode;
            thisNode.nodeNum = i;
            thisNode.LocalCoords.push_back(t*segmentLength);
            nodesToAdd.insert(thisNode);
        }
    }
    //zip through the sorted set of nodes on the line and add them to the nodeSet
    set<NodeLocalCoordData,NodeLocalCoordComparator>::iterator it;
    for(it=nodesToAdd.begin(); it != nodesToAdd.end(); ++it){
        addNode(&(mesh->node[it->nodeNum]));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NodeSet::addOrderedNodesOnFiniteRectangularPlane(BasicMesh* mesh,
                                                      double point0[],
                                                      double point1[],
                                                      double point2[])
{
    //Note - implicitly assumes that nodes in 2D and 1D models have extra coordinates (z or y,z)
    //set to zero, which they should based on Node::initialize()

    //Store the geometry data into geometryCoords
    geometryCoords.push_back(vector<Point>(0));
    geometryCoords.back().push_back(Point(point0));
    geometryCoords.back().push_back(Point(point1));
    geometryCoords.back().push_back(Point(point2));

    //First it must be clearly understood how the finite plane is defined from the input coordinates.
    //The finite plane is a rectangle with one edge running from point0 to point1 (point0 and point1
    //are vertices of the rectangle).  Point2 is colinear to the opposite edge of the rectangle (but
    //is not necessarily a vertex).

    //Note that on the rectangle, direction 1 indicates the direction from point0 to point1.  
    //Direction 2 indicates the direction perpendicular to direction 1 going from point0 towards
    //(but not necessarily directly to) point2

    //Indices of nodes sorted by local coordinates
    set<NodeLocalCoordData,NodeLocalCoordComparator> nodesToAdd;

    //Determine the componenets of the plane's normal vector
    double vecEdge1[3], vecEdge2[3], vec0to2[3], planeNormal[3];

    for(int i=0;i<3;++i){
        vecEdge1[i] = point1[i]-point0[i];
        vec0to2[i] = point2[i]-point0[i];
    }

    crossProduct3D(vecEdge1,vec0to2,planeNormal);
    normalizeVector3D(planeNormal);

    //Now create a vector running from the point0 along the edge perpendicular to 
    //the 0 - 1 edge (this will more clearly define the rectangle)
    crossProduct3D(planeNormal,vecEdge1,vecEdge2);
    normalizeVector3D(vecEdge2);

    //Determine the rectangle's dimensions
    double rectangleDim1 = sqrt(dotProduct3D(vecEdge1,vecEdge1));
    double rectangleDim2 = dotProduct3D(vecEdge2,vec0to2);

    //From this point on, vecEdge1 and vecEdge2 will be normalized (their magnitudes
    //are the dimensions of the rectangle).  This will simplify the process to determine
    //if a point is on the rectangle or not
    normalizeVector3D(vecEdge1);

    //This will be a vector from point0 to the node
    double vec0toNode[3];

    //These doubles store the location of a node relative to point0 in the rectangle's local
    //coordinate system.  These values will be used to determine the sorting of the nodes
    double nodeDir1Position, nodeDir2Position;

    double distance;
    
    for(int i=0;i<mesh->numNodes;++i){
        //First, determine if the node is within the specified distance of the plane
        if(nodeOnPlane(&(mesh->node[i]),planeNormal,point0,nodeSearchTolerance,distance)){
            //Node is on plane, determine if it is within the finite rectangle.
            vec0toNode[0] = mesh->node[i].x - point0[0];
            vec0toNode[1] = mesh->node[i].y - point0[1];
            vec0toNode[2] = mesh->node[i].z - point0[2];

            nodeDir1Position = dotProduct3D(vecEdge1,vec0toNode);
            nodeDir2Position = dotProduct3D(vecEdge2,vec0toNode);

            //A little logic to determine if the node lies within all four edges plus or minus the tolerance.
            if(    ( nodeDir1Position > -nodeSearchTolerance && nodeDir1Position < rectangleDim1 + nodeSearchTolerance ) 
                && ( nodeDir2Position > -nodeSearchTolerance && nodeDir2Position < rectangleDim2 + nodeSearchTolerance ) )
            {
                //Add the node if it isn't in sortedNodes already
                NodeLocalCoordData thisNode;
                thisNode.nodeNum = i;
                thisNode.LocalCoords.push_back(nodeDir1Position);
                thisNode.LocalCoords.push_back(nodeDir2Position);
                nodesToAdd.insert(thisNode); // Adds the node in order by dir1, then dir2
            }
        }
    }//end loop over nodes in mesh

    //zip through the sorted set of nodes on the plane and add them to the nodeSet
    set<NodeLocalCoordData,NodeLocalCoordComparator>::iterator it;
    for(it=nodesToAdd.begin(); it != nodesToAdd.end(); ++it){
        addNode(&(mesh->node[it->nodeNum]));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int NodeSet::findNode(const Node* nodeToFind) const
{
    //Returns the first index where nodeToFind is found, or -1 if it is not found
    vector<Node*>::const_iterator result = find(sortedNodes.begin(),sortedNodes.end(),nodeToFind);
    if(result == sortedNodes.end()) {return -1;}
    else return (int) (result-sortedNodes.begin());
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NodeSet::print(ostream &OS) const
{
    //prints a summary of the node set
    OS << "      # of Nodes: " << sortedNodes.size() << "     Nodes:" << endl;
    for(int i=0;i<sortedNodes.size();i++){
        if(i && i%10==0) OS << endl;
        OS << sortedNodes.at(i)->nodeNum << '\t';
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//int NodeSet::nodeIndexFromCoordinates(const double& coords) const
//{
//    //Returns the first index of the node in the nodeset with the given coordinates 
//    //(within nodeSearchTolerance).  Returns -1 if no node exists at those coordinates
//
//    vector<Node*>::const_iterator node;
//    return 0;
//}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Non-class Helper Functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool nodeOnPlane(const Node* nodeToCheck, 
                 const double normal[], 
                 const double pointOnPlane[], 
                 const double tolerance, 
                 double& distance)
{
    //Returns true if the node is within tolerance distance of the plane defined by normal and pointOnPlane
    //See http://mathworld.wolfram.com/Point-PlaneDistance.html 

    distance = normal[0]*(nodeToCheck->x - pointOnPlane[0]) +
               normal[1]*(nodeToCheck->y - pointOnPlane[1]) +
               normal[2]*(nodeToCheck->z - pointOnPlane[2]);

    return(fabs(distance)<=tolerance);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void crossProduct3D(const double vector1[], 
                    const double vector2[], 
                    double result[])
{
    //performs vector1 x vector2 and stores the solution in result
    result[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
    result[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
    result[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double dotProduct3D(const double vector1[], 
                    const double vector2[])
{
    //performs vector1 . vector2 and returns the scalar result
    return vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void normalizeVector3D(double vector[])
{
    //modifies the input vector so its magnitude is unity
    double magnitude = sqrt(dotProduct3D(vector,vector));
    vector[0] /= magnitude;
    vector[1] /= magnitude;
    vector[2] /= magnitude;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void scaleVector3D(double vector[], 
                   const double& scaleFactor)
{
    //scales the magnitude of vector by scaleFactor
    vector[0] *= scaleFactor;
    vector[1] *= scaleFactor;
    vector[2] *= scaleFactor;
}