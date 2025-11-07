#include "stdafx.h"

#include "utility/utility.h"
#include "utility/FileManager.hpp"
#include "BoundaryNodeList.hpp"

void BoundaryNodeList::printBoundaryNodeList(const char *aName)
{
	ofstream os(aName);
	printBoundaryNodeList(&os);
	os.close();
}

void BoundaryNodeList::printBoundaryNodeList(ostream *os)
{
	int i, j=0;
	Node *n;
	*os << numNodesInList << "   Nodes \n";
	for (i=0; i<numNodesInList; i++) {
		n=nList[i];
		*os <<n->nodeNum <<"\t";
		j++;
		if (j>=10) {
			*os <<endl;
			j=0;
		}
	}
	*os << endl;
	*os << "mesh information for a substructure super element with the given boundary nodes" << endl;
	*os << numNodesInList << " 1 3" << endl;
	for (i=0; i<numNodesInList; i++) {
		n=nList[i];
		*os <<i <<"\t";
		*os << n->x << " " << n->y << " " << n->z << endl;
	}
	*os << "0 " << numNodesInList <<"\t";;
	for (i=0; i<numNodesInList; i++) {
		*os << i <<"\t";
		j++;
		if (j>=10) {
			*os <<endl;
			j=0;
		}
	}
}

void BoundaryNodeList::printBoundaryNodeListOnly(ostream *os)
{
	int i, j=0;
	Node *n;
	//*os << numNodesInList << "   Nodes \n";
	for (i=0; i<numNodesInList; i++) {
		n=nList[i];
		*os <<n->nodeNum << endl;
		j++;
	}
	*os << endl;
}

