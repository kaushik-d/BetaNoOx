#include "stdafx.h"

#include "models/BasicModel.hpp"
#include "BC.hpp"

void Constraint::ReadDOFConstraint(istream *inStream, BasicModel *model)
{
	int i;
	while(true){
		*inStream >> i;
		if (i<0)return;
		doflist.push_back(i);
	}

}