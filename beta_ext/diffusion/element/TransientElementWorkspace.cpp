#include "stdafx.h"

#include "TransientElementWorkspace.hpp"

//--------------------------------------------------------------------
TransientElementWorkspace::TransientElementWorkspace(void)
{
	F1=0;
	F2=0;
	M_deltaq=0;
}
//--------------------------------------------------------------------
TransientElementWorkspace::~TransientElementWorkspace(void)
{
	if(F1) delete [] F1; F1=0;
	if(F2) delete [] F2; F2=0;
	if(M_deltaq) delete [] M_deltaq; M_deltaq=0;
}
//--------------------------------------------------------------------
void TransientElementWorkspace::PrintLimits(ostream *os)
{
	//print derived class limits first 
	//... in this case, no new limits are set in this class

	//then print parent class limits
	ElementWorkspace::PrintLimits(os);
}
//--------------------------------------------------------------------
