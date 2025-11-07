#include "stdafx.h"
//#include "utility/utility.h"

#include "InterpAndSolDerivatives.hpp"

//====================================================================
extern	int			verboseFlag;
//====================================================================
InterpAndSolDerivatives::InterpAndSolDerivatives(int max_interp, int max_nodes)
{
	allocated=false;
	MAX_INTERP=0;
	MAX_NODES=0;

	p_u1_x1 = p_u1_x2 = p_u1_x3 = 0;
	p_u2_x1 = p_u2_x2 = p_u2_x3 = 0;
	p_u3_x1 = p_u3_x2 = p_u3_x3 = 0;
	p_u4_x1 = p_u4_x2 = p_u4_x3 = 0;

	setArrayLimits(max_interp, max_nodes);
	numberOfInterp = max_interp;//JV072009: initially assume numberOfInterp is max_interp. if it is different, user has to change later when using the object.
}
//====================================================================
InterpAndSolDerivatives::~InterpAndSolDerivatives(void)
{ 
	if(allocated==false)return;

	delete [] S;
	delete [] p_S_Lx1;
	delete [] p_S_Lx2;
	delete [] p_S_Lx3;
	delete [] p_S_x1;
	delete [] p_S_x2;
	delete [] p_S_x3;
	delete [] xCoor;
	delete [] yCoor;
	delete [] zCoor;

//up to 6 fields
	delete [] nodalValues_u1;
	delete [] nodalValues_u2;
	delete [] nodalValues_u3;
	delete [] nodalValues_u4;
	delete [] nodalValues_u5;
	delete [] nodalValues_u6;
}
//====================================================================
void InterpAndSolDerivatives::setArrayLimits(int max_interp, int max_nodes)
{ 
	MAX_INTERP	=max_interp;
	MAX_NODES	=max_nodes;
}
//====================================================================
bool InterpAndSolDerivatives::allocate()
{ 
	if(MAX_INTERP==0 || MAX_NODES==0)
		return false;

	S=new double[MAX_INTERP];
	p_S_Lx1=new double[MAX_INTERP];
	p_S_Lx2=new double[MAX_INTERP];
	p_S_Lx3=new double[MAX_INTERP];
	p_S_x1=new double[MAX_INTERP];
	p_S_x2=new double[MAX_INTERP];
	p_S_x3=new double[MAX_INTERP];
	xCoor=new double[MAX_NODES];	
	yCoor=new double[MAX_NODES];
	zCoor=new double[MAX_NODES];

//up to 6 fields
	nodalValues_u1=new double[MAX_NODES];
	nodalValues_u2=new double[MAX_NODES];
	nodalValues_u3=new double[MAX_NODES];
	nodalValues_u4=new double[MAX_NODES];
	nodalValues_u5=new double[MAX_NODES];
	nodalValues_u6=new double[MAX_NODES];

	allocated=true;

	return true;
}
//====================================================================
