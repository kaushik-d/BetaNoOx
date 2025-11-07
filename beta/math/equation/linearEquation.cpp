#include "stdafx.h"

#include <stdlib.h>
#include "math/matrix.hpp"
#include "math/equation/linearEquation.hpp"
#include "utility/excepts.hpp"

LinearEquations::LinearEquations(const unsignedInt size) : K(size,size)
{
	numberOfEquations = size;
	// Allocate
	if((F = new Double [size])==0) {
		cerr << "Memory allocation error in Linear Equations.\n";
		exit(1);
	}
	// Initialize
	zeroLoadVector();
}
LinearEquations::~LinearEquations()
{
	delete [] F;
}
void LinearEquations::zeroYourself()
{
	zeroMatrix();
	zeroLoadVector();
}

void LinearEquations::zeroMatrix()
{	K = 0; }

void LinearEquations::zeroLoadVector()
{
	for(unsignedInt i=0;i<numberOfEquations;i++) F[i] = 0;
}


