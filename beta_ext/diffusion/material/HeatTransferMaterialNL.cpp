#include "stdafx.h"

#include "HeatTransferMaterialNL.hpp"
#include "math/matrix.hpp"
#include "math/tensor.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

//==========================================================================
void HeatTransferMaterialNL::UpdateGlobalCmat(double normc)
{
//adding non-linearity to the diffusivity
	(*globalCmat)[0][0]=(*globalCmat)[0][0] * normc*normc;
}
//==========================================================================
