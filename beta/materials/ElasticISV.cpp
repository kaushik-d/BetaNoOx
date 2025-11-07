#include "stdafx.h"
#include "ElasticISV.hpp"

ElasticISV::ElasticISV():
BasicISV(),
Temperature(0.0),
Moisture(0.0)
{
}

ElasticISV::ElasticISV(const ElasticISV &Arg):
BasicISV(Arg),
Temperature(Arg.Temperature),
Moisture(Arg.Moisture)
{
}