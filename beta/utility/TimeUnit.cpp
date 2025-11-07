#include "stdafx.h"

#include "utility.h"
#include <string>


string getTimeUnitString(BETA_Time_Unit timeunit)
{
	switch(timeunit){
	case seconds:
		return "seconds";
		break;
	case minutes:
		return "minutes";
		break;
	case hours:
		return "hours";
		break;
	default:
		return "Undefined Time Unit";
		break;
	}

}
