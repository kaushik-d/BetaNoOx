// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#if defined(WIN32) || defined(WIN64)
 #pragma message("Compiling stdafx.h - this should happen just once per project.\n")
 #define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
 #include <tchar.h>
#endif



// TODO: reference additional headers your program requires here
	#include <string>
	#include <algorithm>
	#include <cctype>
	#include <iostream>
	#include <iomanip>
	#include <fstream>
	#include <cstdlib>
	#include <cstdio>
	#include <vector>
