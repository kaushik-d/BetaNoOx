#ifndef utility_defines_h
#define utility_defines_h

using namespace std;


// Macro function definitions============================
#if defined(WIN32) || defined(WIN64)
#define COMPARE _stricmp
#include <string.h>
#else
#define COMPARE strcasecmp
#include <strings.h>
#endif

#define INT_FORMAT setw(12)
#define DOUBLE_FORMAT setw(24)<<setprecision(16)
#define DECIMAL_FORMAT setw(12)<<setprecision(10)

#define SET_SCIENTIFIC(a) a << resetiosflags(ios::floatfield) << setiosflags(ios::scientific)
#define SET_FIXED(a) a << resetiosflags(ios::floatfield) << setiosflags(ios::fixed)


#define ChangeToUpper(command) transform (command.begin(),command.end(), command.begin(), (int(*)(int))toupper);

#define IS_COMMAND(COMMAND) \
    if(COMPARE(token, #COMMAND)==0){                     \
		   banner(#COMMAND, BETA_OUT) ;     \
           COMMAND(currentStream); OK                     \
	}

//=========================================================

#define IS_COMMAND_2(COMMAND) \
    if(COMPARE(token, #COMMAND)==0){                     \
		   banner(#COMMAND, BETA_OUT ) ;     \
           COMMAND(currentStream, BETA_OUT); OK                     \
	}

#define OK  foundMatch = 1;


#endif
