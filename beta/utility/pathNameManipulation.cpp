#include "stdafx.h"

#include <iomanip>
#include <stdlib.h>
#include <string>
using namespace std;

#if defined(WIN32) || defined(WIN64)
  #include <direct.h>  // xtang: 980331, not usable for UNix
#else
   #include <unistd.h>  
#endif


#include <string.h>

//==================================================================
void formOutputFileName(char *outFile, char * baseName, char *tag)
{
// Combines the three strings:  currentDirectory  baseName    tag

 char currentDirectory[_MAX_PATH];

#if defined(WIN32) || defined(WIN64)
 _getcwd( currentDirectory, _MAX_PATH );
#else
  getcwd( currentDirectory, _MAX_PATH );
#endif

 //Check whether file name has entire path ... search for colon

 if( strchr(baseName,':')== NULL){
 strcpy(outFile,currentDirectory);
#if defined(WIN32) || defined(WIN64)
 strcat(outFile,"\\");
#else
 strcat(outFile,"/");
#endif
 }
 strcat(outFile,baseName);
 strcat(outFile,tag);
}

//=================================================================
void directory( const char *filepath,  char *dirpath )
{

	// INPUT : filepath -  relative (or full) path to file
	// OUTPUT: dirpath  -  relative Path to directory 
	unsigned short i;

#if defined(WIN32) || defined(WIN64)
	if(strstr(filepath,"\\")==NULL)
	{strcpy(dirpath,".\\");
    return;
    };

	for( i=strlen(filepath)-1; i>=0; i-- )
		if( filepath[i] == '\\' || filepath[i] == '/' )
		{
			strncpy( dirpath, filepath, i+1 );
			dirpath[i+1] = '\0';
			break;
		}
#else
	if(strstr(filepath,"/")==NULL)
	{strcpy(dirpath,"./");
    return;
    };

	for( i=strlen(filepath)-1; i>=0; i-- )
		if( filepath[i] == '/')
		{
			strncpy( dirpath, filepath, i+1 );
			dirpath[i+1] = '\0';
			break;
		}

#endif
}
//==================================================================
/*
char *TrimName(char *sFileName)
{
// I think this extracts the filename from the full pathname
    int len,i;
    int k=0;
    char *sDest;   //whit
    sDest = sFileName; //whit

    len=i=strlen(sFileName);

#if defined(WIN32) || defined(WIN64)
    while ((k<2) && (i>0)) {  
       i--;
       if (sFileName[i] =='\\') k++;
       else if (sFileName[i] ==':') {
          if (sFileName[i+1] =='\\') i++; // modified here
          break;
       }        
       if (k==1) break;
    }
    if (sFileName[i] =='\\') i++;  // 100197-tang

    for (k=i; k<=len; k++) sDest[k-i]=sFileName[k]; 
    return sDest;
#else
    while ((k<2) && (i>0)) {  
       i--;
       if (sFileName[i] =='/') k++;
       else if (sFileName[i] ==':') {
          if (sFileName[i+1] =='/') i++; // modified here
          break;
       }        
       if (k==2) break;
    }
    if (sFileName[i] =='/') i++;  // 100197-tang

    for (k=i; k<=len; k++) sDest[k-i]=sFileName[k]; 
    return sDest;
#endif
}
//=========================================================================
char *TrimNameOnly(char *sFileName)
{
//JV111702 - this extracts the filename ONLY from the full path name
    int len,i;
    int k=0;
    char *sDest;   //whit
    sDest = sFileName; //whit

    len=i=strlen(sFileName);

#if defined(WIN32) || defined(WIN64)
    while ((k<1) && (i>0)) {  
       i--;
       if (sFileName[i] =='\\') k++;
       else if (sFileName[i] ==':') {
          if (sFileName[i+1] =='\\') i++; // modified here
          break;
       }        
       if (k==1) break;
    }
    if (sFileName[i] =='\\') i++;  // 100197-tang

    for (k=i; k<=len; k++) sDest[k-i]=sFileName[k]; 
    return sDest;
#else
    while ((k<1) && (i>0)) {  
       i--;
       if (sFileName[i] =='/') k++;
       else if (sFileName[i] ==':') {
          if (sFileName[i+1] =='/') i++; // modified here
          break;
       }        
       if (k==1) break;
    }
    if (sFileName[i] =='/') i++;  // 100197-tang

    for (k=i; k<=len; k++) sDest[k-i]=sFileName[k]; 
    return sDest;
#endif
}
*/
//=========================================================================
string TrimNameOnly(string sFileName)
{
//JV111702 - this extracts the filename ONLY from the full path name
	size_t found;
#if defined(WIN32) || defined(WIN64)
	found=sFileName.rfind("\\");
#else
	found=sFileName.rfind("/");
#endif

	if(found==string::npos)
		return sFileName;

	found++;
	string result=sFileName.substr(found);
	return result;
}
//=========================================================================
string removeExtensionFromFilename(string filename)
{
	int loc=filename.rfind('.',filename.length());
	string sub;
	sub.insert(0,filename, 0, loc);
	return sub;

}
//=========================================================================
