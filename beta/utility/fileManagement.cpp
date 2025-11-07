#include "stdafx.h"

#include "utility_defines.h"
#include "directoryManipulation.hpp"

#if defined(WIN32) || defined(WIN64)
	#include <direct.h> 
#else
	#include <sys/stat.h> 
	#include <dirent.h> 
#endif

#include <string>
using namespace std;

int getLineAndTokenize(istream * inStream, char *terminateString,
							  char   *tokenList[], int &numberOfTokens);

//==========================================================
bool doesFileExist(string fileName)
{
	ifstream *checkFile;
	checkFile = new ifstream(fileName.c_str());
	bool status = checkFile->good();
	if(status==true){cout<<"The specified file '"<<fileName<<"' exists\n";}
	else {
		cout<<"The specified file '"<<fileName<<"' does not exist\n";
		char currentDirectory[_MAX_PATH];
		_getcwd( currentDirectory, _MAX_PATH );
		cout<<"CurrentDirectory ="<<currentDirectory<<endl;

	}
	checkFile->close();
	delete checkFile;
	return status;
}
//==========================================================
bool doesFileExist_win(char * fileName)
{
 //Check whether file exists
 char pathbuffer[_MAX_PATH];

#if defined(WIN32) || defined(WIN64)
 _searchenv(fileName, ".", pathbuffer);
#endif

 if( *pathbuffer != '\0' )
	{return true;}
   else
   {return false; }
}
//==========================================================
bool doesDirExist(const char * fileName)
{
 //Check whether file exists (work for directories too !!!)
 char pathbuffer[_MAX_PATH];

#if defined(WIN32) || defined(WIN64)
	_searchenv(fileName, ".", pathbuffer);
	if( *pathbuffer != '\0' )
		{return true;}
	else
		{return false; }
#else
	DIR* dir = opendir(fileName);
	if (!dir)return false;
	else {
		closedir(dir);
		return true;
	}
#endif
}
//==========================================================
bool makeDir(const char * dirPath)
{
	//removeDir(dirPath);  //uses system command
//	deleteDir(dirPath);  //uses windows API function - 
#if defined(WIN32) || defined(WIN64)
	if( _mkdir(dirPath) == 0 ) return true; else return false;
#else
	if( mkdir(dirPath, 0700) == 0 ) return true; else return false;
#endif

}
//==========================================================
void createFile(istream * inStream)
{
	char filename[80];
 	char line[256];
	(*inStream)>> filename;
	ofstream *outStream = new ofstream(filename);
top:
 	(*inStream).getline(line, 80,'\n');
	if(COMPARE(line,"exitCreateFile")==0){outStream->close(); return;}
	(*outStream)<<line<<endl;
 goto top;
}
//=======================================================
ostream * openFile(char * filename)
{ 
	ostream * pointer;

	pointer = new ofstream(filename);
	if(!(*pointer))	{
		cout <<"Cannot open file you selected:"<<filename<<endl;      
		exit(0);
	}
	return(pointer);
}

//=======================================================
ifstream * openFileForInput(char * filename)
{ 
//	ifstream * pointer= new ifstream(filename, ios::nocreate); // old c++ std
	ifstream * pointer= new ifstream(filename, ios_base::in);
//	ifstream * pointer= new ifstream(filename, ios::binary|ios::nocreate); // xtang 030101

	if( (pointer->good()) ){ return(pointer); }
	else{
		cout <<"Cannot open file you selected:  " <<filename<<endl;       
		exit(0);
		return(0);
	}
}
//====================================================================
void combineFiles(istream * inStream)
{
// This function will combine files like those listyed in the output from win3d8h	
char filename[80];
char line[256];
char *tokenList[10];
 int i,exitFlag,numberOfTokens;
 (*inStream)>> filename;
 ofstream *outStream = new ofstream(filename);
top:

exitFlag=getLineAndTokenize(inStream,"end",
			                       tokenList, numberOfTokens);

if(COMPARE(tokenList[0],"end")==0){(*outStream)<<"end"; outStream->close(); return;} 

if(COMPARE(tokenList[0],"openFile")==0){
	ifstream *newInStream = new ifstream(tokenList[1]);
	(*outStream)<<tokenList[2]<<endl;

	while( (*newInStream).getline(line, 255,'\n')  ){
	(*outStream)<<line<<endl;	 
	}
    newInStream->close();
	goto top;
} 


for(i=0; i<numberOfTokens; i++)
{(*outStream)<<tokenList[i]<<"  ";} (*outStream)<<endl;

 goto top;
}
// ===================================================================
void deleteFiles(string files)
{
    string cdump;

#if defined(WIN32) || defined(WIN64)
	cdump = "del " + files;
#else
	cdump = "rm -f " + files; // unix
#endif
	system(cdump.c_str());

}
// ===================================================================

