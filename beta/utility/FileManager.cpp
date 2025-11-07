#include "stdafx.h"

#include <string>
using namespace std;

#include "FileManager.hpp"
#include "utility/directoryManipulation.hpp"

void directory( const char *filepath,char *dirpath ); 
//char *TrimName(char *sFileName);
string TrimNameOnly(string sFileName);
bool doesDirExist(const char * fileName);
bool makeDir(const char * dirPath);
string removeExtensionFromFilename(string filename);

//========================================================================
FileManager* createFileManager(int type)
{
	switch(type){
		case COMMON_OUTPUT_SCHEME:
			return new FileManager();
			break;
		case QUALIFIED_OUTPUT_SCHEME:
			return new QualifiedFileManager();
			break;
		default:
			cout <<"Cannot recognize filemanager type"<<endl;
			exit(1);
	}
}
//========================================================================
FileManager::FileManager()
{
	outStream=&cout;
	filenamingoption=NoPrefix;
//	ostreamListCount=0;
}
//========================================================================
FileManager::~FileManager()
{
	if(outStream != &cout)
		CloseOutputStream((ofstream*)outStream);

}
//========================================================================
bool FileManager::restorePath()
{
	_chdir(inputDirectory.c_str()); 
	return true;
}
//========================================================================
bool FileManager::initialize(string filename) 
{
	cout <<"Initializing FileManager..." << endl;

char pathName[_MAX_PATH];
char inputPath[_MAX_PATH];
directory(filename.c_str(),pathName); //gives relative path to the actual file
_chdir(pathName); //change dir to the actual dir with the input file
_getcwd( inputPath, _MAX_PATH ); // save this directory into inputdirectory
inputDirectory=inputPath;

cout<<"Update current directory based on location of input filename"<<endl;
cout<<"CurrentDirectory is now " <<inputDirectory<<endl;

	char trimmedName[_MAX_PATH];
	strcpy(trimmedName, filename.c_str());
	inputFilename=TrimNameOnly(trimmedName);	

	setOutputScheme();
	setOutputDirectory();
	return true;
}
//========================================================================
bool FileManager::setOutputDirectory() 
{
	restorePath();
	if(!doesDirExist(OutputDirectoryName.c_str()))
		if( ! makeDir( OutputDirectoryName.c_str() )) {
			cout << "unable to create output directory - " << OutputDirectoryName.c_str() << endl;
			exit(1);
		}
	OutputPrefix=OutputDirectoryName;

#if defined(WIN32) || defined(WIN64) //needed because commands like 'del' don't work as intended with '/'
	OutputPrefix += '\\';
#else
	OutputPrefix += '/';
#endif

	return true;
}
//========================================================================
bool FileManager::setDefaultOutputStream()
{
	string ofile = OutputPrefix + outFilename;
	outStream = new ofstream(ofile.c_str());
	return true;
}
//========================================================================
bool FileManager::setDefaultOutputStream(string name)
{
	outFilename=name;
	return setDefaultOutputStream();
}
//========================================================================
ofstream* FileManager::OpenOutputStream(string name)
{
	restorePath();

	//fix to include different naming options...
	//save the filenames to a list ?
	ofstream *outs;
	string ofile = OutputPrefix + name;
	outs = new ofstream(ofile.c_str() );
	return outs;
}
//========================================================================
ifstream* FileManager::OpenInputStream(string name)
{
	ifstream *outs=0;
	restorePath();
	outs = new ifstream(name.c_str());
	if(outs->good())
		return outs;
	else{
		cout << "This input file cannot be opened :" << name << endl;
		cout << "The current working directory is " << inputDirectory << endl;
		cout << "\nListing the contents of the current working directory....\n" << endl;
#if defined(WIN32) || defined(WIN64)
		system("dir");
#endif
		exit(1);
		return 0;
	}
}
//========================================================================
string FileManager::GetOutputFilenameWithRelativePath(string name)
{
	//fix to include different naming options...
	//save the filenames to a list ?
	string ofile = OutputPrefix + name;
	return ofile;
}
//========================================================================
fstream* FileManager::OpenBinaryOutputStream(string name)
{
	restorePath();

	//fix to include different naming options...
	//save the filenames to a list ?
	fstream *outs;
	string ofile = OutputPrefix + name;
	outs = new fstream(ofile.c_str(),ios::trunc |ios::in | ios::out | ios:: binary ); 
	return outs;
}
//========================================================================
bool FileManager::CloseOutputStream(ofstream *os)
{
	//fix to include different naming options...
	//save the filenames to a list ?
	if(os!=NULL && os != &cout){	
		os->close(); 
		delete os;
		os=0;
	}
	return true;
}
//========================================================================
bool FileManager::CloseOutputStream(ostream *os)
{
	return CloseOutputStream( (ofstream*) os);
}
//========================================================================
bool FileManager::CloseInputStream(ifstream* is)
{
	if(is!=NULL){	
		is->close(); 
		delete is;
	}
	return true;
}
//========================================================================
/*
bool FileManager::OpenOutputFile(char* name)
{
	//fix to include different naming options...
	ostreamList[ostreamListCount] = new ofstream(name);
	ostreamListCount++;
	return true;
}
*/
//========================================================================
bool FileManager::setOutputScheme()
{
	outFilename="Output.txt";
	OutputDirectoryName="results";
	return true;
}
//========================================================================
bool QualifiedFileManager::setOutputScheme()
{
	string scriptname=removeExtensionFromFilename(inputFilename);
	outFilename="Output.txt";
	OutputDirectoryName=scriptname+"_results";
	return true;
}
//========================================================================
string FileManager::GetFullOutputDirectoryPath()
{
	string fullpath;
#if defined(WIN32) || defined(WIN64)
	fullpath += inputDirectory + '\\' + OutputPrefix;
#else
	fullpath += inputDirectory + '/' + OutputPrefix;
#endif
	return fullpath;
}
//========================================================================
string FileManager::GetOutputPrefix()
{
	return OutputPrefix;
}
//========================================================================
void FileManager::ChangeCurrentDirectory(string newPath)
{
	_chdir(newPath.c_str()); 
}
//========================================================================
void FileManager::ChangeCurrentDirectoryToInputDirectory()
{
	_chdir(inputDirectory.c_str()); 
}
//========================================================================
bool FileManager::AddDirectoryToSearchPath(string newPath)
{
	SearchPath.push_back(newPath);
	return true;
}
//========================================================================
ifstream* FileManager::OpenInputFileFromSearchPath(string filename)
{
	ifstream *is=0;
//	if(SearchPath.size() == 0)
//		(*outStream) << "No directories in the BETA search path. Searching in the current working directory." << endl;

	for(int i=0; i<(int)SearchPath.size(); i++){
		ChangeCurrentDirectory(SearchPath[i]);
		is=new ifstream(filename.c_str());
		if(is->good()){
			restorePath();
			return is;
		}
		else{
			delete is;
		}
	}

	restorePath();
	return 0;
}
//========================================================================
