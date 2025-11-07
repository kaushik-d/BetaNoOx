#include "stdafx.h"

#include "excepts.hpp"
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

int ErrorHandler::counter=0;
ErrorHandler::ErrorHandler()
{
	myID = counter++;
	className[0]='\0';
	methodName[0]='\0';
}
ErrorHandler::ErrorHandler(char *cN)
{
	myID = counter++;
	SetClassName(cN);
	methodName[0]='\0';
}
ErrorHandler::ErrorHandler(char *cN, char *mN)
{
	myID = counter++;
	SetClassName(cN);
	SetMethodName(mN);
}
ErrorHandler::~ErrorHandler()
{
}
void ErrorHandler::_FatalError(char *filename, int linenumber, char *message) const
{
	cerr << "\nFatal Error:  ";
	if(strlen(className)) {
		cerr << className;
		if(strlen(methodName))
			cerr << "::" << methodName;
		cerr << " - ";
	}
	cerr << message << endl;
	cerr << filename << " : " << linenumber << endl << endl << "Exiting\n\n";
	exit(1);
}	
void ErrorHandler::_FatalError(char *filename, int linenumber, string message) const
{
	cerr << "\nFatal Error:  ";
	if(strlen(className)) {
		cerr << className;
		if(strlen(methodName))
			cerr << "::" << methodName;
		cerr << " - ";
	}
	cerr << message << endl;
	cerr << filename << " : " << linenumber << endl << endl << "Exiting\n\n";
	exit(1);
}	
void ErrorHandler::_Warning(char *filename, int linenumber,char *message) const
{
	cerr << "\nWarning:  ";
	if(strlen(className)) {
		cerr << className;
		if(strlen(methodName))
			cerr << "::" << methodName;
		cerr << " - ";
	}
	cerr << message << endl;
	cerr << filename << " : " << linenumber << endl;
}
void ErrorHandler::SetClassName(char *str) 
{
	strncpy(className,str,79);
}
void ErrorHandler::SetMethodName(char *str) 
{
	strncpy(methodName,str,79);
}
void ErrorHandler::ShowClassMethod() const
{
	cout << "==> " << className << "::" << methodName << " <==" << endl;
}
