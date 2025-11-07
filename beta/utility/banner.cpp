#include "stdafx.h"
#include "FileManager.hpp"

using namespace std;

void banner( char * string, ostream &outStream)
{int i,length;
	length = (int)strlen(string);
	length += 2;
	outStream<<"*";
	for(i=0; i<length; i++)
		outStream<<".";
	outStream<<"*"<<endl;

	 outStream<<"* "<< string << " * "<<endl;

	outStream<<"*";
	for(i=0; i<length; i++)
		 outStream<<".";
	outStream<<"*"<<endl;
}
//=====================================================
void banner( char * string, ostream *outStream)
{int i,length;
	length = (int)strlen(string);
	length += 2;
	(*outStream)<<"*";
	for(i=0; i<length; i++)
		(*outStream)<<".";
	(*outStream)<<"*"<<endl;

	 (*outStream)<<"* "<< string << " * "<<endl;

	(*outStream)<<"*";
	for(i=0; i<length; i++)
		 (*outStream)<<".";
	(*outStream)<<"*"<<endl;
}

