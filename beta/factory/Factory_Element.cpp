#include "stdafx.h"

#include "utility/utility.h"
#include "Factory.hpp"

/*========================================================
To remove a element from this project, do the following
	* delete the include file
	* delete the corresponding "IF_CREATE( ....." statement below
========================================================*/
int getNumberOfElementsFromMesh(char* filename);
//====================================================================
#include "elements/BasicElement.hpp"
#include "elements/1D/RodElas.hpp"
#include "elements/3D/ElasticityElement3D.hpp"
#include "elements/2D/ElasticityElement2D.hpp"
//====================================================================
void Factory::ElementFactory(istream *inStream,
							ElementGroup *element, int &numElements, ostream* outStream)
{
	int i,first,last,increment;
	bool *created;
	char name[80];

	char *localTokenList[20];
	int  numberOfTokens;

	getLineAndTokenize(inStream,"end",localTokenList, numberOfTokens);
	if(numberOfTokens==2){
		if(COMPARE(localTokenList[0],"getNumberOfElementsFromMesh")==0)	{
			if(doesFileExist(localTokenList[1])==false) {
				cout << localTokenList[1] << " : this file does not exist" << endl;
				exit(1);
			}
			numElements=getNumberOfElementsFromMesh(localTokenList[1]);
		}else{
			cout << "incorrect command - use getNumberOfElementsFromMesh <meshfilename>" << endl;
			exit(1);
		}
	}else{
		numElements=atoi(localTokenList[0]);
	}
	(*outStream) << "Mesh has <" << numElements << "> elements.\n";
	if(numElements<=0) {	banner("numElements must be > zero.\n", *outStream); }

	created = new bool [numElements];
	for(i=0;i<numElements;i++) created[i] = false;

	element->reserve(numElements);
	
//=============================
top:
	getLineAndTokenize(inStream,"end",localTokenList, numberOfTokens);
	strcpy(name, localTokenList[0]);
	(*outStream) <<"ElementType = "<< name<<endl;
	if(COMPARE(name,"exitCreateElements")==0)goto _EndSetElementType;

	while (true) {
		if( getLineAndTokenize(inStream,"end",localTokenList, numberOfTokens)==1) {
			goto _EndSetElementType;
		}
		if(numberOfTokens==1 && COMPARE(localTokenList[0],"all")==0 ){
			first=0;
			last=numElements-1;
			increment=1;
		}else if(numberOfTokens==3){
			first=atoi(localTokenList[0]);
			last=atoi(localTokenList[1]);
			increment=atoi(localTokenList[2]);
		}else{
			first=atoi(localTokenList[0]);
			if(first<0){(*outStream) <<"\t Exit loop input\n";goto top; }
		}
		(*outStream) <<"\t Loop Format: First = "
			<<first<<"  Last= "<<last
			<<"  Increment= "<<increment<<endl;

		//Error checking...
		if(first < 0 || first > numElements || first > last
		             || last  > numElements || increment<1)
		   {banner("Bad parameters...", *outStream);  exit(1);}

		//Create elements...
		if (!createElements(first, last, increment, element, name, created))  {
			//No match found !
			(*outStream) <<"You specified type= "<<name<<endl;
			(*outStream) << "This type is not defined " << endl;
			exit(1);
		}

	}//end of while

_EndSetElementType:
	for(i=0;i<numElements;i++) {
		if(!created[i]) {
			(*outStream) << "Element number '" << i << "' not specified.\n";
			banner("Bad element specification.\n", *outStream);
			exit(1);
		}
	}
	delete [] created;

	return;
}  // end of CreateElements

bool Factory::createElements(int first, int last, int increment, 
							ElementGroup *element, char* name, bool* created)
{
    BasicElement *e=0;
    Create_element_ifTokenIs(BasicElement)	else
    Create_element_ifTokenIs(RodElas)		else
    Create_element_ifTokenIs(ElasticityElement3D)	else 
    Create_element_ifTokenIs(ElasticityElement2D)	else
    return false;

}

//======================================================================
