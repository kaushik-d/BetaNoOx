#include "stdafx.h"

#include "contourData.hpp"
#include "mesh/BasicMesh.hpp"
#include "utility/utility.h"
/////////////////////////////////////////////////////

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cctype>
#include <cmath>

//====================================================================
extern	int			verboseFlag;
//====================================================================
//void resetMaterialNumberer(PlotInterface *);
bool readDisplacements(istream &, int, int, Node **);

int    getLineAndTokenize(istream * inStream, 
								  char *terminateString,
							     char   *tokenList[], 
								  int &numberOfTokens);
//===========================	
void averageDataValues(int numNodes, int numElem, ElementGroup &element, bool averageOnlyActive, int MaxNumberOfMaterialsForContouring);

#define ACTIVE true
#define INACTIVE false

//===================================================================
ContourData::ContourData(void)
{
	elementNodalValuesState=false;
	elementNodalValues=NULL;
	fixed_option=false;group_option=false;
	active_option=false;
	averageOnlyActive = false;
	doNormalize=false;
	doAveraging=false;
	pickColumn = 1;
	numColumns=-1;
	MaxNumberOfMaterialsForContouring=0;
	numElemAlloc    = 0;
	numColumnsAlloc = 0;
	for(int i=0; i<10;i++){normalizingFactor[i]=1.;}
	filename="notDefined";
	filemanager=0;
};
//===================================================================
ContourData::~ContourData(void)
{
//BETA_OUT<<"ContourData::~ContourData"<<endl;
	if(filemanager)
		freeMemory();
};
//===================================================================
void ContourData::setPointers(BasicMesh* mesh_ptr)
{ 
	mesh = mesh_ptr;
//	pi=mesh_ptr->pi;
}
//===========================================================================
void ContourData::readNormalizingFactors(istream &fp)
{
	int i;
	int  numberOfTokens;
	char *tokenList[30];

	getLineAndTokenize(&fp, "exit",tokenList,numberOfTokens);
	for(i=0; i<numberOfTokens; i++){
		normalizingFactor[i] = atof(tokenList[i]);
	}
}
//==========================================================
bool ContourData::readElementNodalValues(istream &fp)
{
	ElementGroup &element=mesh->element;
int numNodes = mesh->numNodes;
int numElem = mesh->numElements;

int	minLocation=-1,maxLocation=1;
	
	double    fixedMin, fixedMax;
	double	minGroup=0.f,maxGroup=0.f,minActive=1e20,maxActive=-1e20;
	
	int	select_group=0;
	char	fileIsBinary = false;
	
int i,j;

/*
for(i=0; i<numElem; i++){
element[i].setMaterialGroup(1); //Fix this !!!!!!
 }
*/
//==========

char * tokenList[10], name[80];
int numberOfTokens,tokNum;
active_option = true;
int exitFlag;
int dispData=0;
ifstream p_dispFile,p_stressFile;
int first,last,inc,group;
enum ContourDataType{stress,displacements,
                     SetElementMaterialAngles,SetElementMaterial};
ContourDataType dataType;

/////////
streampos m_posStartData;

//////////
topOfLoop:
m_posStartData = fp.tellg();   // save current location of pointer into file
exitFlag = getLineAndTokenize(&fp,"end",tokenList,numberOfTokens);
if(exitFlag==1) goto exitCommandLoop;
tokNum=0;

	if (COMPARE(tokenList[tokNum], "displacements") == 0) {
		fp.seekg( m_posStartData ); 
		// This rewinds the file to put "displacements" back in stream
		dispData=1;
	    dataType = displacements;
        tokNum++;
//eventually change this to use method BasicMesh::readDisplacements
	    readDisplacements(fp,mesh->numNodes,mesh->numDims,
			&mesh->displacements);
		//Hardwire for allocating for 6 columns now
		numColumns = 6;
		pickColumn = 1;
        convertNodalSetToElementSet();
		doAveraging=false;
	    goto exitCommandLoop;
	}

   else if (COMPARE(tokenList[tokNum], "SetElementMaterial") == 0) {
	   dataType = SetElementMaterial;
      // Eventually merge with "readMaterialGroups(fp,pi)"
      tokNum++;
      for(i=0;i<100;i++){
			fp>>first>>last>>inc;
			if(first<0)break;
			fp>>group;
			for(j=first;j<last+1;j+=inc)element[j].setMaterialGroup(group);
         }//end of for
	  goto topOfLoop;
   }

	else if (COMPARE(tokenList[tokNum], "SetElementMaterialAngles") == 0) {
	   dataType = SetElementMaterialAngles;
      tokNum++;
      int elementNumber, rotAxis;
	   for(i=0; i<numElem; i++)
			{fp>>elementNumber>>rotAxis;
			for(j = 0; j<element[elementNumber].numNodesPerElement;j++)
				{
             fp>>element[i].nodalValues[j];
             }
	      }
	   goto topOfLoop;
		}
	else if (COMPARE(tokenList[tokNum], "stress") == 0 || 
		     COMPARE(tokenList[tokNum], "strain") == 0) {
	   dataType = stress;
      tokNum++;
		if(!readValues(fp)) 
			{BETA_OUT << "Error: <readElementNodalValues> - Cannot read values." 
               << endl;
				return false;}
      BETA_OUT<<"Finished reading file: command = stress "<<endl;
      goto exitCommandLoop; //Note jump from loop => make last command!
	}

////////////
	else if (COMPARE(tokenList[tokNum], "stressFile") == 0) {
	   dataType = stress;
      tokNum++;
//		fp >> numColumns;
//		fp >> pickColumn;
		fp>>name;
	    p_stressFile.open(name); 
		if(!readValues(p_stressFile)) 
			{cerr << "Error: <readElementNodalValues> - Cannot read values." 
               << endl;
				return false;}
		p_stressFile.close();
      goto exitCommandLoop; //Note jump from loop => make last command!
	}
//////////////
	else if (COMPARE(tokenList[tokNum], "stressFile_noCol") == 0) {
        // pickColumn must have been set earlier to use this option!
		dataType = stress;
        tokNum++;
		fp>>name;
	    p_stressFile.open(name); 
		if(!readValues(p_stressFile)) 
			{BETA_OUT << "Error: <readElementNodalValues> - Cannot read values." 
               << endl;
				return false;}
		p_stressFile.close();
      BETA_OUT<<"Finished reading file: stressFile_noCol "<<endl;
      goto exitCommandLoop; //Note jump from loop => make last command!
	}
///////////////////
	else if (COMPARE(tokenList[tokNum], "fixed") == 0) {
		fixed_option = true;
		fp >> fixedMin >> fixedMax;
      tokNum++;goto topOfLoop;
	}
	else if (COMPARE(tokenList[tokNum], "group") == 0) {
		group_option = true;
		fp >> select_group;
		select_group--;
      tokNum++;goto topOfLoop;
	}
	else if(COMPARE(tokenList[tokNum],"active") == 0)
		{active_option = true; tokNum++;}
	else if(COMPARE(tokenList[tokNum], "actavg") == 0) {
		active_option = true;
		averageOnlyActive = true;
      tokNum++;goto topOfLoop;
	} else {
		return false;
	}

goto topOfLoop;
exitCommandLoop:

//=================================================

if(doAveraging){
	if(MaxNumberOfMaterialsForContouring <1 ) 
		EXIT_BETA("MaxNumberOfMaterialsForContouring is less that 1. Not an acceptable value!");
if(verboseFlag > Basic) 	displayTime("readElementNodalValues: start averaging", BETA_OUT);
	averageDataValues(numNodes,numElem,element,averageOnlyActive, MaxNumberOfMaterialsForContouring);	
if(verboseFlag > Basic) 	displayTime("readElementNodalValues: start normalizing", BETA_OUT);
}
if(doNormalize) {
	normalizeDataValues(fixed_option, group_option, active_option);
if(verboseFlag > Basic) 	displayTime("readElementNodalValues: finished normalizing", BETA_OUT);
}


BETA_OUT<<"averageOnlyActive="<<averageOnlyActive<<endl;
BETA_OUT<<"fixed_option="<<fixed_option<<endl;
BETA_OUT<<"group_option="<<group_option<<endl;
BETA_OUT<<"active_option="<<active_option<<endl;
BETA_OUT<<"min val="<<MinValue<<endl;
BETA_OUT<<"max val="<<MaxValue<<endl;

	return true;
}


//==========================================================================
bool ContourData::switchColumn(int newCol)
{
//JV040109 : column numbering starts from 1
// Need to add flag to check whether data has already been normalized?
// Should save ?
if(numColumns<0) return false; //JV041803
if(newCol > numColumns || newCol < 1) return false; //JV041803


	ElementGroup &element=mesh->element;
int numNodes = mesh->numNodes;
int numElem = mesh->numElements;

int	minLocation=-1,maxLocation=1;
	
	double	minGroup=0.f,maxGroup=0.f,minActive=1e20,maxActive=-1e20;
	
	int	select_group=0;
//	char	fileIsBinary = false;
	pickColumn = newCol;

copyData();

if(verboseFlag > Basic) displayTime("readElementNodalValues: start averaging", BETA_OUT);
	if(doAveraging) averageDataValues(numNodes,numElem,element,averageOnlyActive, MaxNumberOfMaterialsForContouring);
if(verboseFlag > Basic) displayTime("readElementNodalValues: start normalizing", BETA_OUT);
	if(doNormalize) normalizeDataValues(fixed_option, group_option, active_option);

//	getMinMaxValues(mesh->minVal,mesh->maxVal); // decide what to do about this !!!!
//	mesh->statusFlags.haveData = true;
//	mesh->minContourAuto = mesh->minVal;
//	mesh->maxContourAuto = mesh->maxVal;
BETA_OUT<<"averageOnlyActive="<<averageOnlyActive<<endl;
BETA_OUT<<"fixed_option="<<fixed_option<<endl;
BETA_OUT<<"group_option="<<group_option<<endl;
BETA_OUT<<"active_option="<<active_option<<endl;
BETA_OUT<<"min val="<<MinValue<<endl;
BETA_OUT<<"max val="<<MaxValue<<endl;
return true; //JV041803
}
//==========================================================================
void ContourData::copyData()
{
BETA_OUT<<"copyData"<<endl;
int i,n;
	ElementGroup &element=mesh->element;
int numNodes = mesh->numNodes;
int numElem = mesh->numElements;
double factor;
factor = 1./normalizingFactor[pickColumn-1];
	for(i=0; i< numElem; i++) {
		for(n = 0; n<element[i].numNodesPerElement;n++) {
		  element[i].nodalValues[n] = elementNodalValues[pickColumn-1][i][n]*factor;
		}// n		
	} //  i
BETA_OUT<<"copyData...exit"<<endl;
}
//===========================================================================
bool ContourData::readValues(istream &fp)
{
if(verboseFlag > Basic) displayTime("ContourData::readValues", BETA_OUT);
	ElementGroup &element=mesh->element;
int      numElem = mesh->numElements;

	int i,n,elementNumber,materialGroup;
	int elemOffset = mesh->elementOffset;
	/* Read in ASCII Data */

	if(!fp) cerr << "Error: <readValues> - fp is bad on entry..Crash!" << endl;

// Determine number of columns and then rewind to actually read file
int  numberOfTokens;
char *tokenList[30];


//streampos m_posStartData;
//m_posStartData = fp.tellg(); //not going to use tellg and seekg becuase it does not work properly when using UNIX format files

if(verboseFlag > Basic) displayTime("readValues: start", BETA_OUT);

MaxNumberOfMaterialsForContouring=0;
fp>>elementNumber>>materialGroup;
if(materialGroup>MaxNumberOfMaterialsForContouring)
MaxNumberOfMaterialsForContouring=materialGroup;
	
getLineAndTokenize(&fp, "exit",tokenList,numberOfTokens);
numColumns=0;
numColumns = numberOfTokens;

if(verboseFlag > Basic) displayTime("readValues: start memory allocation", BETA_OUT);

allocateMemory();



//==================================================
//fp.seekg( m_posStartData );  // This "rewinds" the file
//not going to use tellg and seekg becuase it does not work properly when using UNIX format files
//instead going to the beginning of the file and skip the first line.
fp.seekg (0, ios::beg);  // This goes back to the beginning of the file
getLineAndTokenize(&fp, "exit",tokenList,numberOfTokens);


//Not compatible with disp. contours !!

/* Let's not mess with binary for now
if(doesFileExist(binaryFilename))
{readValues_binary(numColumns); return true;}
*/
	for(i=0; i< numElem; i++) {
		fp >> elementNumber >> materialGroup;
		if(materialGroup>MaxNumberOfMaterialsForContouring)
			MaxNumberOfMaterialsForContouring=materialGroup;

		if(!fp) {
			cerr << "i = " << i << endl;
			cerr << elementNumber << " " <<materialGroup << endl;
			cerr <<"Error: Material data. Cannot read elemNum or material group.\n";
			return false;
		}

		elementNumber -= elemOffset;
		// Convert material group to its color group
		element[elementNumber].setMaterialGroup(materialGroup);
		//convertMaterialGroupToColorGroup(materialGroup);


//==============
		for(n = 0; n<element[elementNumber].numNodesPerElement;n++) {
         for(int  column =0; column<numColumns; column++)
				{
					if(!(fp >> elementNodalValues[column][i][n])) {
					cerr << "Error: Material data. Cannot read data. d2";
					return false;
					}// if
					}//column loop
			}// n loop ( loop over nodes )
	}// i loop (= element #)
	MaxNumberOfMaterialsForContouring++;//need to add one because c arrays start from 0
//====================================================


mesh->allocateNodalValuesElements();
copyData();

if(verboseFlag > Basic) displayTime("readValues: finished copyData()", BETA_OUT);

//Not compatible with disp. contours !!

// Eventually change writeValues_binary to add enough format info 
//  that ASCII file is not needed  (eg. number of columns & type of file)


/* Let's not mess with binary for now
 
if(! doesFileExist(binaryFilename) )
{writeValues_binary(numColumns); }

*/


if(verboseFlag > Basic) displayTime("ContourData::readValues ...exit", BETA_OUT);
	return true;
}
//==========================================================================

void ContourData::convertNodalSetToElementSet() //istream &fp)
{
BETA_OUT<<"convertNodalSetToElementSet"<<endl;
	ElementGroup &element=mesh->element;
int numElem = mesh->numElements;

int i;

Node *displacements;
displacements = mesh->displacements;

// This depends on the displacements being read in earlier!


allocateMemory();
int nodeNum;
//==============
for(i=0; i<numElem; i++){ 
	for(int n = 0; n<element[i].numNodesPerElement;n++) {
		    nodeNum = element[i].node[n].nodeNum;
			elementNodalValues[0][i][n] = displacements[ nodeNum ].x;
			elementNodalValues[1][i][n] = displacements[ nodeNum ].y;
			elementNodalValues[2][i][n] = displacements[ nodeNum ].z;
			}// n loop ( loop over nodes )
}// i loop (= element #)

copyData();   // This copies from contour database to element structure
BETA_OUT<<"convertNodalSetToElementSet...exit"<<endl;
}
////////////////////////////////////////////////////////////////////////////

//==========================================================================
void ContourData::normalizeDataValues(//	float *MinValue, float *MaxValue,
              	bool fixed_option,bool group_option,bool active_option)
{
	ElementGroup &element=mesh->element;
int numElem = mesh->numElements;

	int n,currentGroup,materialGroup,nnpe,i,j,elementNumber,maxLocation,minLocation;
	double overallMin,overallMax,minActive,maxActive;
	double minGroup,maxGroup,fixedMin=0.0,fixedMax=0.0,minNormal,maxNormal;;
	bool atLeastOne;
	// Now loop through all elements and determine max and min values for
	// all groups

	if(MaxNumberOfMaterialsForContouring==0){
		for( elementNumber =0; elementNumber< numElem; elementNumber++) {
			materialGroup = element[elementNumber].getMaterialNumber();
			if(materialGroup>MaxNumberOfMaterialsForContouring)
				MaxNumberOfMaterialsForContouring=materialGroup;
		}
		MaxNumberOfMaterialsForContouring++;
	}
	
	overallMin =  1.e20f;
	overallMax = -1.e20f;
	for(currentGroup = 0; currentGroup< MaxNumberOfMaterialsForContouring; currentGroup++) {
		atLeastOne = false;
		//if(verbose) printf("cg = %d\n",currentGroup);
		MinValue =  1.e20f;
		MaxValue = -1.e20f;
		for( elementNumber =0; elementNumber< numElem; elementNumber++) {
			//materialGroup = convertMaterialGroupToColorGroup(element[elementNumber].matGroupNumber);
			materialGroup = element[elementNumber].getMaterialNumber();
			if(materialGroup==currentGroup) {
				atLeastOne = true;
				nnpe = element[elementNumber].numNodesPerElement;
				for(i = 0; i<nnpe;i++) {
					if( element[elementNumber].nodalValues[i] < MinValue) {
						MinValue = element[elementNumber].nodalValues[i];
						if(group_option && false) //materialGroup==select_group)
							minGroup = MinValue;
						minLocation = elementNumber;
					}
					if( element[elementNumber].nodalValues[i] > MaxValue) {
						MaxValue = element[elementNumber].nodalValues[i];
						if(group_option && false) //materialGroup==select_group)
							maxGroup = MaxValue;
						maxLocation = elementNumber;
					}
				}
			}
		}
		
		if(atLeastOne) {
			if(MinValue<overallMin)
				overallMin = MinValue;
			if(MaxValue>overallMax)
				overallMax = MaxValue;
			//if(verbose) printf("Min and max values for group %d = %e %e\n",
			//	currentGroup,MinValue,MaxValue);
			//if(verbose) printf("Elements for min and max are %d %d \n",
			//	minLocation+1,maxLocation+1); 
		}
	}

	minActive = 1e20f;
	maxActive = -1e20f;

int elemMax, elemMin;

	for(i=0;i<numElem;i++)
		if(element[i].activeElementFlag==ACTIVE) {
			nnpe = element[i].numNodesPerElement;
			for(j=0;j<nnpe;j++) {
				if(element[i].nodalValues[j] > maxActive)
					maxActive = element[i].nodalValues[j]; elemMax = i;
				if(element[i].nodalValues[j] < minActive)
					minActive = element[i].nodalValues[j]; elemMin = i;
			}
		}
	
//	if(verbose) printf("For entire mesh, the min. and max. values are %e %e\n",
//	overallMin, overallMax);
	MinValue = overallMin; 
	MaxValue = overallMax;	
	
	/************** Override Auto Scaling *******************/
	if (fixed_option == true)	{
		MinValue = fixedMin;
		MaxValue = fixedMax;
	}
	if (group_option && false)	{
		MinValue = minGroup;
		MaxValue = maxGroup;
	}
	if (active_option)	{
//		if(verbose) printf("\nAuto Scaling Override = active mode\n");
		MinValue = minActive;
		MaxValue = maxActive;
	}

// Check for valid limits.
	if(MinValue>MaxValue)	{
		BETA_OUT << "Error: Reading element nodal values.\n\tMinValue >= MaxValue\n";
		BETA_OUT << "MinValue = " << MinValue << "\newMaxValue = " << MaxValue << '\n';
		exit(1);
	} 

	//vout << "\nScaling data to [" << MinValue << ',' << MaxValue << "]\n";

	/*******************************************************/
	
	// Need logic to determine whether to normalize all groups the same or 
	// individually.
	
	// Calculate normalizing factors  (see function normalize)
	//
	//   finish
	minNormal = 0;
	maxNormal = 12; // changed from 11 to 12 cdc
    double  deltaValue;
	deltaValue = MaxValue-MinValue;
	if(deltaValue==0) {deltaValue=1.e-10;}
	double a,b;
	a = ((minNormal*MaxValue)-(maxNormal*MinValue))/(deltaValue);
	b = (maxNormal - minNormal)/(deltaValue);
	
	// Loop through elements and normalize the nodal values
	//
	//  finish
	for(i=0; i< numElem; i++)
		for(n = 0; n<element[i].numNodesPerElement;n++) 
			element[i].nodalValues[n] = a + b * element[i].nodalValues[n];
}

///////////////////////////////////////////////////////////////
void ContourData::rescale(double newMaxValue,double newMinValue)
{
	double minNormal,maxNormal;
	double a,b;
	int i,n;
	int numElems;
	//double *maxValue = ;
    //double *minValue = ;

	ElementGroup &element=mesh->element;
	numElems = mesh->numElements;


	minNormal = 0;
	maxNormal = 12;
	double  deltaValue;
	deltaValue = newMaxValue-newMinValue;
	if(deltaValue==0) {deltaValue=1.e-10;}

	a = ((minNormal*newMaxValue)-(maxNormal*newMinValue))/(deltaValue);
	b = (maxNormal - minNormal)/(deltaValue);
	
	
	// Loop through elements and normalize the nodal values
	//
	//  finish
	for(i=0; i< numElems; i++)
	{ 
		for(n = 0; n<element[i].numNodesPerElement;n++)   
		{
		// Recover old Value
			element[i].nodalValues[n]= MinValue+element[i].nodalValues[n]* 
				(MaxValue - MinValue)/(maxNormal-minNormal);
		// Scale to new Value
			element[i].nodalValues[n] = a + b * element[i].nodalValues[n];
		}
	}
// Store values in PlotInterface object
//	mesh->maxVal = newMaxValue;  //figure out what to do about this!!!!
//	mesh->minVal = newMinValue;  //figure out what to do about this!!!!
//	if(verbose) printf("Done Rescaling\n");

MaxValue = newMaxValue;
MinValue = newMinValue;
}
//=============================================================================

void ContourData::print()
{


}

//===========================================================
void ContourData::setFilename(string name)
{
//BETA_OUT<<"(setFilename)contourDataFilename = "<<filename<<endl;
filename=name;
BETA_OUT<<"(setFilename)contourDataFilename = "<<filename<<endl;

binaryFilename=filename + ".binary";
//BETA_OUT<<"contourDataFilename = "<<filename<<endl;
//BETA_OUT<<"(setFilename)contourDataFilename(binary) = "<<binaryFilename<<endl; 
};

//==========================================================

void ContourData::allocateMemory()
{
if(verboseFlag > Basic) displayTime("ContourData::allocateMemory-begin", BETA_OUT);

int i,j;
ElementGroup &element=mesh->element;
int	numElem = mesh->numElements;

freeMemory();
numElemAlloc = numElem;
numColumnsAlloc = numColumns;

//for(i=0; i< numElem; i++) {
//element[i].nodalValues = (double *)malloc(element[i].numNodesPerElement * sizeof(double));
//};


// Memory allocation for elementNodalValues database
elementNodalValues = (double ***) malloc( numColumns * sizeof(double**) );
for(i=0; i<numColumns; i++){
	elementNodalValues[i] = (double **) malloc( numElem * sizeof(double*) );
    if(elementNodalValues[i]==NULL) EXIT_BETA("memory allocation failure");
	  for(j=0; j<numElem; j++){
		elementNodalValues[i][j]=(double *)malloc( element[j].numNodesPerElement * 
                                               sizeof(double));
        if(elementNodalValues[i][j]==NULL) EXIT_BETA("memory allocation failure");
	  } // j
} // i
elementNodalValuesState = true;

if(verboseFlag > Basic) displayTime("ContourData::allocateMemory-exit", BETA_OUT);
};

//==========================================================

void ContourData::freeMemory()
{

	BETA_OUT<<"freeMemory"<<endl;
if(elementNodalValuesState == false)
 {BETA_OUT<<"already released"<<endl;return;}

int i,j;

//BasicElement *element = mesh->element;

// Element object contour data    ... this really should be handled by element class
//for(i=0; i< numElemAlloc; i++) {
//		free(element[i].nodalValues); element[i].nodalValues=NULL;};

// ElementNodalValues database

if(elementNodalValuesState == true){
for(i=0; i<numColumnsAlloc; i++){
	for(j=0; j<numElemAlloc; j++){
		free(elementNodalValues[i][j]);
      }
	 free(elementNodalValues[i]);
	}
free(elementNodalValues);  //JV093002 fixed leak
elementNodalValues=NULL;
};//end of if

elementNodalValuesState = false;
numElemAlloc    = 0;
numColumnsAlloc = 0;
};
