#include "stdafx.h"

#include "utility_defines.h"
#include "utility.h"
#include <string.h>
#include<sstream>
using namespace std;

#define STATIC_DATA static
//============================
string peekAtNextLine(istream *inputStream)
//JV040609: noticed an issue with this function when the plt file has been created in UNIX environment and is being read in the Windows environment
{
	// Check content of next line w/o advancing position in file
	streampos m_posStartData;
	
	m_posStartData = inputStream->tellg(); //save current location
//	cout<<"m_posStartData ="<<m_posStartData<<endl;
	
	
	string wholeLine;
	getNextLine:  getline((*inputStream),wholeLine);
	if(wholeLine=="") goto getNextLine;
//	cout<<"wholeLine = "<<wholeLine<<endl;
	
	inputStream->seekg( m_posStartData ); //restore original location
//	cout<<"m_posStartData ="<<m_posStartData<<endl;
	
	return(wholeLine);
}//
//=============================================================
void tokenize( char *string,char **tokenList, int &number)
{
 char delimiters[] = " \t";
		//BETA_OUT<<"String ="<<string<<endl;
		tokenList[0] = strtok(string, delimiters);
int i;

		for(i=1; i<30; i++)
		{
		  tokenList[i] = strtok( NULL, delimiters);
		  if(tokenList[i] == NULL) break;
			//BETA_OUT<<tokenList[i]<<endl;
		}
		number = i;
		if(tokenList[0] == NULL) number = 0;
		  //BETA_OUT<<"I found "<<number<<" tokens"<<endl;
 }
//=============================================================
void tokenize( char *string,char **tokenList, int &number, char* moreDelimiters)
{
		char delimiters[] = " \t";
		strcat(delimiters,moreDelimiters);

		//OUT<<"String ="<<string<<endl;
		tokenList[0] = strtok(string, delimiters);

		int i;
		for(i=1; i<30; i++)
		{
		  tokenList[i] = strtok( NULL, delimiters);
		  if(tokenList[i] == NULL) break;
			//OUT<<tokenList[i]<<endl;
		}
		number = i;
		if(tokenList[0] == NULL) number = 0;
		  //OUT<<"I found "<<number<<" tokens"<<endl;
 }
//=============================================================
// New version of getLineAndTokenize: xtang's version!!!

//==========================================================
// xtang 991103: 
// Note: The function "strpbrk" in VC++ 6.0 does not work correctly.
//       It does not return a null pointer which should be the case when there
//       is no match available between the two strings.

size_t _strlen (const char * str)
{
	int length = 0;
	if (str !=0) while( *str++ ) ++length;
	return( length );
}

char * _strpbrk(char *inputString, char * ppp)
{
	int i=0,j;
	char * p;
	bool found=false;

	do {
		if ((inputString[i] == ppp[0])) {
			found=true;
			for (j=0; j < (int) _strlen(ppp); j++) {found &=(inputString[i+j]==ppp[j]);}
			if (!found) {i+=j;} 
			else {p=inputString+i;}
		}
		else i++;
	} while ((!found) && (i < (int) _strlen(inputString)));
	if (found) return p;
	else return 0;
}

char * effectstrpbrk(char *inputString, char * ppp)
{
	char  *p=inputString;
	int i,j, pL, pR;
	bool fL,fR;
top:	
	fL=false;
	fR=false;
	p = _strpbrk(p,ppp);
	if (p==0) { return 0;}
	pL = p-inputString;
	if (pL==0) { return p;}
	pR = (inputString+_strlen(inputString))-(p+_strlen(ppp));
	if (pR <= 0) { return p;}

	i=pL;
	do { 
		i--;
		if (inputString[i]=='\"') {fL = !fL;}
	} while ((i > 0)); 
	if (fL) {
		i=j=_strlen(ppp);
		do { 
			if (p[i]=='\"') { fR = true;}
			i++;
		} while ((i < (j+pR) ) && (!fR)); 
	}
	if ( fL && fR && (i < (j+pR)) ) { p+=i; goto top;}
	return p;
}

void cutString(char *p1, char *p2, char *str)
{
	char ns[400];
	int p=(p1-str);;
	strncpy(ns,str,p);
	ns[p]=0;  
	strcat(ns,p2);
	strcpy(str,ns);
}

char *stripCStyleComment(char * inputString)
{
	char *p1, *p2;
	bool ok;
top:
	ok=false;
	p1=effectstrpbrk(inputString,"/*");  
	p2=effectstrpbrk(inputString,"//");
	if (_strlen(p1) > _strlen(p2)) {
		p2=effectstrpbrk(inputString,"*/");
		if (p2 !=0) {
			p2+=2;
			cutString(p1,p2,inputString);
			ok=true;
		}
	}
	if (ok) goto top;
	if (_strlen(p2)>_strlen(p1)) { return p2; }
	else if (_strlen(p1)>_strlen(p2)){ return p1;}
	return 0;
}

bool moveToNextEffectLine(istream * inStream)
{
	char in[401], * sp;
    int iii;
	do {
		(*inStream).getline(in, 400,'\n');
  	     // xtang 
         iii=(int)_strlen(in);
	     if (in[iii-1]==13) {in[iii-1]=0;}

		sp=effectstrpbrk(in,"*/");
	} while ((sp==0) && (!(*inStream).eof()));

	// xtang 030101
	if (sp!=0) {
		//(*inStream).seekg(-_strlen(sp)+1,ios::cur);  // put pointer to the new effective location
		(*inStream).seekg(-((int)_strlen(sp)),ios::cur);  // Xytang 030201 xxxx?
		return true;
	}
	else return false;
}

bool process_SSSS(istream * inStream, char * inputString)
{
	// case 1:  abcd // efgh,  return  abcd
	// case 2:  //  abcd,      return  null string
	// case 3:  abcd /* efgh,  return  abcd and move pointer after paired */
	// case 4:  /* abce,       return null string and move pointer after paired */
	// case 5   abcd /* efgh */ ijkl,  return abcd ijkl;

	char *sp;
	bool ok=true;

	sp=stripCStyleComment(inputString);

	// xtang 030101
	if ((sp!=0)&&(sp[1]=='*')) ok=moveToNextEffectLine(inStream); // needs to move pointer
	if (sp !=0 ) {
		cutString(sp,inputString+_strlen(inputString),inputString);
	}
	return ok;
}

int getLineAndTokenize(istream * inStream, char *terminateString,
	char   *tokenList[], int &numberOfTokens, char * inputString)
{ 
	const int lLen=300;
	
	bool hasT;
top: 
	if((*inStream).eof() )  return(1);
	(*inStream).getline(inputString, lLen,'\n');
    int iii=(int)_strlen(inputString);
	if (inputString[iii-1]==13) {inputString[iii-1]=0;}

	hasT=process_SSSS(inStream, inputString);
	tokenize( inputString, tokenList, numberOfTokens);
	if(numberOfTokens==0) goto top; 

	if (!hasT) {
		cout<<"Unpaired Comment symbol \"/*\" !!!"<<endl;
		return (2);
	}

	if(COMPARE(tokenList[0],terminateString)==0)  {   
		cout<<"TerminateString: "<<tokenList[0]<<endl;
		return(1);
	}

	return(0);
}
    
int getLineAndTokenize(istream * inStream, char *terminateString,
	char   *tokenList[], int &numberOfTokens)
{ 
	const int lLen=300;
	STATIC_DATA char inputString[lLen+1];
    
	
	bool hasT;
top: 
	if((*inStream).eof() )  return(1);
	(*inStream).getline(inputString, lLen,'\n');
    int iii=(int)_strlen(inputString);
	if (inputString[iii-1]==13) {inputString[iii-1]=0;}
	hasT=process_SSSS(inStream, inputString);
	tokenize( inputString, tokenList, numberOfTokens);
	if(numberOfTokens==0) goto top; 

	if (!hasT) {
		cout<<"Unpaired Comment symbol \"/*\" !!!"<<endl;
		return (2);
	}

	if(COMPARE(tokenList[0],terminateString)==0)  {   
		cout<<"TerminateString: "<<tokenList[0]<<endl;
		return(1);
	}

	return(0);
}

//==========================================================

vector<string> TokenizeStringToVectorOfStrings(string &stringToTokenize){
    vector<string> result(0);
    istringstream ss(stringToTokenize);
    string token;
    while(ss >> token){
        result.push_back(token);
    }
    return result;
}
//==========================================================
vector<int> TokenizeStringToVectorOfIntegers(string &stringToTokenize){
    vector<int> result(0);
    istringstream ss(stringToTokenize);
    int token;
    while(ss >> token){
        result.push_back(token);
    }
    return result;
}
//==========================================================
vector<double> TokenizeStringToVectorOfDoubles(string &stringToTokenize){
    vector<double> result(0);
    istringstream ss(stringToTokenize);
    double token;
    while(ss >> token){
        result.push_back(token);
    }
    return result;
}
//==========================================================


//======================
















/*  xtang: 012401: old JDW version



int getLineAndTokenize(istream * inStream, char *terminateString,
							  char   *tokenList[], int &numberOfTokens)
{ int i;

// Future mod: comparison of // and token should only use first
// two characters (or less) of token

  //char *inputString = tokenSpace;
 
// xtang: 012301: If one wants to read long string, one should change the 
// default length of string from [80] to a longer string, say [300] in xtang's version

  STATIC_DATA char inputString[80];

top: if((*inStream).eof() )  return(1);
     (*inStream).getline(inputString, 80,'\n');
	 tokenize( inputString, tokenList, numberOfTokens);

	  if(numberOfTokens==0)goto top;
	  if(COMPARE(tokenList[0],terminateString)==0)
			{BETA_OUT<<"TerminateString: "<<tokenList[0]<<endl;
			 return(1);
			}
     if((COMPARE(tokenList[0],"/ * ")==0) )
	    {
skip:   (*inStream).getline(inputString, 80,'\n');
        if((COMPARE(tokenList[0],"* /")==0) )
		  {goto top;}
        goto skip;
	    }

     if((COMPARE(tokenList[0],"//")==0) ){ goto top;}

	  if(strncmp(tokenList[0],"//",2)==0) { goto top;}

     //BETA_OUT<<"numberOfTokens= "<<numberOfTokens<<endl;
	  //for(int i=0; i<numberOfTokens; i++){
	//	  BETA_OUT<<tokenList[i]<<" ";} BETA_OUT<<endl;

     for(i=0; i<numberOfTokens; i++)
	  {if((COMPARE(tokenList[i],"//")==0) )		 
	  {//BETA_OUT<<"Reduced # of tokens from "<<numberOfTokens
	    //   <<" to "<<i<<endl;
		  numberOfTokens = i; return(0);}
	  }
	  return(0);
}
*/
//==========================================================
