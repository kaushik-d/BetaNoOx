#include "stdafx.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
using namespace std;
typedef double REAL;

//===============================================
char *CutEZero(char *sNumStr)
{
    int len,i;
    int k=0, p=0;
    char *sDest;
	sDest= sNumStr;
       
    len=strlen(sNumStr);
    for (i=0; i<len; i++)
      if ((sNumStr[i]=='e') ||(sNumStr[i]=='E')) {
        k++;
        p=i;
    }
    p=len-p;
    if ((k==1) && (p >3)) { 
        i=k=0;
        while (i<=len) {
            sDest[k]=sNumStr[i];
            if ((sNumStr[i]=='e') ||(sNumStr[i]=='E')) {
               i++;
               k++;
               sDest[k]=sNumStr[i];
               i+=(p-4); // 100197-tang
            }
        i++;
        k++;
        }
    }
    else sDest=sNumStr;
    return sDest;
}
//====================================================================
/*
void printRow( int n, int stride, REAL * list, int width, ofstream *stream)
{int index = 0;
 (*stream)<<setiosflags(ios::scientific);
 for(int i =0; i<n; i++)
  {(*stream)<<setw(width)<<list[index];
	index += stride;
  } ; (*stream)<<endl;
}

//====================================================================

void printRow( int n, int stride, REAL * list, int width, 
				  int precision, ostream *stream)
{int index = 0;
 (*stream)<<setiosflags(ios::scientific);
 for(int i =0; i<n; i++)
  {(*stream)<<setw(width)<<setprecision(precision)<<list[index];
	index += stride;
  } ; (*stream)<<endl;
}


//====================================================================


void printRow( int n, int stride, REAL * list, int width)
{int index = 0;
 OS<<setiosflags(ios::scientific);
 for(int i =0; i<n; i++)
  {OS<<setw(width)<<list[index];
	index += stride;
  } ; OS<<endl;
}

//====================================================================
void printTable( int numRows, int numCols, REAL * list, int width)
{
 int row,  j;
 for( row=0; row< numRows; row++)
	  {
		for(j=0; j<numCols; j++)
			{
			 OS<<setw(width)<<list[ row + j * numRows ];
			} ; OS<<endl;
	  }
}
//====================================================================

void spaces()
{
OS<<endl<<endl;
}
//====================================================================
*/