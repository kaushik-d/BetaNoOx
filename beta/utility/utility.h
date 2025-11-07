#ifndef utility_h
#define utility_h

#include "utility_defines.h"

#include "array.hpp"
#include<string>
#include<vector>

using namespace std;

//void banner( char * string);
void banner( char * string, ostream &outStream);

int getLineAndTokenize(istream *inStream,  char *terminateString, char *tokenList[], int  &numberOfTokens);
int getLineAndTokenize(istream * inStream, char *terminateString, char *tokenList[], int &numberOfTokens, char * inputString);
vector <string> TokenizeStringToVectorOfStrings(string &stringToTokenize);
vector <double> TokenizeStringToVectorOfDoubles(string &stringToTokenize);
vector <int> TokenizeStringToVectorOfIntegers(string &stringToTokenize);

void tokenize( char *string,char **tokenList, int &number);
void tokenize( char *string,char **tokenList, int &number, char* moreDelimiters);



string peekAtNextLine(istream *inputStream);
ifstream *openFileForInput(char * filename);

void displayTime(const char * message, ostream &outStream);

bool doesFileExist(string fileName);
bool doesDirExist(const char * fileName);
bool makeDir(const char * dirPath);

void releaseGaussQuadratureStaticMemory();
void releaseExtrapolationStaticMemory();

//Math Utilities
double dotProduct( const double * a, const double * b, const int &length);
void copyVector( const double * copyFrom, double * copyTo, int length);
void sumVector( double *c, double * a, double * b, int length);
void initializeVector( double * a, double value, int length);
void ConvertNumericalZeroToZero(double *a, int length, double numericalZero);
int ConvertNegativeNumbersToZero(double *a, int length);

void WriteDoubleVectorToFile(string filename, double *list, int L);
bool ReadDoubleVectorFromFile(string filename, double *&list, int L);

enum VerboseLevel {Off, Basic, Min, Max};
enum BETA_Time_Unit {seconds, minutes, hours};

#endif
