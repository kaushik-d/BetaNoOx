#include "stdafx.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include "sparseMatrixInterface.hpp"

CreateErrorHandler(SparseMatrixInterface);

SparseMatrixInterface::SparseMatrixInterface() 
{
	type=Undefined;
}

SparseMatrixInterface::~SparseMatrixInterface()
{
	DeleteWhoMethod("~SparseMatrixInterface()");
//	cerr << "Deleting SparseMatrixInterface" << endl;
	int i;
	const int L=numEquations;

	if(L>0) {
		for(i=0;i<L;i++) 
			delete [] mat[i];
		delete [] mat;
		delete [] row;
		delete [] col;
	}
}

char* SparseMatrixInterface::getType(int i)
{
	switch(i){
		case 0: return("Undefined");break;
		case 1: return("Lower Triangular Matrix");break;
		case 2: return("Upper Triangular Matrix");break;
		case 3: return("Unsymmetric Matrix");break;
		default : cout << "type unrecognized" << endl; exit(1);
	}
}

bool SparseMatrixInterface::ReadMatrix(char* filename)
{
	bool answer=false;
	ifstream is;
	is.open(filename);
	answer=ReadMatrix(&is);
	is.close();
	return answer;
}
bool SparseMatrixInterface::ReadMatrix(istream* instream)
{
	int i;
	//reads the RCV format file
	(*instream) >> i;
	if(i>3 || i<1){
		cerr << "incorrect format for globalstiffness.txt" << endl;
		exit(1);
	}

	type=(MatrixType)i;
	cout << "matrix type = " << getType(type) << endl;
	(*instream) >> numEquations >> numNonZeros;
	const int L=numEquations;

	Assert(row = new SortedList<int>[L]);
	Assert(col = new SortedList<int>[L]);

	//store the position of the stream buffer
	streampos m_posStartData;
	m_posStartData = instream->tellg(); 

	int ii,jj;
	double value;

	for(i=0;i<numNonZeros;i++) {
		(*instream) >> ii >> jj >> value;
		specifyNonZeroLocation(ii,jj);
	}

	//allocate memory
	Assert(mat = new double *[L]);

	for(i=0;i<L;i++) {
		if(row[i].getNum()==0) row[i].add(i);
		if(col[i].getNum()==0) col[i].add(i);
		Assert(mat[i] = new double [row[i].getNum()]);

	}

	//rewind back to the start of RCV values
	instream->seekg(m_posStartData);

	//store the actual values into the matrix
	for(i=0;i<numNonZeros;i++) {
		*instream >> ii >> jj >> value;
		operator()(ii,jj) = value;
	}

	return true;
}

void SparseMatrixInterface::specifyNonZeroLocation(int ii, int jj)
{
	row[ii].add(jj);
	col[jj].add(ii);

}

double &SparseMatrixInterface::operator()(const int ii,
		const int jj)
{
	int j;
	if(ii>=numEquations || jj >=numEquations ) {
//_operatorError:
		cerr << "numberOfEquations = " << numEquations << endl;
		WhoMethod("operator()(const unsigned Int,const unsigned Int)");
		cerr << "ii,jj = " << ii << " , " << jj << endl;
		FatalError("Bad operand...");
	}
	j = row[ii].find(jj);
	if(j==row[ii].getNum()) {
		cerr << jj << " not on row " << ii << endl;
		cerr << "Row " << ii << "-" << row[ii].getNum();
		row[ii].print(cerr) << endl;
		FatalError("Bad operand...");
	}
	return mat[ii][j];
}
