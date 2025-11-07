#include "stdafx.h"

#include <stdlib.h>
#include <cmath>
#include <string.h>
#include "math/matrix.hpp"
#include "linearEquation.hpp"
#include "equation.hpp"
#include "utility/utility_defines.h"
#include "utility/excepts.hpp"
#include "utility/formWriter.hpp"
#undef PRINT
//#define PRINT(a) cout << __LINE__ << ": " << #a << " = " << a << endl;
#define PRINT(a)

CreateErrorHandler(Equations);

extern "C" {void MKLGetVersionString(char * buffer, int len);}
extern "C" {int omp_get_max_threads();}
//====================================
//handle solver/platform checks
#if !defined(WIN32) && defined(OLAF)
	#undef OLAF
#endif
//====================================
#define FORMAT setw(30)<<setprecision(17)
//====================================

void Equations::setStorageMethod(istream * inStream)
{
	char method[80];
	(*inStream) >> method;
	setStorageMethod(method);
}

void Equations::setStorageMethod(char * method)
{
char defaultmethod[80];
StorageMethod defaultsm=sm;
if (defaultsm==Sparse) strcpy(defaultmethod,"SPARSE");
if (defaultsm==Profile) strcpy(defaultmethod,"PROFILE");
if (defaultsm==Osparse) strcpy(defaultmethod,"OLAF");
if (defaultsm==mklgbmat) strcpy(defaultmethod,"MKLGBMAT");
if (defaultsm==LapackGB) strcpy(defaultmethod,"LAPACKGB");
if (defaultsm==PetSc) strcpy(defaultmethod,"PETSC");
if (defaultsm==PetScUnsymm) strcpy(defaultmethod,"PETSCUNSYMM");
if (defaultsm==MKLPardiso) strcpy(defaultmethod,"MKLPARDISO");
if (defaultsm==MKLPardisoUnsymm) strcpy(defaultmethod,"MKLPARDISOUNSYMM");
if (defaultsm==MKLCG) strcpy(defaultmethod,"MKLCG");
if (defaultsm==ParCG) strcpy(defaultmethod,"PARCG");
if (defaultsm==Pardiso) strcpy(defaultmethod,"PARDISO");
if (defaultsm==WSMPSymm) strcpy(defaultmethod,"WSMPSYMM");

//	(*outStream)<<"Storage method="<<method<<endl;
	int len=198;
	char buf[198];
#ifdef MKL
	MKLGetVersionString(buf, len);
#endif

#ifdef SPARSE
	if(COMPARE(method, "SPARSE")==0){
		(*outStream) << "Storage Method: SPARSE" << endl;
		sm =Sparse;} 
	else
#endif
#ifdef PROFILE
	if(COMPARE(method, "PROFILE")==0){
		(*outStream) << "Storage Method: PROFILE" << endl;
		sm =Profile;}
	else
#endif
#ifdef OLAF
	if(COMPARE(method, "OLAF")==0){
		(*outStream) << "Storage Method: OLAF" << endl;
		sm =Osparse;}
	else
#endif
#ifdef LAPACKGB
	if(COMPARE(method, "LAPACKGB")==0){
		(*outStream) << "Storage Method: LAPACKGB" << endl;
		sm =LapackGB;}
	else
#endif
#ifdef PETSC
	if(COMPARE(method, "PETSC")==0){
		(*outStream) << "Storage Method: PETSC" << endl;
		sm =PetSc;}
	else
	if(COMPARE(method, "PETSCUNSYMM")==0){
		(*outStream) << "Storage Method: PETSCUNSYMM" << endl;
		sm =PetScUnsymm;}
	else
#endif
#ifdef MKL
	if(COMPARE(method, "MKLGBMAT")==0){
		(*outStream) << "Storage Method: MKLGBMAT" << endl;
		(*outStream) << "MKL Version : " << buf << endl;
		sm =mklgbmat;}
	else
	if(COMPARE(method, "MKLPARDISO")==0){
		(*outStream) << "Storage Method: MKLPARDISO" << endl;
		(*outStream) << "MKL Version : " << buf << endl;
		sm =MKLPardiso;}
	else
	if(COMPARE(method, "MKLPARDISOUNSYMM")==0){
		(*outStream) << "Storage Method: MKLPARDISOUNSYMM" << endl;
		(*outStream) << "MKL Version : " << buf << endl;
		sm =MKLPardisoUnsymm;}
	else
	if(COMPARE(method, "MKLCG")==0){
		(*outStream) << "Storage Method: MKLCG" << endl;
		(*outStream) << "MKL Version : " << buf << endl;
		sm =MKLCG;}
	else
#endif
#ifdef PARCG
	if(COMPARE(method, "PARCG")==0){
		(*outStream) << "Storage Method: PARCG" << endl;
		sm =ParCG;}
	else
#endif
#ifdef PARDISO
	if(COMPARE(method, "PARDISO")==0){
		(*outStream) << "Storage Method: PARDISO" << endl;
		sm =Pardiso;}
	else
#endif
#ifdef WSMP
	if(COMPARE(method, "WSMPSYMM")==0){
		(*outStream) << "Storage Method: WSMPSYMM" << endl;
		sm =WSMPSymm;}
	else
#endif
	{
		(*outStream)<<"Solver Method: " << method << " selected does not exist for this platform."<<endl;
		exit(1);
	}

}
/////////////////////////////////////////////////////////////////////////////
//
//  S E T   A N D   G E T   E Q U A T I O N   I N F O . . .  
//
bool Equations::readMPCNeutralFile(char *filename)
{
	ifstream inStream;
	inStream.open(filename);
	if(!inStream.good()) return false;
	bool result=true;

	AdditionalEquation *eqn;

	
	int					numberOfTokens;
	char				*tokenList[20];
	int					exitFlag;

    while (true) {
		exitFlag=getLineAndTokenize(&inStream,"exitReadMultiPointConstraints",
			                         tokenList, numberOfTokens);
		if(exitFlag==1) goto _EndreadMPCNeutralFile;

		if(numberOfTokens%2!=0){
			cout <<"incorrect additional equations" << endl;
			result=false;
			goto _EndreadMPCNeutralFile;
		}

		eqn = new AdditionalEquation;
		eqn->setSlave( atoi(tokenList[0]) );

		for (int i=0; i<numberOfTokens/2 - 1; i++) {
			int mdof		=atoi(tokenList[2*i + 1]);
			double mfactor	=atof(tokenList[2*i+1 + 1]);
			eqn->addAdditionalTerm( mdof, mfactor);
		}

		eqn->setConstant(atof(tokenList[numberOfTokens-1]));
		if(!addAdditionalEquation(eqn)){
			delete eqn;
			cout <<"bad equation" << endl;
			exit(1);
		}

	}
_EndreadMPCNeutralFile:
	inStream.close();
	return result;
}
//=============================================================
bool Equations::addAdditionalEquation(AdditionalEquation *eqnToAdd)
{
	WhoMethod("addAdditionalEquation(AdditionalEquation *)");
	if(allocated) 
		FatalError("No additional equations may added after storage has been allocated.");
	if(eqnToAdd->slave < 0) {
		cerr << "Slave cannot be negative " << eqnToAdd->slave << '\n';
		FatalError("Bad slave unknown.");
	}
	return addAdditionalEquation_Kernel(eqnToAdd);
}
//=============================================================
bool Equations::addAdditionalEquation_Kernel(AdditionalEquation *eqnToAdd)
{
	WhoMethod("addAdditionalEquation_Kernel(AdditionalEquation *)");

//the logic for this function is explained very well in Clint Chapman's dissertation page 26 and 27

// Procedure for adding a new equation:
// 1. Replace all dof on rhs which are already slaved.
// 2. If new slave dof is already a slave, switch with one of the masters on rhs
// 2a. Replace the slave dof which is now on rhs
// 3. This equation should be unique.
// 4. Replace the new slave dof in all existing equations.
	
	//consolidates the addnl. equation to the form : slaveDOF= [masterDOF*masterFactor] + constant
	AdditionalEquation::Result result;
	PRINT(*eqnToAdd);
	result = eqnToAdd->collectTerms();
	switch(result) {
		case (AdditionalEquation::Good):
			break;
		case (AdditionalEquation::Bad):
			cerr << *eqnToAdd << " is impossible to satisfy." << endl;
			Warning("Bad equation.");
			return false;
		case (AdditionalEquation::Trivial):
			delete eqnToAdd;
			return true;
		default:
			FatalError("Don't understand AdditionalEquation::Result");
			break;
	}

	//step 1: (from pg 26 in chapman dissertation)
	if(eqnToAdd->numTerms > 0) {
		if(aeList.begin()) {
			do {
				result = eqnToAdd->apply(*aeList());
				PRINT(*eqnToAdd);
				PRINT(result);
				PRINT(aeList.AreYouEmpty());
				switch(result) {
					case (AdditionalEquation::Good):
						break;
					case (AdditionalEquation::Bad):
						cerr << *eqnToAdd << " is impossible to satisfy." << endl;
						Warning("Bad equation.");
						return false;
					case (AdditionalEquation::Trivial):
						delete eqnToAdd;
						return true;
					default:
						FatalError("Don't understand AdditionalEquation::Result");
						break;
				}
			} while (aeList++);
		}
	}

	//step 2: (from pg 26 in chapman dissertation)
	if(aeList.begin()) {
		do {
			result = aeList()->apply(*eqnToAdd);
			switch(result) {
				case (AdditionalEquation::Good):
					break;
				case (AdditionalEquation::Bad):
					cerr << "New Equation:      " << *eqnToAdd << '\n';
					cerr << "Existing Equation: " << *aeList() << '\n';
					Warning("Impossible to satisfy together.");
					return false;
				case (AdditionalEquation::Trivial):
					delete eqnToAdd;
					return true;
				default:
					FatalError("Don't understand AdditionalEquation::Result");
					break;
			}
		} while (aeList++);
	}

	//add resulting unique equation to existing list of equations
	if(!aeList.insert(eqnToAdd))
		FatalError("Cannot add AdditionEquation to aeList!");
	return true;
}

// xtang 03102003
void Equations::restrainHangingDOF()
{
   int i;
   for(i=0;i<numberOfEquations;i++) 
   {
	   cout <<i <<"\t"<<(*a)(i,i)<<endl;
		if(fabs((*a)(i,i)) <1e-5) 
			setRestraint(i,true);
   }

}
void Equations::emptyRestraintList()
{
	restraintList.DestroyList();
}
void Equations::setRestraint(const int eqnNum, const bool val)
{
	//JV070608 : if memory has been allocated, there is no use of the 'restrainList'
	if(allocated) 
		restraint[eqnNum] = val;
	else{
		if(val) {
			if(restraintList.find(eqnNum)==0) restraintList.insert(new int(eqnNum));
		} else {
			restraintList.erase(eqnNum);
		}
	}
}

void Equations::setSolution(const double &val, const int eqnNum)
{
	// A problem exists here. If one wants to set the value of a slaved
	// unknown, what should be done?
	// 1. Make that slave a master. This requires making an unrestrained master a slave.
	// 2. The new slave equation should then be applied to all existing additional eqns.
	WhoMethod("setSolution(const double &, const int)");
	if(!allocated) FatalError("Preliminary allocation must be performed before a value can be set.");
	RangeCheck(0,eqnNum,numberOfEquations-1);
	if(restraint[eqnNum]) {
		if(additionalEquation[eqnNum]) 
		{cout<<"Equation# "<<eqnNum<<" is restrained"<<endl;
         cout<<" I think there is an additional equation for this dof...???"<<endl;
         cout<<" JV060805 - Be sure you understand what this is doing. check if this was expected or OK for your model. "<<endl;
         cout<<" JV060805 - It might be expected for some hierarchical models. "<<endl;
//JV063009 commented to handle hierarchical models
//			FatalError("Cannot set value of slaved unknown this way. Use additionalEquations.");
		}
	}
	restraint[eqnNum] = true;
	solution[eqnNum] = val;
}

void Equations::updatePseudoLoadVectorUsingGK(const double &val, const int eqnNum, LargeMatrix *gK)
//JV101805 :this function was written initially when generating multi-field macro element stiffness matrices
// IF multiple models (SAME MESH) need to be solved AND the restrained dofs do not change (their values may change)
// then this function can be used to avoid factorizing the global matrix each time. gK has to be an assembled global K that 
// does not have any restraints. see documentation for explanation
{
	WhoMethod("updatePseudoLoadVectorUsingGK(const double &val, const int eqnNum, LargeMatrix *gK)");
	if(!allocated) FatalError("Preliminary allocation must be performed before a value can be set.");
	RangeCheck(0,eqnNum,numberOfEquations-1);

	int i;
	double *F = pseudoLoadVector;

	if(restraint[eqnNum]) {
		if(additionalEquation[eqnNum]) 
		{cout<<"Equation# "<<eqnNum<<" is restrained with an additional equation"<<endl;
         cout<<" this function does not handle such a case "<<endl;
		 FatalError("this function cannot be used for this model because of additional equations.");
		 return;
		}
		for(i=0;i<numberOfEquations;i++) {
			if(!restraint[i]){
				if(gK->checkNonZeroLocation(i,eqnNum)){ //taking non-zero location information from unrestrained gK (passed thru parameters)
//				if(a->checkNonZeroLocation(i,eqnNum)){  //taking non-zero location information from gK in the equation class
//the solvers are coded in the assumption that they are symmetric and only lower triangle is stored
//it crashes when an element in the upper triangle is accessed. there are two ways to fix this:
//either the solver should accept the upper triangle (i,j) and return the corresponding value in the lower triangle i.e. (j,i)
//or make correct the parameters before it goes to the solver's function. 
//I am chosing the latter option. but one needs to be careful when using this function 
// on a non-symmetric matrix.
					if(i < eqnNum)
						F[i] -= val * (*gK)(eqnNum,i);
					else
						F[i] -= val * (*gK)(i,eqnNum);
				}
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
//
//  Z E R O I N G   M E T H O D S 
//
void Equations::zeroMatrix(void)
{
	WhoMethod("zeroMatrix(void)");
	if(!allocated) FatalError("Equations not allocated yet.");
	(*a).zeroYourself();
	int i=0,limit=numberOfEquations;
	double *F = pseudoLoadVector;
	for(i=0;i<limit;i++) F[i]=0;

}

void Equations::zeroLoadVector(void)
{
	WhoMethod("zeroLoadVector(void)");
	if(!allocated) FatalError("Equations not allocated yet.");
	int i=0,limit=numberOfEquations;
	double *F = loadVector;
	for(i=0;i<limit;i++) F[i]=0;

}

void Equations::zeroPseudoLoadVector(void)
{
	WhoMethod("zeroPseudoLoadVector(void)");
	if(!allocated) FatalError("Equations not allocated yet.");
	int i=0,limit=numberOfEquations;
	double *F = pseudoLoadVector;
	for(;i<limit;i++) F[i]=0;
}

void Equations::zeroResultantVector(void)
{
	WhoMethod("zeroResultantVector(void)");
	if(!allocated) FatalError("Equations not allocated yet.");
	int i=0,limit=numberOfEquations;
	double *F = resultantVector;
	for(;i<limit;i++) F[i]=0;
}

// xtang 121700
// xtang 07072003 Make it initializing with a prescribed val.
void Equations::zeroSolution(double val)
{
	WhoMethod("zeroResultantVector(void)");
	if(!allocated) FatalError("Equations not allocated yet.");
	int i=0,limit=numberOfEquations;
	double *F = solution;
	for(;i<limit;i++) F[i]=val;
}


/////////////////////////////////////////////////////////////////////////////
//
//  A D D / S E T   M E T H O D S
//
void Equations::addLinearEquations(LinearEquations &eqns, const int *map)
{ addMatrix(eqns.K,map,eqns.numberOfEquations,eqns.F); }

void Equations::addMatrix(const Matrix &matrix, const int *map, const int size, 
		const double *force, int threadNum)
{
	WhoMethod("addMatrix(const Matrix &, const int *, const int, const double *)");
	if(!allocated) FatalError("Equations not allocated yet.");

	//JV091104 - diverting to another function if the solver is MKLGBMAT or LAPACKGB
	if (sm==mklgbmat || sm==LapackGB) {
		addMatrixMKL(matrix, map, size, force);
		return;
	}

	//JV091104 - diverting to another function if the solver is PetSc
	if (sm==PetSc || sm==PetScUnsymm) {
		addMatrix2(matrix, map, size, force);
		return;
	}

	// This assembly process is fairly complicated. Be sure to make sure of what
// you are doing if you make any changes!
// 'matrix' should be in upper triangluar format.
// 'Kg' Matrix should behave like normal lower triangular matrix.
// Copying object variables to local variables.

	int i,j,k,m,I,J,ip,jp;
	double matrix_ij=0;
	double *F =0;
	bool *rest = restraint;

	F = pseudoLoadVector;

	LargeMatrix *Mat;
	Mat = a;
	LargeMatrix &Kg = *Mat;

// Assemble upper triangle of Ke and any force contribution
	for(i=0;i<size;i++) {
		I = map[i];
		RangeCheck(0,I,numberOfEquations-1);
		// Add in force if force vector exists and dof not restrained
		if(force!=0) {
			if(!rest[I]) F[I] += force[i];
			else if(additionalEquation[I]) {
				AdditionalEquation &aMi = *additionalEquation[I];
				for(k=0;k<aMi.numTerms;k++)
					F[aMi.master[k]] += force[i] * aMi.factor[k];
			}
		}
		for(j=0;j<size;j++) {
			if(i<j) matrix_ij = matrix(i,j);
			else matrix_ij = matrix(j,i);
			J=map[j];
			RangeCheck(0,J,numberOfEquations-1);
			// If J restrained move to lhs and all J restrained cases
			if(rest[J]) {
				if(I==J) Kg(I,I)=1;
				if(!rest[I]) {
					if(additionalEquation[J]) F[I] -= matrix_ij * additionalEquation[J]->constant;
					else F[I] -= matrix_ij * solution[J];
				} else if(additionalEquation[I]) {
					AdditionalEquation &aMi = *additionalEquation[I];
					for(k=0;k<aMi.numTerms;k++) {
						if(additionalEquation[J])
							F[aMi.master[k]] -= matrix_ij * aMi.factor[k] *
								additionalEquation[J]->constant;
						else
							F[aMi.master[k]] -= matrix_ij * aMi.factor[k] *
								solution[J];
					}
				}
				if(additionalEquation[J]) { 
					AdditionalEquation &aMj = *additionalEquation[J];
					if(rest[I]) {
						if (additionalEquation[I]) {
							AdditionalEquation &aMi = *additionalEquation[I];
							for(k=0;k<aMj.numTerms;k++) {
								ip = aMj.master[k];
								if(rest[ip]) Kg(ip,ip) = 1;
								else for(m=0;m<aMi.numTerms;m++) {
									jp = aMi.master[m];
									if(rest[jp]) {
//										cout << "master " <<jp<<" is restrained."<<endl;
										F[ip] -= matrix_ij * aMj.factor[k] *
												aMi.factor[m] * solution[jp];
									} else if(ip >= jp)
										Kg(ip,jp) += matrix_ij * aMj.factor[k] * 
												aMi.factor[m];
								}
							}
						}
					} else {
						for(k=0;k<aMj.numTerms;k++) {
							jp = aMj.master[k];
							if(rest[jp]) {
//								cout << "master " <<jp<<" is restrained."<<endl;
								F[I] -= matrix_ij * aMj.factor[k] * solution[jp];
							} else if(I >= jp) Kg(I,jp) += matrix_ij * aMj.factor[k];
						}
					}
				} 
			} else if(rest[I]) {  // If I is restrained and J is not
				if(additionalEquation[I]) {
					AdditionalEquation &aMi = *additionalEquation[I];
					for(k=0;k<aMi.numTerms;k++) {
						ip = aMi.master[k];
						if(rest[ip]) Kg(ip,ip) = 1;
						else if(aMi.master[k] >= J)
							Kg(aMi.master[k],J) += matrix_ij * aMi.factor[k];
					}
				}
			} else if(I>=J) // If neither restrained do normal assembly
				Kg(I,J) += matrix_ij;
		}
	}
#ifdef DEBUG_Level3
	Kg.printMatrix("Kg = ");
#endif
}
//JV101106 - added function to assemble matrix that does not use the operator() function
// this cannot be used with PETSC library
void Equations::addMatrix2(const Matrix &matrix, const int *map, const int size, 
		const double *force)
{
	WhoMethod("addMatrix2(const Matrix &, const int *, const int, const double *)");
	if(!allocated) FatalError("Equations not allocated yet.");

	//JV091104 - diverting to another function if the solver is MKLGBMAT or LAPACKGB
	if (sm==mklgbmat || sm==LapackGB) {
		addMatrixMKL(matrix, map, size, force);
		return;
	}
	/*
	//fix this function, if the matrix is converted to general, then the symmertic petsc solver does not work!
	//changes made in 3 places - convertSymmetricToGeneral(), assigning matrix_ij and in normal assembly Kg.addtoaij(I,J,val);
	if(matrix.isSymmetric)
		matrix.convertSymmetricToGeneral();
	*/

// This assembly process is fairly complicated. Be sure to make sure of what
// you are doing if you make any changes!
// 'matrix' should be in upper triangluar format.
// 'Kg' Matrix should behave like normal lower triangular matrix.
// Copying object variables to local variables.

	int i,j,k,m,I,J,ip,jp;
	double matrix_ij=0;
	double *F = pseudoLoadVector;
	bool *rest = restraint;
	LargeMatrix &Kg = *a;
	double	val;

// Assemble upper triangle of Ke and any force contribution
	for(i=0;i<size;i++) {
		I = map[i];
		RangeCheck(0,I,numberOfEquations-1);
		// Add in force if force vector exists and dof not restrained
		if(force!=0) {
			if(!rest[I]) F[I] += force[i];
			else if(additionalEquation[I]) {
				AdditionalEquation &aMi = *additionalEquation[I];
				for(k=0;k<aMi.numTerms;k++)
					F[aMi.master[k]] += force[i] * aMi.factor[k];
			}
		}
		for(j=0;j<size;j++) {
			if(i<j) matrix_ij = matrix(i,j);
			else matrix_ij = matrix(j,i);
			//since it is 'unsymmetric', we are not just assembling using upper triangle of Ke
			//matrix_ij = matrix(i,j);

			J=map[j];
			RangeCheck(0,J,numberOfEquations-1);
			// If J restrained move to lhs and all J restrained cases
			if(rest[J]) {
				//if(I==J) Kg(I,I)=1;
				if(I==J) Kg.setaij(I,I,1.0);
				if(!rest[I]) {
					if(additionalEquation[J]) F[I] -= matrix_ij * additionalEquation[J]->constant;
					else F[I] -= matrix_ij * solution[J];
				} else if(additionalEquation[I]) {
					AdditionalEquation &aMi = *additionalEquation[I];
					for(k=0;k<aMi.numTerms;k++) {
						if(additionalEquation[J])
							F[aMi.master[k]] -= matrix_ij * aMi.factor[k] *
								additionalEquation[J]->constant;
						else
							F[aMi.master[k]] -= matrix_ij * aMi.factor[k] *
								solution[J];
					}
				}
				if(additionalEquation[J]) { 
					AdditionalEquation &aMj = *additionalEquation[J];
					if(rest[I]) {
						if (additionalEquation[I]) {
							AdditionalEquation &aMi = *additionalEquation[I];
							for(k=0;k<aMj.numTerms;k++) {
								ip = aMj.master[k];
								//if(rest[ip]) Kg(ip,ip) = 1;
								if(rest[ip]) Kg.setaij(ip,ip, 1.0);
								else for(m=0;m<aMi.numTerms;m++) {
									jp = aMi.master[m];
									if(rest[jp]) {
//										cout << "master " <<jp<<" is restrained."<<endl;
										F[ip] -= matrix_ij * aMj.factor[k] *
												aMi.factor[m] * solution[jp];
									} else if(ip >= jp){
										val = matrix_ij * aMj.factor[k] * 
												aMi.factor[m];
										Kg.addtoaij(ip,jp,val);
									}
								}
							}
						}
					} else {
						for(k=0;k<aMj.numTerms;k++) {
							jp = aMj.master[k];
							if(rest[jp]) {
//								cout << "master " <<jp<<" is restrained."<<endl;
								F[I] -= matrix_ij * aMj.factor[k] * solution[jp];
							} else if(I >= jp){
								val = matrix_ij * aMj.factor[k];
								Kg.addtoaij(I,jp,val);
							}
						}
					}
				} 
			} else if(rest[I]) {  // If I is restrained and J is not
				if(additionalEquation[I]) {
					AdditionalEquation &aMi = *additionalEquation[I];
					for(k=0;k<aMi.numTerms;k++) {
						ip = aMi.master[k];
						if(rest[ip]) Kg.setaij(ip,ip,1.0);
						else if(aMi.master[k] >= J){
							val = matrix_ij * aMi.factor[k];
							Kg.addtoaij(aMi.master[k],J,val);
						}
					}
				}
			} else if(I>=J){ // If neither restrained do normal assembly
				//since it is 'unsymmetric', we are not just assembling using upper triangle of Ke
				val = matrix_ij;
				Kg.addtoaij(I,J,val);
			}
		}
	}
#ifdef DEBUG_Level3
	Kg.printMatrix("Kg = ");
#endif
}



//JV 091104 - function to add matrix to general (unsymmetric) Kg (MKL)
void Equations::addMatrixMKL(const Matrix &matrix, const int *map, const int size, 
		const double *force)
{
	WhoMethod("addMatrixMKL(const Matrix &, const int *, const int, const double *)");
	if(!allocated) FatalError("Equations not allocated yet.");

// This assembly process is fairly complicated. Be sure to make sure of what
// you are doing if you make any changes!
// Copying object variables to local variables.

	int i,j,I,J;
	double matrix_ij=0;
	double *F = pseudoLoadVector;
	bool *rest = restraint;
	LargeMatrix &Kg = *a;

	double big=10e30;

/*
// Assemble Ke and any force contribution
	for(i=0;i<size;i++) {
		I = map[i];
		RangeCheck(0,I,numberOfEquations-1);
		// Add in force if force vector exists and dof not restrained
		if(force!=0) {
			if(!rest[I]) F[I] += force[i];
		}
		for(j=0;j<size;j++) {
			matrix_ij = matrix(i,j);
			J=map[j];
			RangeCheck(0,J,numberOfEquations-1);
			if(rest[J] && I==J) {
				Kg(I,I)=big;
				solution[J] *= big;
			} else // If neither restrained do normal assembly
				Kg(I,J) += matrix_ij;
		}
	}
*/

//had to change the algo to the one above cos decomposition method gave a singular method.
//so, not we are going for penalty method

// Assemble Ke and any force contribution
	for(i=0;i<size;i++) {
		I = map[i];
		RangeCheck(0,I,numberOfEquations-1);
		// Add in force if force vector exists and dof not restrained
		if(force!=0) {
			if(!rest[I]) F[I] += force[i];
		}
		for(j=0;j<size;j++) {
			matrix_ij = matrix(i,j);
			J=map[j];
			RangeCheck(0,J,numberOfEquations-1);
			// If J restrained move to lhs and all J restrained cases
			if(rest[J]) {
				if(I==J) Kg(I,I)=1;
				if(!rest[I]) {
					F[I] -= matrix_ij * solution[J];
				} 
			} else if(rest[I]) {  // If I is restrained and J is not
				continue;
			} else // If neither restrained do normal assembly
				Kg(I,J) += matrix_ij;
		}
	}

#ifdef DEBUG_Level3
	Kg.printMatrix("Kg = ");
#endif
}


void Equations::subtractLinearEquations(LinearEquations &eqns, const int *map)
{
	WhoMethod("subtractLinearEquations(LinearEquations &, const int)");
	if(!allocated) FatalError("Equations not allocated yet.");
	int i,j;
	for(i=0;i<eqns.numberOfEquations;i++) {
		eqns.F[i] = -eqns.F[i];
		for(j=0;j<eqns.numberOfEquations;j++) 
			eqns.K(i,j) = - eqns.K(i,j);
	}
	addMatrix(eqns.K,map,eqns.numberOfEquations,eqns.F, 0); //any additions or subtractions assigned to a
	for(i=0;i<eqns.numberOfEquations;i++) {
		eqns.F[i] = -eqns.F[i];
		for(j=0;j<eqns.numberOfEquations;j++) 
			eqns.K(i,j) = - eqns.K(i,j);
	}
}

void Equations::assembleResultantVector(const int *map, const int size, 
	const double *force)
{
	int i,I,k;
	double *F = resultantVector;
	WhoMethod("assembleResultantVector(const int *, const int, const double *)");
	if(!allocated) FatalError("Equations not allocated yet.");
	for(i=0;i<size;i++) {
		I = map[i];
		RangeCheck(0,I,numberOfEquations-1);
		if(additionalEquation[I]) {
			AdditionalEquation &aMi = *additionalEquation[I];
			for(k=0;k<aMi.numTerms;k++)
				if(!restraint[aMi.master[k]])
					F[aMi.master[k]] += force[i] * aMi.factor[k];
		} else if(!restraint[I]) F[I] += force[i];
	}
}
void Equations::addToPseudoLoadVector(const int *map, const int size, 
	const double *force)
{
	int i,I,k;
	double *F = pseudoLoadVector;
	WhoMethod("addToPseudoVector(const int *, const int, const double *)");
	if(!allocated) FatalError("Equations not allocated yet.");
	for(i=0;i<size;i++) {
		I = map[i];
		RangeCheck(0,I,numberOfEquations-1);
		if(additionalEquation[I]) {
			AdditionalEquation &aMi = *additionalEquation[I];
			for(k=0;k<aMi.numTerms;k++)
				if(!restraint[aMi.master[k]])
					F[aMi.master[k]] += force[i] * aMi.factor[k];
		} else if(!restraint[I]) F[I] += force[i];
	}
}

void Equations::addToLoadVector(const int *map, const int size, 
	const double *force, int loadVecNum)
{
	if(force==0) return;
	int i,I,k;

	double *F = NULL;
	F = loadVector;


	WhoMethod("addToLoadVector(const int *, const int, const double *)");
	if(!allocated) FatalError("Equations not allocated yet.");
	for(i=0;i<size;i++) {
		I = map[i];
		RangeCheck(0,I,numberOfEquations-1);
		if(additionalEquation[I]) {
			AdditionalEquation &aMi = *additionalEquation[I];
			for(k=0;k<aMi.numTerms;k++)
				if(!restraint[aMi.master[k]])
					F[aMi.master[k]] += force[i] * aMi.factor[k];
		} else if(!restraint[I]) F[I] += force[i];
	}
}

void Equations::addLoadToEqn(const double &load, const int eqnNum, 
		bool ignoreIfRestrained)
{
	WhoMethod("addLoadToEqn(const double &,const int)");
	double *F = loadVector;
	if(!allocated) FatalError("Equations not allocated yet.");
	RangeCheck(0,eqnNum,numberOfEquations-1);
	if(restraint[eqnNum]) {
		if(!additionalEquation[eqnNum]) {
			if(ignoreIfRestrained) return;
			FatalError("Cannot add load to a restrained degree of freedom.");
		}
		AdditionalEquation &aM=*additionalEquation[eqnNum];
		for(int i=0;i<aM.numTerms;i++) 
			if(!restraint[aM.master[i]])
				F[aM.master[i]] += load * aM.factor[i];
	} else F[eqnNum] += load;
}

void Equations::setLoadForEqn(const double &load, const int eqnNum)
{
	WhoMethod("setLoadForEqn(const double &,const int)");
	if(!allocated) FatalError("Equations not allocated yet.");
	RangeCheck(0,eqnNum,numberOfEquations-1);
	loadVector[eqnNum] = load;
} 

/////////////////////////////////////////////////////////////////////////////
//
//  P R I N T I N G   M E T H O D S
//
void Equations::printResultantVector(FormWriter &fr)
{
	int i;
	WhoMethod("printResultantVector(FormWriter &)");
	fr.h2("Resultant Load Vector");
	if(!allocated) {fr << "\tEquations not allocated.\n"; return;}
	fr << "numberOfEquations = " << numberOfEquations << '\n';
	fr << "Eqn.\tFr\tRestraint\n";
	for(i=0;i<numberOfEquations;i++)
		fr << i << '\t' << resultantVector[i] << '\t' << restraint[i] << endl;
}
void Equations::printEquivalentLoadVector(FormWriter &fr)
{
	int i;
	WhoMethod("printLoadVector(FormWriter &)");
	fr.h2("Equivalent Load Vector");
	if(!allocated) {fr << "\tEquations not allocated.\n"; return;}
	fr << "numberOfEquations = " << numberOfEquations << '\n';
	fr << "Eqn.\tFeq\tRestraint\n";
	for(i=0;i<numberOfEquations;i++)
		fr << i << '\t' << loadVector[i] << '\t' << restraint[i] << endl;
}
void Equations::printLoadVector(FormWriter &fr)
{
	int i;
	WhoMethod("printLoadVector(FormWriter &)");
	fr.h2("Load Vector");
	if(!allocated) {fr << "\tEquations not allocated.\n"; return;}
	fr << "numberOfEquations = " << numberOfEquations << '\n';
	fr << "Eqn.\tF\tRestraint\n";
	for(i=0;i<numberOfEquations;i++)
		fr << i << '\t' << loadVector[i] << '\t' << restraint[i] << endl;
}

void Equations::printMatrix(FormWriter &fr)
{
	WhoMethod("printMatrix(FormWriter &)");
	fr.h2("Coefficient Matrix");
	if(!allocated) {fr << "\tNot allocated.\n"; return;}
//	(*a).printMatrix("Kg = ");
	(*a).printMatrix("Kg = ",fr.getOutputStream());
}

void Equations::printSolution(FormWriter &fr)
{
	WhoMethod("printSolution(FormWriter &)");
	fr.h2("Solution");
	if(!allocated) {fr << "\tNot allocated.\n"; return;}
	if(!solved) {fr << "\tNot solved.\n"; return;}
	for(int i=0;i<numberOfEquations;i++)
		fr << i << '\t' << solution[i] << '\n';
	fr << "End Solution." << endl;
}

void Equations::printAdditionalEquations(FormWriter &fr)
{
	WhoMethod("printAdditionalEquations(FormWriter &)");
	fr.h2("AdditionalEquations");
	if(aeList.begin()) {
		do {
			fr << *aeList() << '\n';
		} while(aeList++);
	} else {
		fr << "There are no additional equations." << endl;
	}
}

void Equations::printAdditionalEquations(ostream *os)
{
	WhoMethod("printAdditionalEquations(ostream *os)");
	*os << "AdditionalEquations" << endl;
	if(aeList.begin()) {
		do {
			*os << *aeList() << '\n';
		} while(aeList++);
	} else {
		*os << "There are no additional equations." << endl;
	}
}

void Equations::printRestraints(FormWriter &fr)
{
	WhoMethod("printRestraints(FormWriter &)");
	fr.h2("Restraints");
	if(restraintList.begin()) {
		int i=0;
		fr << "The following unknowns are restrained.";
		do {
			if((i%=8)==0) fr << "\n\t"; else fr << " : ";
			fr << *restraintList();
			i++;
		} while(restraintList++);
		fr << '\n';
	} else {
		fr << "No unknowns are restrained." << endl;
	}
}

/////////////////////////////////////////////////////////////////////////////
//
//  M E M O R Y   A L L O C A T I O N 
//
#ifdef SPARSE
#include "sparse.hpp"
#endif
#ifdef OLAF   // VSS -OLAF
#include "oSparse.hpp"
#endif
#ifdef PROFILE
#include "promat.hpp"
#endif
/* this is an MKL function, so group under the same ifdef statement
#ifdef MKLGBMAT
#include "mklgbmat.hpp"
#endif
*/
#ifdef LAPACKGB
#include "Lapackgb.hpp"
#endif
#ifdef PETSC
#include "petscsparse.hpp"
#include "petscsparseunsymm.hpp"
#endif
#ifdef MKL
#include "mklgbmat.hpp"
#include "MKLPardiso.hpp"
#include "MKLPardisoUnsymm.hpp"
#include "MKLCG.hpp"
#endif
#ifdef PARCG
#include "ParCG.hpp"
#endif
#ifdef PARDISO
#include "Pardiso.hpp"
#endif
#ifdef WSMP
#include "WSMPSymm.hpp"
#endif

void Equations::CreateMatrix(StorageMethod storageMethod, int n)
{
	sm = storageMethod;
	switch (sm) {
#ifdef SPARSE
		case Sparse:
			a = new SymmetricSparseMatrix(n);
			((SymmetricSparseMatrix*)a)->setMAX_NUMBER_ITERATIONS(MaxIterations);
			((SymmetricSparseMatrix*)a)->setSOLUTION_EPSILON(sol_epsilon);
			break;
#endif
#ifdef PROFILE
		case Profile:
			a = new ProfileMatrix(n);
			break;
#endif
#ifdef OLAF
		case Osparse:
			a = new SymmetricOSparseMatrix(n);
			break;
#endif
#ifdef LAPACKGB
		case LapackGB:
			a = new LapackGBMatrix(n);
			break;
#endif
#ifdef PETSC
		case PetSc:
			a = new PETScSparseMatrix(n);
			break;
		case PetScUnsymm:
			a = new PETScUnsymmSparseMatrix(n);
			break;
#endif
#ifdef MKL
		case mklgbmat:
			a = new MKLGeneralMatrix(n);
			break;
		case MKLPardiso:
			a = new MKLPardisoSymmMatrix(n);
			break;
		case MKLPardisoUnsymm:
			a = new MKLPardisoUnsymmMatrix(n);
			break;
		case MKLCG:
			a = new MKLCGSymmMatrix(n);
			((MKLCGSymmMatrix*)a)->setMAX_NUMBER_ITERATIONS(MaxIterations);
			((MKLCGSymmMatrix*)a)->setSOLUTION_EPSILON(sol_epsilon);
			break;

#endif
#ifdef PARCG
		case ParCG:
			a = new ParCGMatrix(n);
			((ParCGMatrix*)a)->setMAX_NUMBER_ITERATIONS(MaxIterations);
			((ParCGMatrix*)a)->setSOLUTION_EPSILON(sol_epsilon);
			break;
#endif
#ifdef PARDISO
		case Pardiso:
			a = new PardisoSymmMatrix(n);
			break;
#endif
#ifdef WSMP
		case WSMPSymm:
			a = new WSMPSymmMatrix(n);
			break;
#endif
		default:
			cerr << "StorageMethod = " << sm << endl;
			FatalError("Unrecognized storage method.");
			break;
	}

	if (a != NULL ) {
		a->setReplaceZeroDiagonal_Option(zeroDiagonalReplacementValue);
		a->setVerboseFlag(verbose);
		a->setUsePreviousFillInOrdering(UsePreviousFillInOrdering);
		a->setUsePreviousFactorizedMatrix(UsePreviousFactorizedMatrix);
		a->setReorderingScheme(reorderingScheme);
		a->setNum_threads_for_solver(num_threads_for_solver);
	}
}
void Equations::CreateAndInitializeArrays(int n)
{
	int i;
	additionalEquation = new AdditionalEquation *[n];
	for(i=0;i<n;i++) additionalEquation[i] = 0;
	restraint = new bool [n];
	for(i=0;i<n;i++) restraint[i] = false;
	solution = new double [n];
	for(i=0;i<n;i++) solution[i] = 0;
	lastSolution = new double [n];
	for(i=0;i<n;i++) lastSolution[i] = 0;
	pseudoLoadVector = new double [n];
	for(i=0;i<n;i++) pseudoLoadVector[i] = 0;
	loadVector = new double [n];
	for(i=0;i<n;i++) loadVector[i] = 0;
	resultantVector = new double [n];
	for(i=0;i<n;i++) resultantVector[i] = 0;
	residualVector = new double [n];
	for(i=0;i<n;i++) residualVector[i] = 0;

}
void Equations::InitializeArrays(int n)
{
	int i;
	for(i=0;i<n;i++) additionalEquation[i] = 0;
	for(i=0;i<n;i++) restraint[i] = false;
	for(i=0;i<n;i++) solution[i] = 0;
	for(i=0;i<n;i++) lastSolution[i] = 0;
	for(i=0;i<n;i++) pseudoLoadVector[i] = 0;
	for(i=0;i<n;i++) loadVector[i] = 0;
	for(i=0;i<n;i++) resultantVector[i] = 0;
	for(i=0;i<n;i++) residualVector[i] = 0;

}
void Equations::cleanRestraints()
{
	int i;
	for(i=0;i<numberOfEquations;i++) restraint[i] = false;

	if(restraintList.AreYouEmpty())
		restraintList.DestroyList();

	if(restraintList.AreYouEmpty() >0 ) {
		cout << "Could not empty the restraintList !" << endl;
		exit(1);
	}

}

void Equations::setupStorage(StorageMethod storageMethod, int n)
{
	WhoMethod("setStorage(StorageMethod, int)");
	if(allocated) releaseMemory();
	numberOfEquations = n;
	CreateAndInitializeArrays(n);
	CreateMatrix(storageMethod, n);

	if(restraintList.begin()) {
		do {
			restraint[*restraintList()] = true;
		} while (restraintList++);
	}
	emptyRestraintList(); // since this list is never going to be used again during this analysis, its better to release the memory taken up by this list.

	if(aeList.begin()) {
		do {
			if(additionalEquation[aeList()->slave] == 0) {
				additionalEquation[aeList()->slave] = aeList();
				restraint[aeList()->slave] = true;
			} else {
				//cerr << "slave dof: " << aeList()->slave << endl; //JV082707 : added more debugging info
				if( aeList()->slave >= n ){
					cerr << "You have additional equations slaving a dof that larger than the number of equations." << endl;
					cerr << "Check your additional equations specifications." << endl;
					FatalError("Panic! Conflicting equations should not exist by this point!");
				}
				if(*aeList()!=*additionalEquation[aeList()->slave]) {
					cerr << *aeList() << endl; 
					cerr << *additionalEquation[aeList()->slave] << endl;
					cerr << "slave dof: " << aeList()->slave << endl; //JV082707 : added more debugging info
					FatalError("Panic! Conflicting equations should not exist by this point!");
				}
			}
		} while (aeList++);
	}
	allocated = true;
}

void Equations::initializeStorage(int n)
{
	WhoMethod("initializeStorage(int)");
	if(numberOfEquations == n)
		InitializeArrays(n);
	else{
		cout << "you are trying to change the number of equations!!!" << endl;
		exit(1);
	}

	deleteLargeMatrix();
	CreateMatrix(sm, n); //this function creates only one matrix
						 //this does not work if multiple matrices are used for parallel assembly

	allocated = true;
}

LargeMatrix* Equations::createCopyOfMatrix()
{
	LargeMatrix *m;
	m=a->clone();
	return(m);
}

void Equations::specifyNonZeroLocation(const int *map, const int size)
//assumes the symmetric global stiffness matrix is storing the lower diagonal info.
//if the global stiffness matrix is actually storing the upper diagonal info,
//then the 'large matrix' object should handle that specifically.
{
	int I,J,i,j,k,klimit,m,mlimit,ip,jp;
	LargeMatrix &Kg = *a;
	bool *rest = restraint;
	for(i=0;i<size;i++) {
		I=map[i];
		for(j=0;j<size;j++) {
			J=map[j];
			if(rest[J]) {
				if(I==J) Kg.specifyNonZeroLocation(I,I);
				if(additionalEquation[J]) {
					AdditionalEquation &aMj = *additionalEquation[J];
					if(rest[I]) {
						if(additionalEquation[I]) {
							AdditionalEquation &aMi = *additionalEquation[I];
							for(k=0,klimit=aMj.numTerms,mlimit=aMi.numTerms;k<klimit;k++) {
								ip = aMj.master[k];
								if(rest[ip]) 
									Kg.specifyNonZeroLocation(ip,ip);
								else for(m=0;m<mlimit;m++) {
									jp = aMi.master[m];
									if(!rest[jp] && ip >= jp) 
										Kg.specifyNonZeroLocation(ip,jp);
								}
							}
						}
					} else
						for(k=0;k<aMj.numTerms;k++)
							if(!rest[aMj.master[k]] && I >= aMj.master[k]) 
								Kg.specifyNonZeroLocation(I,aMj.master[k]);
				}
			} else if(rest[I]) {  // If I is restrained and J is not
				if(additionalEquation[I]) {
					AdditionalEquation &aMi = *additionalEquation[I];
					for(k=0,klimit=aMi.numTerms;k<klimit;k++) {
						ip = aMi.master[k];
						if(rest[ip]) Kg.specifyNonZeroLocation(ip,ip);
						else if(ip >= J) Kg.specifyNonZeroLocation(ip,J);
					}
				}
			} else if(I>=J) { // If neither restrained do normal assembly
				Kg.specifyNonZeroLocation(I,J);
				//assumes, symm. matrix is lower diagonal. large matrix object should handle it specifically if it is upper diagonal.
				//if the matrix is general, the matrix should handle specifying (J,I) as non-zero location
			}
		}
	}
}

void Equations::allocate(void)
{
	int i,L=numberOfEquations;
	double *u = solution;
	// Copy the additional Equation constants into solution. 
	// After this no more addtional equations may be added....
	for(i=0;i<L;i++) 
		if(additionalEquation[i])
			u[i]=additionalEquation[i]->constant;
	(*a).allocate();

}

/////////////////////////////////////////////////////////////////////////////
//
//  S O L V I N G   M E T H O D S
//

#include <iomanip>
double *Equations::solve(double *initGuess)
{
	//JV091104 - diverting to another function if the solver is MKLGBMAT
	if (sm==mklgbmat || sm==LapackGB) {
		return solveMKL(initGuess);
	}

	WhoMethod("solve(double *)");
	int L=numberOfEquations;
	if(!allocated) FatalError("Equations not allocated yet.");
	bool cannotContinue = false;

//whitcomb
	if(verbose >= Basic){
		displayTime("", *outStream);
		cout<<"Solve"<<endl;
		cout<<"Number of equations : "<< L << endl; 
	}

	ofstream SolverTestFile;

	if(verbose >= Min){
		SolverTestFile.open("PRESOLVEpseudoLoadVector.txt");
		SolverTestFile.setf(ios::scientific, ios::floatfield);
		for(int i=0;i<L;i++) SolverTestFile << DOUBLE_FORMAT << pseudoLoadVector[i] << endl;
		SolverTestFile.close();
		SolverTestFile.open("PRESOLVEloadVector.txt");
		SolverTestFile.setf(ios::scientific, ios::floatfield);
		for(int i=0;i<L;i++) SolverTestFile << DOUBLE_FORMAT << loadVector[i] << endl;
		SolverTestFile.close();
	}


	int i;
	// Create the real loadVector in solution vector except for restrained values, 
	//	because the solution vector already has the solution for restrained dof's.
	for(i=0;i<L;i++) 
		if(!restraint[i]) solution[i] = loadVector[i] + pseudoLoadVector[i];

	// Copy the additional Equation constants in and check to see if master dofs
	// are restrained.  If master restrained, complain and leave.
	// This is final check before solving....
	for(i=0;i<L;i++) {
		if(additionalEquation[i]) {
			solution[i]=additionalEquation[i]->constant;
#ifdef OBSOLETE
			for(j=0;j<additionalEquation[i]->numTerms;j++)
				if(restraint[additionalEquation[i]->master[j]]) {
					cerr << "Additional Equation Conflict: " << (*additionalEquation[i]) << endl;
					Warning("A master dof is restrained.  Program not setup for this yet.");
					cannotContinue = true;
				}
#endif
		}
	}
	if(cannotContinue) FatalError("Cannot continue, due to previous warnings.");

//	displayTime("Time to prepare load vector", *outStream);

	if(verbose==Max) {
		SolverTestFile.open("GlobalStiffness.txt");
		(*a).printMatrixAddressFormat("Global Stiffness",&SolverTestFile);
		SolverTestFile.close();
		//SolverTestFile.open("GlobalStiffness_CSR.txt");
		//(*a).printCompressedRowFormat("Global Stiffness",&SolverTestFile);
		//SolverTestFile.close();

		(*a).printCSRFormatASCII("CSR_ASCII.txt", cout);
		(*a).printCSRFormatBinary("CSR_Binary.dat", cout);
	}
	if(verbose >= Min){
		SolverTestFile.open("LoadVector.txt");
		SolverTestFile.setf(ios::scientific, ios::floatfield);
		for(i=0;i<L;i++) SolverTestFile << DOUBLE_FORMAT << solution[i] << endl;
		SolverTestFile.close();
	}


	// Not all solvers can use lastSolution, but provide it anyway just in case.
	if(verbose >= Basic){
		cout << "Solving:"; cout.flush();
	}
	double *tmpGuess = new double [L];
	if(initGuess) {
		for(i=0;i<L;i++) {
			if(restraint[i] || additionalEquation[i]) tmpGuess[i] = solution[i];
			else tmpGuess[i] = initGuess[i];
//			cout << i << '\t' << initGuess[i] << '\t' << tmpGuess[i] << '\n';
		}
		(*a).solve(solution,tmpGuess);
	} else {
		(*a).solve(solution,lastSolution);
	}
	if(verbose >= Basic){
		displayTime("Time to finish solve step", *outStream);
		cout << "\nDone Solving:" << endl;
	}

	if(verbose >= Min){
		SolverTestFile.open("Solution.txt");
		SolverTestFile.setf(ios::scientific, ios::floatfield);
		SolverTestFile << L << endl;
		for(i=0;i<L;i++) SolverTestFile << DOUBLE_FORMAT << solution[i] << endl;
		SolverTestFile.close();
	}


	// Apply the additional equations and save the solution in lastSolution
	for(i=0;i<L;i++) {
		if(fabs(tmpGuess[i]-solution[i])>1e-9)
//		cout << i << '\t' << tmpGuess[i] << '\t' << solution[i] << '\n';
		lastSolution[i] = solution[i];
		if(additionalEquation[i]) 
			solution[i] = additionalEquation[i]->resolve(solution);
	}
	delete [] tmpGuess;
	return solution;
}

double *Equations::solveMKL(double *initGuess)
//using penalty method
{
	WhoMethod("solveMKL(double *)");
	int L=numberOfEquations;
	if(!allocated) FatalError("Equations not allocated yet.");
	bool cannotContinue = false;

//whitcomb
cout<<"Solve"<<endl;

	int i;
	// Create the real loadVector in solution vector except for restrained values, 
	//	because the solution vector already has the solution for restrained dof's.
	for(i=0;i<L;i++) 
		if(!restraint[i]) solution[i] = loadVector[i] + pseudoLoadVector[i];
//		if(!restraint[i]) solution[i] = loadVector[i];

	if(verbose==Max) (*a).printMatrix("Global Stiffness");

	// Copy the additional Equation constants in and check to see if master dofs
	// are restrained.  If master restrained, complain and leave.
	// This is final check before solving....
	// <-- snipped code -->
	if(cannotContinue) FatalError("Cannot continue, due to previous warnings.");

	// Not all solvers can use lastSolution, but provide it anyway just in case.
	cout << "Solving:"; cout.flush();
	double *tmpGuess = new double [L];
	if(initGuess) {
		for(i=0;i<L;i++) {
			if(restraint[i] || additionalEquation[i]) tmpGuess[i] = solution[i];
			else tmpGuess[i] = initGuess[i];
//			cout << i << '\t' << initGuess[i] << '\t' << tmpGuess[i] << '\n';
		}
		(*a).solve(solution,tmpGuess);
	} else (*a).solve(solution,lastSolution);
	cout << "\nDone Solving:" << endl;
	// Apply the additional equations and save the solution in lastSolution
	for(i=0;i<L;i++) {
		if(fabs(tmpGuess[i]-solution[i])>1e-9)
//		cout << i << '\t' << tmpGuess[i] << '\t' << solution[i] << '\n';
		lastSolution[i] = solution[i];
		if(additionalEquation[i]) 
			solution[i] = additionalEquation[i]->resolve(solution);
	}
	delete [] tmpGuess;
	return solution;
}
//===================================================================================================
//==== JV091808 : this function (solveUsingResidual) does not seem to be used anywhere!!!============
//===================================================================================================
double *Equations::solveUsingResidual(double *initGuess)
{
	WhoMethod("solveUsingResidual(double *)");
	int L=numberOfEquations;
	if(!allocated) FatalError("Equations not allocated yet.");
	bool cannotContinue = false;
	double *deltaSolution = residualVector;

	int i;
	// Copy difference in the applied loads into solution vector except for restrained values, because
	// the solution vector already has the solution for restrained dof's.
	for(i=0;i<L;i++) 
		if(!restraint[i]) deltaSolution[i] = loadVector[i]-resultantVector[i];
		else deltaSolution[i] = 0;

	if(verbose==Max) (*a).printMatrix("Global Stiffness");

	// Copy the additional Equation constants in and check to see if master dofs
	// are restrained.  If master restrained, complain and leave.
	// This is final check before solving....
	for(i=0;i<L;i++) {
		if(additionalEquation[i]) {
			deltaSolution[i]=0;
#ifdef OBSOLETE
			for(j=0;j<additionalEquation[i]->numTerms;j++)
				if(restraint[additionalEquation[i]->master[j]]) {
					cerr << "Additional Equation Conflict: " << (*additionalEquation[i]) << endl;
					Warning("A master dof is restrained.  Program not setup for this yet.");
					cannotContinue = true;
				}
#endif
		}
	}
	if(cannotContinue) FatalError("Cannot continue, due to previous warnings.");

	// Not all solvers can use lastSolution, but provide it anyway just in case.
	cout << "Solving:"; cout.flush();
	double *tmpGuess = new double [L];
	if(initGuess) {
		for(i=0;i<L;i++) 
			if(restraint[i] || additionalEquation[i]) tmpGuess[i] = solution[i];
			else tmpGuess[i] = initGuess[i];
		(*a).solve(deltaSolution,tmpGuess);
	} else (*a).solve(deltaSolution,lastSolution);
	cout << "\nDone Solving:" << endl;
	// Save deltaSolution in lastSolution
	for(i=0;i<L;i++) lastSolution[i] = deltaSolution[i];
	// Update the actual solution (restrained values are zero at this point)
	for(i=0;i<L;i++) if(!restraint[i]) solution[i] += deltaSolution[i];
	// Apply the additional equations 
	for(i=0;i<L;i++)
		if(additionalEquation[i])
			solution[i] = additionalEquation[i]->resolve(solution);
	delete [] tmpGuess;
	return solution;
}

double Equations::calculateResidualNorm(void)
{
	WhoMethod("calculateResidualNorm(void)");
	int i,L=numberOfEquations;
	double L2Norm=0,tmp;
	double *r = residualVector;
	for(i=0;i<L;i++) {
		if(!restraint[i]) {
			tmp = r[i] = loadVector[i]-resultantVector[i];
			L2Norm += tmp*tmp;
		}
	}
	return sqrt(L2Norm);
}

void Equations::vectorProduct(double *result,double *vector)
{
	(*a).vectorProduct(result,vector);
}


/////////////////////////////////////////////////////////////////////////////
//
//  C O N S T R U C T O R S / D E S T R U C T O R S
//
Equations::Equations() 
{ 
	WhoClass("Equations"); WhoMethod("Equations(...)");
	a = 0;
	zeroDiagonalReplacementValue=0;//+ DG_Feb2004
	numberOfEquations = 0;
	additionalEquation = 0;
	restraint = 0;
	solution = 0;
	lastSolution = 0;
	loadVector = 0;
	pseudoLoadVector = 0;
	resultantVector = 0;
	residualVector = 0;
	allocated = false;
	solved=false;
	sm = Sparse;  
	verbose = Off;
	outStream = &cout;
	UsePreviousFillInOrdering=false;
	UsePreviousFactorizedMatrix=false;
	num_threads_for_solver=omp_get_max_threads();

	
	//iterative solver properties
	MaxIterations	=10000;
	sol_epsilon		=1e-4;
	//direct solver properties
	reorderingScheme=METIS;

	numberOfOpenMPThreads=1;
	//a2=0;
}

Equations::~Equations()
{ 
	WhoMethod("~Equations()");
	releaseMemory();
	numberOfEquations=0;
}
void Equations::releaseMemory()
{
	if(allocated) { // ie. is storage already allocated
		delete [] additionalEquation;
		delete [] restraint;
		delete [] solution;
		delete [] lastSolution;
		delete [] pseudoLoadVector;
		delete [] loadVector;
		delete [] resultantVector;
		delete [] residualVector;
		deleteLargeMatrix();

	}

	allocated=false;
}
void Equations::deleteLargeMatrix()
{
	if(a) {
		delete a;a=0;
	}
}

void Equations::setVerboseFlag(VerboseLevel level)
{
	verbose=level;
}
void Equations::setMaxIterations(int max_iter)
{
	MaxIterations=max_iter;
}
void Equations::setSOLUTION_EPSILON(double epsilon)
{
	sol_epsilon=epsilon;
}
void Equations::setReorderingScheme(string scheme)
{
	ChangeToUpper(scheme)
	if(scheme=="METIS"){
		reorderingScheme=METIS;
	}else if(scheme=="MMD"){
		reorderingScheme=MMD;
	}else if(scheme=="PAR"){
		reorderingScheme=PAR;
	}else{
		*outStream << "Not a valid reordering scheme" << endl;
		exit(1);
	}
	*outStream << "Setting reordering scheme to " << scheme << endl;

}


/////////////////////////////////////////////////////////////
// Whitcomb

double *Equations::incrementalSolutionUsingResidual(double *initGuess)
{
	WhoMethod("incrementalSolutionUsingResidual(double *)");
	int L=numberOfEquations;
	if(!allocated) FatalError("Equations not allocated yet.");
	double *deltaSolution = residualVector;

	if(verbose >= Basic){
		displayTime("", *outStream);
		cerr << "Solving:"; cerr.flush();
	}

	int i;
	// Copy difference in the applied loads into solution vector except for restrained values, because
	// the solution vector already has the solution for restrained dof's.
	for(i=0;i<L;i++) 
		if(!restraint[i]) deltaSolution[i] = loadVector[i]-resultantVector[i];
		else deltaSolution[i] = 0;

	//JV092208: why do we have this for loop??? if there is an additional eq,  the corresponding slave is already restrained, so this is alredy handled in the previous for loop.
	for(i=0;i<L;i++) {
		if(additionalEquation[i]) {
			deltaSolution[i]=0;
		}
	}

	ofstream SolverTestFile;
	if(verbose==Max) {
//		SolverTestFile.open("GlobalStiffnessNL.txt");
		SolverTestFile.open("GlobalStiffness.txt");
		(*a).printMatrixAddressFormat("Global Stiffness",&SolverTestFile);
		SolverTestFile.close();
	}
	if(verbose>=Min){
//		SolverTestFile.open("deltaSolution.txt");
		SolverTestFile.open("LoadVector.txt");
		SolverTestFile.setf(ios::scientific, ios::floatfield);
		for(i=0;i<L;i++) SolverTestFile << DOUBLE_FORMAT << deltaSolution[i] << endl;
		SolverTestFile.close();
	}


	// Not all solvers can use lastSolution, but provide it anyway just in case.
	double *tmpGuess = new double [L];
	if(initGuess) {
		for(i=0;i<L;i++) 
			if(restraint[i] || additionalEquation[i]) tmpGuess[i] = solution[i];
			else tmpGuess[i] = initGuess[i];
		(*a).solve(deltaSolution,tmpGuess);
	} else (*a).solve(deltaSolution,lastSolution);

	if(verbose >= Min){
		SolverTestFile.open("Solution.txt");
		SolverTestFile.setf(ios::scientific, ios::floatfield);
		SolverTestFile << L << endl;
		for(i=0;i<L;i++) SolverTestFile << DOUBLE_FORMAT << deltaSolution[i] << endl;
		SolverTestFile.close();
	}
	if(verbose >= Basic){
		displayTime("Time to finish solve step", *outStream);
		cout << "\nDone Solving:" << endl;
	}
	
	delete [] tmpGuess;
	return deltaSolution;
}  // end of incrementalSolutionUsingResidual

////////////////////////////////

void Equations::addIncrementalDisplacements(double *deltaSolution)
{int L=numberOfEquations;
 int i;

	// Save deltaSolution in lastSolution
	for(i=0;i<L;i++) lastSolution[i] = deltaSolution[i];
	// Update the actual solution (restrained values are zero at this point)
	for(i=0;i<L;i++) if(!restraint[i]) solution[i] += deltaSolution[i];
	// Apply the additional equations 
	for(i=0;i<L;i++)
		if(additionalEquation[i])
			solution[i] = additionalEquation[i]->resolve(solution);

}//end of addIncrementalDisplacements


//////////////////////////// xtang 990822

void Equations::applyAdditionalEquations()
{
	int L=numberOfEquations;
    int i;
    // Apply the additional equations 
	for(i=0;i<L;i++)
		if(additionalEquation[i])
			solution[i] = additionalEquation[i]->resolve(solution);

}//end of applyAdditionalEquations

//--------------------------------------------------------------------
#undef BETA_OUT
#define BETA_OUT (*outStream)
bool Equations::readCommands(ifstream *defaultStream)
{
	char *localTokenList[20];
	int  numberOfTokens;
	int foundMatch;

	ifstream *currentStream;

    BETA_OUT<<"Current stream in 'Equations::readCommands'="<<(*defaultStream)<<endl;
	currentStream = defaultStream;

	while(2==2) {
		foundMatch=0;		
		currentStream = defaultStream;		
		if( getLineAndTokenize(currentStream,"ExitSolverSettings",localTokenList, numberOfTokens)==1) {
			goto _ExitChoose;
		}
	
		foundMatch=processCommand(localTokenList,currentStream, numberOfTokens); // decendent class first

		if(foundMatch ==0) {	exit(1);} 
		
	} //end of while

	_ExitChoose:
	BETA_OUT << "Exiting Solver Settings reader"<<endl; 
	//defaultStream->close();
	return true;
}
//===================================================================================
int Equations::processCommand(char** localTokenList,ifstream *currentStream, int numToken)
{
	ifstream *originalStream = currentStream;
	char newFileName[80];

	char token[200];
	int  foundMatch=0;

strcpy(token,localTokenList[0]);
BETA_OUT<<"Command processor = 'Equations::processCommand'          "<<token<<endl;

	if( COMPARE(token,"readFromCommandFile")==0 ){
		strcpy(newFileName,localTokenList[1]);
		BETA_OUT<<"Read commands in file:  "<<newFileName<<endl;
		ifstream  newStream(newFileName);
		readCommands(&newStream);  
		OK
		return foundMatch;
	}

	if( COMPARE(token,"openFile")==0 ) {   //This is for this pass only
		//Should be 3 tokens
		strcpy(newFileName,localTokenList[1]);
		currentStream = openFileForInput(newFileName);
		strcpy(token,localTokenList[2]);
		BETA_OUT<<"Input stream for command ***"<<token<<"*** = "<<newFileName<<endl;
		BETA_OUT<<"Stream(openFile) = "<<currentStream<<endl;
		BETA_OUT<<"Stream(openFile) = "<<(*currentStream)<<endl;  
		//OK //think about this
	}

	if(COMPARE(token,"setStorageMethod") ==0) 
	{ 
		setStorageMethod( localTokenList[1] ); 
		OK
	}

	if(COMPARE(token,"setSparseSolverMaxIterations") ==0) 
	{ 
		setMaxIterations(atoi(localTokenList[1]));
		OK
	}

	if(COMPARE(token,"setSparseSolverTolerance") ==0) 
	{ 
		setSOLUTION_EPSILON(atof(localTokenList[1]));
		OK
	}
	
	//+ DG_Feb2004
	if(COMPARE(token,"ReplaceZeroDiagonal") ==0) 
	{ 
		setReplaceZeroDiagonal_Option( atof(localTokenList[1]) );
		OK
	}
	//- DG_Feb2004

	if(COMPARE(token,"setSolverVerboseFlag") ==0) 
	{ 
		setVerboseFlag( (VerboseLevel)(atoi(localTokenList[1])) );
		OK
	}
	if(COMPARE(token,"setReorderingScheme") ==0) 
	{ 
		setReorderingScheme( localTokenList[1] );
		OK
	}

	if(COMPARE(token,"setUsePreviousFillInOrdering") ==0) 
	{
		if((atoi(localTokenList[1]))==0)
			setUsePreviousFillInOrdering(false);
		else
			setUsePreviousFillInOrdering(true);
		OK
	}

	if(COMPARE(token,"setUsePreviousFactorizedMatrix") ==0) 
	{
		if(sm != MKLPardiso){
			cout<<"Can not used Previous Factorized Matrix with this solver (only MKLPardiso)"<<endl;
			exit(1);}
		if((atoi(localTokenList[1]))==0)
			setUsePreviousFactorizedMatrix(false);
		else
			setUsePreviousFactorizedMatrix(true);
		OK
	}

	if(COMPARE(token,"UseMultiCore") ==0) 
	{ 
		BETA_OUT << "Requesting solver to use " << localTokenList[1] << " cores" << endl;
		int numCoresRequestedForSolver=1;
		if(COMPARE(localTokenList[1],"max") ==0)
			numCoresRequestedForSolver=omp_get_max_threads();
		else
			numCoresRequestedForSolver=atoi(localTokenList[1]);

		if( numCoresRequestedForSolver < 1){
			BETA_OUT << "The requested number of cores," << numCoresRequestedForSolver << " is not valid. Exiting..." << endl;
			exit(1);
		}else{
			int MaxPossibleThreads=omp_get_max_threads();
			BETA_OUT << "This machine has " << MaxPossibleThreads << " core(s) that can be potentially used." << endl;
			if(numCoresRequestedForSolver != MaxPossibleThreads){
				BETA_OUT << "The requested number of cores is not the same as the max possible cores." << endl;
			}
			if(numCoresRequestedForSolver > MaxPossibleThreads){ 
				BETA_OUT << "you are asking for more cores than max possible... exiting!";
				exit(1);
			}
			num_threads_for_solver=numCoresRequestedForSolver;
		}
		BETA_OUT << "Solver will attempt to use " << num_threads_for_solver << " thread(s)." << endl;
		OK
	}

	//========
	if(currentStream!=originalStream) {
		BETA_OUT<<"Close alternate input stream"<<endl;
		(*currentStream).close();
		delete currentStream;
		currentStream=originalStream;
		BETA_OUT<<"Automatically reverts to default stream "<< *currentStream
			<<" upon return to readCommands"<<endl;
	}

	if(foundMatch ==0){	
		BETA_OUT << "Command: " << token << endl;
		BETA_OUT << "Did not recognize command: " <<token<<endl;
		exit(1);
	}

	return foundMatch;
}
//========================================================================
/////////////////////////////////////////////////////////////////////////////
//
//  S T U F F	F O R	P A R A L L E L
//

void Equations::assignElementsToThreads(const int *map, const int size,Array<threadWorkspace> &tList,int elNum)
{ 
	int I,J,i,j,k,klimit,m,mlimit,ip,jp;
	LargeMatrix &Kg = *a;
	bool *rest = restraint;
	for(i=0;i<size;i++) {
		I=map[i];
		for(j=0;j<size;j++) {
			J=map[j];
			if(rest[J]) {
				if(I==J) processDofForThreadAssignment(tList,I,elNum);
				if(additionalEquation[J]) {
					AdditionalEquation &aMj = *additionalEquation[J];
					if(rest[I]) {
						if(additionalEquation[I]) {
							AdditionalEquation &aMi = *additionalEquation[I];
							for(k=0,klimit=aMj.numTerms,mlimit=aMi.numTerms;k<klimit;k++) {
								ip = aMj.master[k];
								if(rest[ip]) {
									processDofForThreadAssignment(tList,ip,elNum);
								}
								else for(m=0;m<mlimit;m++) {
									jp = aMi.master[m];
									if(!rest[jp] && ip >= jp) {
										processDofForThreadAssignment(tList,ip,elNum);
								//		processDofForThreadAssignment(tList,jp,elNum);
									}
								}
							}
						}
					} else
						for(k=0;k<aMj.numTerms;k++)
							if(!rest[aMj.master[k]] && I >= aMj.master[k]) {
								processDofForThreadAssignment(tList,I,elNum);
							//	processDofForThreadAssignment(tList,aMj.master[k],elNum);
							}
				}
			} else if(rest[I]) {  // If I is restrained and J is not
				if(additionalEquation[I]) {
					AdditionalEquation &aMi = *additionalEquation[I];
					for(k=0,klimit=aMi.numTerms;k<klimit;k++) {
						ip = aMi.master[k];
						if(rest[ip]) processDofForThreadAssignment(tList,ip,elNum);
						else if(ip >= J) {processDofForThreadAssignment(tList,ip,elNum);
							//			processDofForThreadAssignment(tList,J,elNum);

						}
					}
				}
			} else if(I>=J) { // If neither restrained do normal assembly
				processDofForThreadAssignment(tList,I,elNum);
			//	processDofForThreadAssignment(tList,J,elNum);
				//assumes, symm. matrix is lower diagonal. large matrix object should handle it specifically if it is upper diagonal.
				//if the matrix is general, the matrix should handle specifying (J,I) as non-zero location
			}
		}
	}

}

void Equations::processDofForThreadAssignment(Array<threadWorkspace> &tList,int dof,int elNum)
{
	int numThreads = numberOfOpenMPThreads;

	for(int i=0;i<numThreads;i++)
	{
		int min = tList[i].minEquationNumber;
		int max = tList[i].maxEquationNumber;

		if (dof >= min && dof <=max)
		{
			tList[i].elementList.add(elNum);
			tList[i].secondaryElementList.add(elNum);
			return;
		}
		
	}
}


void Equations::addToLoadVector_Parallel(const int *map, const int size, 
	const double *force, threadWorkspace *thread)
{
	if(force==0) return;
	int i,I,k;

	double *F = NULL;
	F = loadVector;

	int maxEqNum = thread->maxEquationNumber;
	int minEqNum = thread->minEquationNumber;

	WhoMethod("addToLoadVector(const int *, const int, const double *)");
	if(!allocated) FatalError("Equations not allocated yet.");
	for(i=0;i<size;i++) {
		I = map[i];
		RangeCheck(0,I,numberOfEquations-1);
//JV took out code related to addnl eqns to make code simpler
		if(additionalEquation[I]) {
			AdditionalEquation &aMi = *additionalEquation[I];
			for(k=0;k<aMi.numTerms;k++){
				if(!restraint[aMi.master[k]]){
					if(aMi.master[k] >=minEqNum && aMi.master[k]<=maxEqNum){
						F[aMi.master[k]] += force[i] * aMi.factor[k]; }
				}
			}
		}
		else if(!restraint[I]) 
		{
			if(I >= minEqNum && I <=maxEqNum) 
				{F[I] += force[i];}
		}
	}
}


//======================================================================================


void Equations::addMatrix_Parallel(const Matrix &matrix, const int *map, const int size, 
		const double *force, threadWorkspace *thread)
{
	WhoMethod("addMatrix(const Matrix &, const int *, const int, const double *)");
	if(!allocated) FatalError("Equations not allocated yet.");

	//JV091104 - diverting to another function if the solver is MKLGBMAT or LAPACKGB
	if (sm==mklgbmat || sm==LapackGB) {
		addMatrixMKL(matrix, map, size, force);
		return;
	}

	//JV091104 - diverting to another function if the solver is PetSc
	if (sm==PetSc || sm==PetScUnsymm) {
		addMatrix2(matrix, map, size, force);
		return;
	}

	// This assembly process is fairly complicated. Be sure to make sure of what
// you are doing if you make any changes!
// 'matrix' should be in upper triangluar format.
// 'Kg' Matrix should behave like normal lower triangular matrix.
// Copying object variables to local variables.

	int i,j,k,m,I,J,ip,jp;
	double matrix_ij=0;
	double *F =0;
	bool *rest = restraint;

	F = pseudoLoadVector;


	LargeMatrix *Mat;
	Mat = a;
	LargeMatrix &Kg = *Mat;

	int minEq = thread->minEquationNumber;
	int maxEq = thread->maxEquationNumber;


// Assemble upper triangle of Ke and any force contribution
	for(i=0;i<size;i++) {
		I = map[i];                                  //selects global DOF corresponding to local DOF i

		RangeCheck(0,I,numberOfEquations-1);

		// Add in force if force vector exists and dof not restrained
		if(force!=0) {
            if(!rest[I]){ 
				if (I <= maxEq && I >= minEq) {
					F[I] += force[i];   }                              //If DOF I is not restrained simply add in force contribution from Fe}
            }else if(additionalEquation[I]) {							   //IF DOF I is an MPC modify the force component of the master DOF equations ... why do we multiply a force by a factor??
				AdditionalEquation &aMi = *additionalEquation[I];          //IF DOF I is a simple restraint nothing is done
				for(k=0;k<aMi.numTerms;k++){
						if (aMi.master[k] <= maxEq && aMi.master[k] >= minEq){
							F[aMi.master[k]] += force[i] * aMi.factor[k];} 
				}
			} 
		}



		for(j=0;j<size;j++) {
			if(i<j) matrix_ij = matrix(i,j);             // this conditional statements ensure the upper triangular portion of Ke is used
			else matrix_ij = matrix(j,i);
			J=map[j];									 //selects global DOF corresponding to local DOF j
			RangeCheck(0,J,numberOfEquations-1);


			// If J restrained move to lhs and all J restrained cases
			if(rest[J]) {
				if(I==J) if (I <= maxEq && I >= minEq) {
					Kg(I,I)=1; }     //the diagonal of Equation/Row J to 1, all other entries of the row in Kg are 0

				if(!rest[I]) {           //if J is restrained but I is not
					if(additionalEquation[J]){ if (I <= maxEq && I >= minEq) {
						F[I] -= matrix_ij * additionalEquation[J]->constant;}
					//cout<<"[[[[[[[[[[[[[[[||||"; //added by bco
					}//if J has additional eq. move (K_IJ * MPC Const) to Fg 
					else 
						if (I <= maxEq && I >= minEq){
							F[I] -= matrix_ij * solution[J];}     //if J does not have additional eq. simply    move K_IJ*u_J to Fg               
				}
				
				else if(additionalEquation[I]) {   //if J is restrained and I has additional eq. (MPC)
					AdditionalEquation &aMi = *additionalEquation[I];
					if (aMi.numTerms==0) //cout<<"ZERO!!!!"; //added by bco
					for(k=0;k<aMi.numTerms;k++) {
						if (aMi.master[k] <= maxEq && aMi.master[k] >= minEq){
						if(additionalEquation[J])  //if J and I both have additional eqs.
						   //cout<<"<<<<<<<<<<<<||";  //added by bco
							F[aMi.master[k]] -= matrix_ij * aMi.factor[k] *
								additionalEquation[J]->constant;
						else
							F[aMi.master[k]] -= matrix_ij * aMi.factor[k] *
								solution[J];
						}
					}
				}

				if(additionalEquation[J]) {
					AdditionalEquation &aMj = *additionalEquation[J];
					if(rest[I]) {
						if (additionalEquation[I]) {
							//cout<<"||>>>>>>>>>>>>>>"<<endl; //added by bco
							AdditionalEquation &aMi = *additionalEquation[I];
							for(k=0;k<aMj.numTerms;k++) {
								ip = aMj.master[k];
								if(rest[ip]) {
									if (ip <= maxEq && ip >= minEq){
										Kg(ip,ip) = 1;
									}
								}
								else 
								{
									for(m=0;m<aMi.numTerms;m++) {
									jp = aMi.master[m];
									if(rest[jp]) {
//										cout << "master " <<jp<<" is restrained."<<endl;
										if (ip <= maxEq && ip >= minEq){
										F[ip] -= matrix_ij * aMj.factor[k] *
												aMi.factor[m] * solution[jp];
										}
									}
									
									else if(ip >= jp){
										if (ip <= maxEq && ip >= minEq){
										Kg(ip,jp) += matrix_ij * aMj.factor[k] * 
												aMi.factor[m];
										//cout<<"WRITING TO LOCATION(a): "<<ip<<"    "<<jp<<endl;
										}
									}
									}
								}
							}
						}
					}
					
					else {
						//cout<<"||||]]]]]]]]]]]]]]]"<<endl; //add by bco
						for(k=0;k<aMj.numTerms;k++) {
							jp = aMj.master[k];
							if(rest[jp]) {
//								cout << "master " <<jp<<" is restrained."<<endl;
								if (I <= maxEq && I >= minEq){
									F[I] -= matrix_ij * aMj.factor[k] * solution[jp];}
							}
							
							else if(I >= jp) 
								if (I <= maxEq && I >= minEq)
								{ Kg(I,jp) += matrix_ij * aMj.factor[k];
								//cout<<"WRITING TO LOCATION(b): "<<I<<"    "<<jp<<endl;
								}
						}
					}
				} 
			}
			//--------------------------------------------------------------------------------------------------
			else if(rest[I]) {  // If I is restrained and J is not
				if(additionalEquation[I]) {
					AdditionalEquation &aMi = *additionalEquation[I];
					for(k=0;k<aMi.numTerms;k++) {
						ip = aMi.master[k];
						if(rest[ip]) {
							if (ip <= maxEq && ip >= minEq){ Kg(ip,ip) = 1;}
						}
						else if(aMi.master[k] >= J){
							if (aMi.master[k] <= maxEq && aMi.master[k] >= minEq){
							Kg(aMi.master[k],J) += matrix_ij * aMi.factor[k];
							//cout<<"WRITING TO LOCATION(c): "<<aMi.master[k]<<"    "<<J<<endl;
							}
						}
					}
				} 
			}
			
			else if(I>=J) // If neither restrained do normal assembly
				if (I <= maxEq && I >= minEq)
				{Kg(I,J) += matrix_ij;		}
		}
	
}

	#ifdef DEBUG_Level3
	Kg.printMatrix("Kg = ");
	#endif

}



