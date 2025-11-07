#include "stdafx.h"

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "utility/excepts.hpp"
#include "equation.hpp"
#undef PRINT
//#define PRINT(a) cout << __FILE__ << '[' << __LINE__ << "]: " << #a << " = " << a << endl;
#define PRINT(a)

CreateErrorHandler(AdditionalEquation);
AdditionalEquation::AdditionalEquation()
{
	WhoClass("AdditionalEquation");
	WhoMethod("AdditionalEquation()");
	numTerms=0;
	constant=0;
	slave=-1; //JV112807 - initializing ALL variables
}

AdditionalEquation::~AdditionalEquation()
{
	WhoMethod("~AdditionalEquation()");
}

bool AdditionalEquation::operator!=(const AdditionalEquation &ae) const
{
	WhoMethod("operator!=(const AdditionalEquation &) const");
	int i,j;
	if(numTerms!=ae.numTerms) return true;
	if(fabs(ae.constant-constant) > (fabs(ae.constant)+fabs(constant))/20000) return true;
	for(i=0;i<numTerms;i++) {
		for(j=0;j<ae.numTerms;j++) {
			if(ae.master[j] == master[i]) {
				if(fabs(ae.factor[j]-factor[i]) > 
						(fabs(ae.factor[j])+fabs(factor[i]))/20000)
					return true;
				else j = ae.numTerms;	
			} 
		}
	}
	return false;
}

void AdditionalEquation::divideRHSby(const double &val)
{
	WhoMethod("divideRHSby(const double &)");
	if(val==0) {
		Warning("Divide by zero.");
		return;
	}
	for(int i=0,limit=numTerms;i<limit;i++) factor[i]/=val;
	constant /= val;
}

AdditionalEquation::Result AdditionalEquation::collectTerms(void)
{
	WhoMethod("collectTerms(void)");
	// Terms on rhs have been collected as they were added
	int i,limit = numTerms;
	for(i=0;i<limit;i++) {
		if(factor[i]==0.0) {
			removeTerm(i);
			limit = numTerms;
		}
	}
	double slaveFactor=1.0;
	for(i=0;i<limit;i++) {
		if(master[i] == slave) {
			slaveFactor -= factor[i];
			removeTerm(i);
			limit = numTerms;
			if(slaveFactor == 0.0) { // Pick a new slave
				if(numTerms == 0) {
					if(fabs(constant) > 1e-8) {slave=-1; return Bad;}
					else return Trivial;
				}
				slave = master[0];
				slaveFactor = -factor[0];
				removeTerm(0);
				limit = numTerms;
			}
			divideRHSby(slaveFactor);
			return Good;
		}
	}
	return Good;
}

AdditionalEquation::Result AdditionalEquation::apply(AdditionalEquation &eqn)
{
	//this function plugs 'eqn' into 'this' equation (if possible) 
	WhoMethod("apply(AdditionalEquation &)");
	int i,limit=numTerms,j;
	AdditionalEquation::Result result;
	// Replace masters if they are slaved
	for(i=0;i<limit;i++) {
		if(master[i]==eqn.slave) {
			double oldFactor = factor[i];
			if(eqn.numTerms) {
				replaceTerm(i,eqn.master[0],eqn.factor[0]*oldFactor);
				for(j=1;j<eqn.numTerms;j++)
					addAdditionalTerm(eqn.master[j],eqn.factor[j]*oldFactor);
			} else {
				removeTerm(i);
			}
			constant += eqn.constant*oldFactor;
			limit = numTerms;
			i=limit;
		}
	}
	//consolidate 'this' equation
	result = collectTerms();
	if(result!=Good) {
		if(result==Trivial) return Trivial;
		return Bad;
	}
	PRINT(eqn);
	// If the eqn.slave is the same as slave then rewrite our equation if we can
	if(eqn.slave==slave) {
		if(numTerms) {
			double oldFactor = factor[0];
			slave = master[0];
			factor[0] = -1.0;
			divideRHSby(-oldFactor);
			oldFactor = factor[0];
			if(eqn.numTerms) {
				PRINT(*this);
				replaceTerm(0,eqn.master[0],eqn.factor[0]*oldFactor);
				PRINT(*this);
				for(j=1;j<eqn.numTerms;j++)
					addAdditionalTerm(eqn.master[j],eqn.factor[j]*oldFactor);
			} else {
				removeTerm(0);
			}
			constant += eqn.constant*oldFactor;
			result = collectTerms();
			if(result!=Good) {
				if(result==Trivial) return Trivial;
				return Bad;
			}
			if(numTerms==0 && fabs(constant)<1e-8) return Good;
		} else { // Modify eqn to be (eqn.slave = this->constant) and this to 
					// eqn.master[0] = eqn.master[1]*eqn.factor[1]/-eqn.factor[0] ...
			if(eqn.numTerms>0) {
				double tmp = constant;
				slave = eqn.master[0];
				for(j=1;j<eqn.numTerms;j++) 
					addAdditionalTerm(eqn.master[j],eqn.factor[j]/(-eqn.factor[0]));
				constant = (eqn.constant-constant)/(-eqn.factor[0]);
				eqn.numTerms=0;
				eqn.master.clear();
				eqn.factor.clear();
				eqn.constant = tmp;
			} else {
				if(fabs(eqn.constant-constant)>1e-8) return Bad;
			}
		}
	}
	return Good;
}

bool AdditionalEquation::operator==(const AdditionalEquation &ae) const
{ return ((operator!=(ae))==true)?false:true;}

bool AdditionalEquation::addAdditionalTerm(const int newMaster, const
		double newFactor)
{
	WhoMethod("addAdditionalTerm(const int, const double)");
	int i;

	if(fabs(newFactor)<1e-8) 
		return true;
	
	//cerr << "numTerms " << numTerms <<endl;
	if(numTerms<0){
		cerr << "numTerms<0 in AdditionalEquations! major bug!"<<endl;
		exit(1);
	}

	// If master already exists just add in factor
	for(i=0;i<numTerms;i++)
		if(master[i] == newMaster) {
			factor[i] += newFactor;
			if(fabs(factor[i])<1e-8) {
//				cout << "(fabs(factor[i])<1e-8)" << endl;
				removeTerm(i);      //Fix!! <-- JV: who put this comment??
				return true;//JV052008 need to return after removing term. continuing will lead to reading invalid information
			}
            //if(fabs(factor[i])>1e-8 && fabs(factor[i])<1e-4 ) { //there is a change the previous removeTerm(i) 
			// will release the memory from factor so, we need to add a check to catch that
            if(factor.size() && fabs(factor[i])>=1e-8 && fabs(factor[i])<1e-4 ) {
				cout<<"Missed term: "<<factor[i]<<endl ;
			    cout<<"newMaster = " <<newMaster<<endl;
			}
			return true;
		}
	numTerms++;
	master.push_back(newMaster);
	factor.push_back(newFactor);
	return true;
}

void AdditionalEquation::replaceTerm(const int termNum, const int newMaster, 
		const double &newFactor)
{
	WhoMethod("replaceTerm(const int, const int, const double &)");
	int i;
	// If master already exists just add in factor
	for(i=0;i<numTerms;i++) {
		if(i==termNum) continue;
		if(master[i] == newMaster) {
			factor[i] += newFactor;
			removeTerm(termNum);
			return;
		}
	}
	master[termNum] = newMaster;
	factor[termNum] = newFactor;
}

void AdditionalEquation::freeMasters(void)
{
	master.clear();
	factor.clear();
	numTerms=0;
}

void AdditionalEquation::removeTerm(int termNum)
{
	WhoMethod("removeTerm(int)");
	static bool firsttime=false;
	if(!firsttime){
		//this is added to avoid flooding cout with this statement. this way it is printed only once.
		cout << "suppressing AdditionalEquation::removeTerm(int termNum) output...uncomment code in source if needed." <<endl;
		//cout << "removeTerm: " << termNum << " from " << (*this) << endl;
		firsttime=true;
	}
//	cout << "removeTerm: " << termNum << " from " << (*this) << endl;

	if(termNum >= numTerms) {
		cout << "invalid input parameter to AdditionalEquation::removeTerm(int termNum)" << endl;
		exit(1);
		return;
	}
	numTerms--;
	master.erase(master.begin()+termNum);
	factor.erase(factor.begin()+termNum);
}

double AdditionalEquation::resolve(const double *solution) const
{
	int i;
	double val=constant;
	for(i=0;i<numTerms;i++) val += factor[i] * solution[master[i]];
	return val;
}
/*
double AdditionalEquation::resolveResidual(const double *solution) const
{
	int i;
	double val=0;
	for(i=0;i<numTerms;i++) val += factor[i] * solution[master[i]];
	return val;
}
*/

ostream &operator << (ostream &o, const AdditionalEquation &ae)
{
	int i;
	o << '[' << ae.slave/3 << '.' << ae.slave%3 << "] = ";
	for(i=0;i<ae.numTerms;i++) o << ae.factor[i] << " * [" << ae.master[i]/3 << '.' << ae.master[i]%3 << "] + ";
	o << ae.constant;
	return o;
}
