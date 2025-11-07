#include "stdafx.h"

#include "threadWorkspace.hpp"

//===================================================================================
	threadWorkspace::threadWorkspace()
	{
		eWorkspace=0;
		filemanager=0;
	};
//===================================================================================
	threadWorkspace::~threadWorkspace()
	{

		if(filemanager)
			delete filemanager;
			
		if(eWorkspace)
			delete eWorkspace;

	};
//===================================================================================
	void threadWorkspace::initialize()
	{
		eWorkspace = 0;
		filemanager = 0;
	};
//===================================================================================
	void threadWorkspace::setThreadFileName(int i)
	{
		stringstream number;
		number << i;
		threadFileName = "thread" + number.str() + ".txt";
	};
//===================================================================================
	void threadWorkspace::setDefaultOutputStreamForFileManager()
	{
		filemanager->setDefaultOutputStream(threadFileName);
	};
//===================================================================================
	void threadWorkspace::printThreadAssignments()
	{
//		(*filemanager->outStream)<<"Elements Assigned to this Thread for Calculation:     "<<minElementCalculationNumber<<"  to  "<<maxElementCalculationNumber<<endl<<endl;
		(*filemanager->outStream)<<"Equations Assigned to this Thread:     "<<minEquationNumber<<"  to  "<<maxEquationNumber<<endl<<endl;
		(*filemanager->outStream)<<"Elements Assigned to this Thread:"<<endl<<endl;
		int numElementInList = elementList.getNum();
		for(int i=0; i<numElementInList; i++) {
			(*filemanager->outStream) << elementList[i] <<"     "<< secondaryElementList[i] << endl;
		}
	};
//==================================================================================
	void threadWorkspace::printThreadCalculationAssignments(){

		(*filemanager->outStream)<<"Elements Assigned to this Thread for Calculation:     "<<minElementCalculationNumber<<"  to  "<<maxElementCalculationNumber<<endl<<endl;
		
	};
//===================================================================================
	bool threadWorkspace::isThreadFinishedF()
	{
		for(int i=0; i<secondaryElementList.getNum();i++)
		{
			if(secondaryElementList[i] != -1)
				return false;
		}

		return true;
	};
//===================================================================================
	bool threadWorkspace::isThreadFinishedK()
	{
		for(int i=0; i<elementList.getNum();i++)
		{
			if(elementList[i] != -1)
				return false;
		}

		return true;
	};
//===================================================================================