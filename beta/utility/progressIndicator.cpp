#include "stdafx.h"

#include <string.h>
#include "progressIndicator.hpp" 
using namespace std;


void ProgressIndicator::start(const char *message)
{
	int i;
	column = 0;
	cerr << '\n' << message;
	cerr << '\n' << '[';
	for(i=0;i<50;i++) cerr << ' ';
	cerr << ']';
	for(i=0;i<51;i++) cerr << (char)0x08;
	cerr << ends;
	startTime = time(NULL);
	return;
}

void ProgressIndicator::notifyPercentDone(const float val)
{
	int i,limit,newColumn = (int)val/2;
	if(newColumn>50) newColumn=50;
	for(i=newColumn,limit=column;i<limit;i++) cerr << (char)0x08 << ends;
	for(i=column;i<newColumn;i++) cerr << '.' << ends;
	column = newColumn;
}
	
void ProgressIndicator::done(const char *message, bool reportTime)
{
	double elapsedTime;
	endTime = time(NULL);
	notifyPercentDone(100.0f);
	cerr << '\n' << message;
	elapsedTime=difftime(endTime,startTime);
	if(reportTime) cerr << " (" << elapsedTime << " secs.)";
	if (elapsedTime > 60)
		cerr << " (" << elapsedTime/60 << " mins.)";
	if (elapsedTime > 3600)
		cerr << " (" << elapsedTime/3600 << " hours.)";
	cerr << '\n';
}

#ifdef DEBUG_CommandReader

#include <fstream>
#include <stdio.h>
void main(void)
{
	cout << "Starting..." << endl;
	int i,j;
	CommandReader *cr;
   
	ifstream fin("command.dat");
	cr = new CommandReader(fin);

	StrCC s = " 1       2 3   4   ";
	StrCC *s1;

	cout << "Split Test:\n";
	j = s.split(s1);
	for(i=0;i<j;i++)
		cout << s1[i] << endl;

	cout << endl;
	for (j=0;j<4;j++) {
		cout << "Read Status: "; cout.flush(); cout << cr->getCommand() << endl;
		cout << "Command: " << cr->command << "\nArguments: " << endl;
		for(i=0;i<cr->args.numArgs;i++) {
			cout << i << " : " << cr->args(i) << endl;
		}
		cout << endl;
		getchar();
	}
	for (j=0;j<3;j++) {
		if (cr->getList() == commandReader_error) 
			cout << "Error: List Read" << endl;
		else {
			cout << "List: \n";
			for(i=0;i<cr->list.numItems;i++) {
				cout << i << " : " << cr->list(i) << '\n';
			}
		}
		cout << endl;
		getchar();
	}
	cout << endl;
	delete cr;
	CommandReader *cr2;
	cr2 = new CommandReader("command2.dat");
	for (j=0;j<4;j++) {
		cout << "Read Status: " << cr2->getCommand() << '\n';
		cout << "Command: " << cr2->command << "\nArguments: \n";
		for(i=0;i<cr2->args.numArgs;i++) {
			cout << i << " : " << cr2->args(i) << '\n';
		}
		cout << endl;
		getchar();
	}
	for (j=0;j<3;j++) {
		if (cr2->getList() == commandReader_error) 
			cout << "Error: List Read" << endl;
		else {
			cout << "List: \n";
			for(i=0;i<cr2->list.numItems;i++) {
				cout << i << " : " << cr2->list(i) << '\n';
			}
		}
		cout << endl;
		getchar();
	}
	
	cout << endl;

	float val;
	cr2->notifyPercentDone(-1);
	for(val=0;val<100;val+=1) {
		cr2->notifyPercentDone(val);
		for(i=0;i<200000;i++);
	}
	delete cr2;

}
#endif

