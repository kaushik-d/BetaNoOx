#include "stdafx.h"
#include "FileManager.hpp"

#include <iomanip>
#include <sys/types.h>
#include <sys/timeb.h>

using namespace std;
ofstream *timelog=0;

//#define FORMAT setw(10) << setprecision(3)
#define FORMAT ""

class Stopwatch {
    public:
      timeb starttime, lasttime;
      Stopwatch() {reset();}
	  void reset(){
		  ftime(&lasttime); 
		  ftime(&starttime); 
	  }
};

Stopwatch stopwatch;

//======================================================================
void OpenTimeLogFile(FileManager *filemanager)
{
	if(timelog==0){
		timelog=filemanager->OpenOutputStream("timelog.txt");
		timelog->setf(ios::fixed);
		stopwatch.reset();
	}
}
//======================================================================
void displayTime(const char * message, ostream &outStream)
 {
	 if(timelog==0){
		 cerr<<"timelogfile not opened! Use void OpenTimeLogFile(FileManager *filemanager)" << endl;
	 }

	 timeb newtime;
	 ftime(&newtime);
	 double t1, t2;

	if(strcmp(message, "") != 0 )
	{
	t1= double(newtime.time - stopwatch.starttime.time) + double(newtime.millitm - stopwatch.starttime.millitm)/1000.0;
	t2= double(newtime.time - stopwatch.lasttime.time)  + double(newtime.millitm - stopwatch.lasttime.millitm)/1000.0;

	outStream << message;
	outStream << resetiosflags(ios::floatfield);
	outStream << setiosflags(ios::fixed);
	outStream << FORMAT;
	outStream << "....Elapsed, total time (seconds) = "
						<< t2 << "  " << t1 << endl;
	outStream << resetiosflags(ios::floatfield);
	 if(timelog!=0){
		*timelog << message;
		*timelog << FORMAT;
		*timelog << "....Elapsed, total time (seconds) = "
						<< t2 << "  " << t1 << endl;
	 }

	}
 stopwatch.lasttime = newtime;
}
//======================================================================

/*
#include <time.h>
//======================================================================
void displayTime(char * message)
 {
 static time_t lastTime = time(NULL);
 static time_t startTime = time(NULL);

 time_t newTime=time(NULL);
 
 if(strcmp(message, "") != 0 )
 {
 BETA_OUT<< message;
 //Logfile << "the time is " << ctime(&newTime)<<endl ;
 BETA_OUT << setiosflags(ios::fixed);
 BETA_OUT << "....Elapsed, total time (seconds) = "
	                    << FORMAT <<difftime( newTime, lastTime ) <<"  "
	                    << FORMAT <<difftime( newTime, startTime )<<endl;
 }
 lastTime = newTime;
 }
//======================================================================
*/

