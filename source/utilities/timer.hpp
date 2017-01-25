#ifndef TIMER_HPP
#define TIMER_HPP

#include <time.h>
#include <sys/time.h>

class Timer 
{
	private:
		double startTime;
		struct timeval time;
		
		double getTimeInSec(void);
	
	public:		
		//constructor		
		Timer(void){};
		
		//functions
		void start(void);
		double getTimeElapsed(void);
		clock_t getStartTime(void);
		void setStartTime(clock_t start);
};

#endif
