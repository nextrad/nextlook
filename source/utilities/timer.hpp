#ifndef TIMER_HPP
#define TIMER_HPP

#include <time.h>

class Timer 
{
	private:
		clock_t startTime;
	
	public:		
		//constructor		
		Timer(void){};
		
		//functions
		void start(void);
		float getTimeElapsed(void);
		clock_t getStartTime(void);
		void setStartTime(clock_t start);
};

#endif
