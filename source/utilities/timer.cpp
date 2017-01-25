#include "timer.hpp"


double Timer::getTimeInSec(void)
{
	gettimeofday(&time,NULL);    
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void Timer::start(void)
{
	startTime = getTimeInSec();
}


double Timer::getTimeElapsed(void)
{    
    gettimeofday(&time,NULL);         
	return getTimeInSec() - startTime;
}


clock_t Timer::getStartTime(void)
{
	return startTime;
}

void Timer::setStartTime(clock_t start)
{
	startTime = start;
}
