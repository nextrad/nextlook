#include "timer.hpp"


void Timer::start(void)
{
	startTime = clock();
}


float Timer::getTimeElapsed(void)
{
	return ((float)clock() - (float)startTime)/CLOCKS_PER_SEC;
}


clock_t Timer::getStartTime(void)
{
	return startTime;
}

void Timer::setStartTime(clock_t start)
{
	startTime = start;
}
