#include "logger.hpp"


Logger::Logger()
{
	
}


void Logger::write(string message)
{
	cout << message << endl;	
}


void Logger::write(string message, Timer &timer)
{
	cout << setprecision (3) << fixed << timer.getTimeElapsed() << "s: \t" << message << endl;	
}
