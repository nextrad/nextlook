#include "logger.hpp"


Logger::Logger()
{
	
}


void Logger::write(std::string message)
{
	std::cout << message << std::endl;	
}


void Logger::write(std::string message, Timer &timer)
{
	std::cout << std::setprecision(3) << std::fixed << timer.getTimeElapsed() << "s: \t" << message << std::endl;	
}
