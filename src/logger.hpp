#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>
#include <iostream>
#include <iomanip> 

#include "timer.hpp"

class Logger
{
	private:
		
		
	public:
		Logger(void);
		void write(std::string message);
		void write(std::string message, Timer &timer);
		
};

#endif
