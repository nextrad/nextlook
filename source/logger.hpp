#ifndef LOGGER_HPP
#define LOGGER_HPP

#include "timer.hpp"
#include "includes.hpp"

class Logger
{
	private:
		
		
	public:
		Logger(void);
		void write(string message);
		void write(string message, Timer &timer);
		
};

#endif
