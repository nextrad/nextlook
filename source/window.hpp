#ifndef WINDOW_HPP
#define WINDOW_HPP

#include "includes.hpp"

enum WindowFunction{HANNING, HAMMING, UNIFORM, BLACKMAN};

class Window
{
	private:
		WindowFunction function;
		float* window;		
		
	public:		
		Window(WindowFunction windowFunction, int windowSize);		
		float getSample(int sampleNumber);
        double sinc(double x);
};

#endif
