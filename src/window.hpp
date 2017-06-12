#ifndef WINDOW_HPP
#define WINDOW_HPP

#include <stdlib.h>
#include <math.h>

enum WindowFunction {HANNING, HAMMING, UNIFORM, BLACKMAN};

class Window
{
	private:
		WindowFunction function;
		float* window;		
		
	public:		
		Window(void);		
		void init(WindowFunction windowFunction, int windowSize);
		float getSample(int sampleNumber);
        double sinc(double x);
};

#endif
