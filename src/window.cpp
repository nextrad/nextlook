#include "window.hpp"

Window::Window(void)
{

}

void Window::init(WindowFunction windowFunction, int windowSize)
{
	function = windowFunction;
    window = (float*)malloc(windowSize*sizeof(float));
    double commonTerm;
    const double pi = 3.1415926;

	for (int i = 0; i < windowSize; i++)
	{

        float common = (2*pi*i)/(windowSize - 1);

        switch(windowFunction)
		{
        case HANNING  : window[i] = 0.50 - 0.50*cos(common); break;
        case HAMMING  : window[i] = 0.54 - 0.46*cos(common); break;
        case BLACKMAN : window[i] = 0.42659 - 0.49656*cos(common) + 0.076849*cos(2*common); break;
        case UNIFORM  : window[i] = 1;
		}	
	}
}

float Window::getSample(int sampleNumber)
{
	return window[sampleNumber];
}

double Window::sinc(double x)
{
   if (x == 0) return 1;
   return sin(x)/x;
}
