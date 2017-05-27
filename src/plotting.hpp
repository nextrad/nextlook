#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <opencv2/opencv.hpp>
#include <fftw3.h>

#include "parameters.hpp"
#include "logger.hpp"

#define COLOUR_MAP_WIDTH 30
#define COLOUR_MAP_PATH "../colour_maps/"

//Functions
void plotWaterfall(void);
void plotDoppler(void);
void updateWaterfall(int rangeLine, uint8_t  *imageValues);
void updateDoppler(uint8_t  *imageValues);
void initOpenCV(void);
static void updatePlots(int, void*);


extern int dopplerThresholdSlider;
enum plotType {NORMAL, FFT_SHIFT};
enum plotStyle {AMPLITUDE, IQ};

class Plot 
{
	private:
		std::string title;
		std::string xlable, ylable;
		Logger logger;
	public:
		Plot(void);
		void gnuPlot(fftw_complex *array, char const *plotTitle, plotType type = NORMAL, plotStyle style = AMPLITUDE);	
		void gnuPlot(uint8_t *array, char const *plotTitle);
};

#endif
