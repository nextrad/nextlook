#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <opencv2/opencv.hpp>
#include <boost/thread.hpp>
#include <fftw3.h>
#include <string.h>

#include "logger.hpp"
#include "experiment.hpp" 

/*
#define COLOUR_MAP_WIDTH 30
#define COLOUR_MAP_PATH "../colour_maps/"*/

enum plotType {NORMAL, FFT_SHIFT};
enum plotStyle {AMPLITUDE, IQ};

class OpenCVPlot 
{
	private:
		const int cmap_width = 30;
		const char* cmap_path = "../colour_maps/";	
		
		const int slowMax = 500;		
		const int histogramMax = 1;
		const int colourMapMax = 11;	
		const int thresholdMax = 255;
		const int averagingMax = 100;

		int slowSlider;		
		int	thresholdSlider;
		int histogramSlider;	
		int	averagingSlider;	
		int dopplerColourMapSlider;
		int waterfallColourMapSlider;
		
		//holds the index of the current doppler plot 
		int doppler_plot_index;
		
		cv::Mat waterImage;
		cv::Mat doppImage;
		
		cv::Mat resizedWaterImage;
		cv::Mat processedWaterImage;	
		
		std::vector<cv::Mat> dopplerMatrix;
		cv::Mat averagedDopplerImage;
		cv::Mat scaledDopplerImage;
		cv::Mat processedDopplerImage;
		
		cv::Size waterSize;
		cv::Size doppSize;

		Experiment* experiment;
	public:
		OpenCVPlot(Experiment* exp);
		void plotWaterfall(void);
		void plotDoppler(void);
		void addToWaterPlot(int rangeLine, double  *imageValues);
		void addToDopplerPlot(int dopplerLine, double *imageValues);
		void initOpenCV(void);
		void savePlots(void);
		static void updatePlots(int, void*);
};

class GNUPlot 
{
	private:
		std::string title;
		std::string xlable, ylable;
		Logger logger;		
		Experiment* experiment;
	public:
		GNUPlot(Experiment* exp);
		void gnuPlot(fftw_complex *array, char const *plotTitle, Experiment* exp, plotType type = NORMAL, plotStyle style = AMPLITUDE);	
		void gnuPlot(uint8_t *array, char const *plotTitle, Experiment* exp);
};

#endif
