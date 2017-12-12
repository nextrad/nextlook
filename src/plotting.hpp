#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <opencv2/opencv.hpp>
#include <boost/thread.hpp>
#include <fftw3.h>
#include <string.h>

#include "logger.hpp"
#include "experiment.hpp" 


// Abbreviations List
// rt: range-time
// rd: range-doppler

enum plotType {NORMAL, FFT_SHIFT};
enum plotStyle {LOG, MAG, IQ};

class OpenCVPlot 
{
	private:
		const int cMapWidth = 30;
		std::string cMapRoot = "../cmaps/";	
		
		const int slowMax = 500;		
		const int histMax = 1;
		const int cMapMax = 11;	
		const int thrsMax = 255;
		
		//set by n_range_lines/update_rate i.e. max number of updates possible.
		int avrgMax;

		int slowSldr;		
		int	thrsSldr;
		int histSldr;	
		int	avrgSldr;	
		int cmapSldr;
		
		//holds the index of the current doppler plot 
		int rdIndex;
		
		cv::Mat rtImage;
		cv::Mat rtImageResize;
		cv::Mat rtImage8bit;	
		
		cv::Mat rdImage;
		cv::Mat rdImageAvg;
		cv::Mat rdImage8bit;
		std::vector<cv::Mat> rdVector;
		
		cv::Mat spImage;
		cv::Mat spImageResize;
		cv::Mat spImage8bit;
		
		cv::Size rtSize;
		cv::Size rdSize;
		cv::Size spSize;

		Experiment* experiment;
	public:
		OpenCVPlot(Experiment* exp);
		void plotRTI(void);
		void plotRD(void);
		void plotSP(void);
		void addRTI(int rangeLine, double  *imageValues);
		void addRD(double *imageValues);
		void addSP(double *imageValues);
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
	public:
		void plot(fftw_complex *array, int ncs_padded, char const *plotTitle, plotType type = NORMAL, plotStyle style = MAG, std::string storDir = "");	
		void plot(double *array, int ncs_padded, char const *plotTitle, std::string storDir = "");
		
		double mag(fftw_complex value){return sqrt(pow(value[0], 2) + pow(value[1], 2));};
};

#endif
