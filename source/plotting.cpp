#include "plotting.hpp"

//globals
cv::Mat waterImage, doppImage;
int waterfallColourMapSlider = 0;
int dopplerColourMapSlider = 0;
int	thresholdSlider = 0;

const int thresholdMax = 255;
const int colourMapMax = 11;


Plot::Plot(void)
{

}


void Plot::gnuPlot(fftw_complex *array, char const *plotTitle, plotType type, plotStyle style)
{
	if (DEBUG_MODE)
	{
		FILE *pipe_gp = popen("gnuplot", "w");	
	
		char writeTitle[100];
		strcpy(writeTitle, "set title '");
		strcat(writeTitle, plotTitle);
		strcat(writeTitle, "'\n");
		
		char writeFileName[100];
		strcpy(writeFileName, "set output '");
		strcat(writeFileName, plotTitle);
		strcat(writeFileName, ".eps'\n");	

		fputs("set terminal postscript eps enhanced color font 'Helvetica,20' linewidth 0.5\n", pipe_gp);
		fputs(writeTitle, pipe_gp);	
		//fputs("set yrange [-1:1] \n", pipe_gp);
		//fputs("set xrange[500:600] \n", pipe_gp);
		//fputs("set bmargin at screen 0.13 \n", pipe_gp);
		//fputs("set lmargin at screen 0.13 \n", pipe_gp);
		//fputs("set xtics ('0' 0, '10' 205, '20' 410, '30' 615, '40' 820, '50' 1024) \n", pipe_gp);
		//fputs("set xtics ('-500' 0, '-376' 32, '-252' 64, '-128' 96, '0' 128, '128' 160, '252' 192, '376' 224, '500' 255) \n", pipe_gp);
		//fputs("set xlabel 'Frequency [Hz]' \n", pipe_gp);
		//fputs("set xlabel 'Sample Number' \n", pipe_gp);
		//fputs("set ylabel 'Magnitude' \n", pipe_gp);
		fputs(writeFileName, pipe_gp);
		
		switch(type)
		{
			case NORMAL:
			{
				switch (style)
				{	
					case AMPLITUDE: 
					{
						fputs("plot '-' using 1:2 with lines notitle\n", pipe_gp);
						for (int i = 0; i < COMPLEX_SAMPLES_AFTER_ZERO_PADDING; i++) 
						{
							float magnitude = sqrt(pow(array[i][0], 2) + pow(array[i][1], 2));							
							fprintf(pipe_gp, "%i %f\n", i, magnitude);
						}
						break;
					}
					
					case IQ:
					{
						fputs("plot '-' title 'I Samples' with lines, '-' title 'Q Samples' with lines\n", pipe_gp);
						
						for (int i = 0; i < COMPLEX_SAMPLES_AFTER_ZERO_PADDING; i++) 
						{
							fprintf(pipe_gp, "%i %f\n", i, array[i][0]);
						}
						fflush(pipe_gp);
						fprintf(pipe_gp, "e\n");						

						for (int i = 0; i < COMPLEX_SAMPLES_AFTER_ZERO_PADDING; i++) 
						{
							fprintf(pipe_gp, "%i %f\n", i, array[i][1]);
						}
						break;
					}										
				}	
				break;
			}
			
			case FFT_SHIFT:
			{
				fputs("plot '-' using 1:2 with lines notitle\n", pipe_gp);
				for (int i = 0; i < COMPLEX_SAMPLES_AFTER_ZERO_PADDING; i++) 
				{
					if (i < (COMPLEX_SAMPLES_AFTER_ZERO_PADDING/2 + 1))
					{
						fprintf(pipe_gp, "%i %f\n", i, (abs(sqrt(array[i + (COMPLEX_SAMPLES_AFTER_ZERO_PADDING/2 - 1)][0]*array[i + (COMPLEX_SAMPLES_AFTER_ZERO_PADDING/2 - 1)][0] +
																 array[i + (COMPLEX_SAMPLES_AFTER_ZERO_PADDING/2 - 1)][1]*array[i + (COMPLEX_SAMPLES_AFTER_ZERO_PADDING/2 - 1)][1]))));
					}
					else
					{
						fprintf(pipe_gp, "%i %f\n", i, (abs(sqrt(array[i - (COMPLEX_SAMPLES_AFTER_ZERO_PADDING/2 + 1)][0]*array[i - (COMPLEX_SAMPLES_AFTER_ZERO_PADDING/2 + 1)][0] + 
																 array[i - (COMPLEX_SAMPLES_AFTER_ZERO_PADDING/2 + 1)][1]*array[i - (COMPLEX_SAMPLES_AFTER_ZERO_PADDING/2 + 1)][1]))));
					}
				}
				break;	
			}
		}	

		fputs("e\n", pipe_gp);
		pclose(pipe_gp);	
		
		char writeMessage[100];
		strcpy(writeMessage, "Generated Plot: ");
		strcat(writeMessage, plotTitle);	
		
		logger.write(writeMessage);
	}	
}


void initOpenCV(void)
{	
	cv::namedWindow("Waterfall Plot");
	cv::namedWindow("Doppler Plot");
	cv::namedWindow("Control Window", cv::WINDOW_NORMAL);

	cv::moveWindow("Waterfall Plot", 0, 0);					//trackbar is 54 units in height
	cv::moveWindow("Doppler Plot", 500, 0); 
	cv::moveWindow("Control Window", 100, 0); 

	cv::resizeWindow("Control Window", 300, 3*54);

	cv::createTrackbar( "Threshold Value", "Control Window", &thresholdSlider, thresholdMax);
	cv::createTrackbar( "Doppler Colour Map", "Control Window", &dopplerColourMapSlider, colourMapMax);
	cv::createTrackbar( "Waterfall Colour Map", "Control Window", &waterfallColourMapSlider, colourMapMax);
}

void updateWaterfall(int rangeLine, uint8_t  *imageValues)
{
	cv::Mat row = cv::Mat(1, COMPLEX_SAMPLES_AFTER_ZERO_PADDING/2 + 1, CV_8U, imageValues);
	waterImage.push_back(row);
			
	if (((rangeLine%(UPDATE_LINE - 1) == 0) || rangeLine == (NUMBER_OF_RANGE_LINES - 1)) && rangeLine != 0)
		plotWaterfall();
}

void updateDoppler(uint8_t  *imageValues)
{
	cv::Mat row = cv::Mat(1, RANGE_LINES_PER_DOPPLER_CPI, CV_8U, imageValues);
	doppImage.push_back(row);
}

void plotWaterfall(void)
{
		cv::Mat resizedImage;	

		cv::Size size(500, 500);		
		cv::resize(waterImage, resizedImage, size);	
		
		cv::equalizeHist(resizedImage, resizedImage);	
		cv::threshold(resizedImage, resizedImage, thresholdSlider, thresholdMax, 3);
		cv::applyColorMap(resizedImage, resizedImage, waterfallColourMapSlider);	
		
		cv::transpose(resizedImage, resizedImage);
		cv::flip(resizedImage, resizedImage, 0);

		cv::imshow("Waterfall Plot", resizedImage);
		cv::imwrite("waterfall_plot.png", resizedImage);
		cv::waitKey(1);	
		resizedImage.release();
		//waterImage.release();
}

void plotDoppler(void)
{
		cv::Mat resizedImage;	

		cv::Size size(250, 500);		
		cv::resize(doppImage, resizedImage, size);	
		
		doppImage.release();
		
		cv::equalizeHist(resizedImage, resizedImage);
		cv::threshold(resizedImage, resizedImage, thresholdSlider, thresholdMax, 3);		
		cv::applyColorMap(resizedImage, resizedImage, dopplerColourMapSlider);
		
		cv::flip(resizedImage, resizedImage, 0);
	
		cv::imshow("Doppler Plot", resizedImage);
		//printf("Doppler Plot:\t\t\tOK\t%fs\n", getTime());
		cv::imwrite("doppler_plot.png", resizedImage);
		cv::waitKey(1);
		resizedImage.release();	
		doppImage.release();
}

