#include "plotting.hpp"

GNUPlot::GNUPlot(Experiment* exp)
{
	experiment = exp;
}

void GNUPlot::gnuPlot(uint8_t *array, char const *plotTitle, Experiment* exp)
{
	if (exp->is_debug)
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
		fputs(writeFileName, pipe_gp);
		
		fputs("plot '-' using 1:2 with lines notitle\n", pipe_gp);
		for (int i = 0; i < exp->ncs_padded; i++) 
		{
			//float magnitude = 10*log(sqrt(pow(array[i][0], 2) + pow(array[i][1], 2)));							
			fprintf(pipe_gp, "%i %i\n", i, array[i]);
		}
		
		fputs("e\n", pipe_gp);
		pclose(pipe_gp);	
		
		char writeMessage[100];
		strcpy(writeMessage, "Generated Plot: ");
		strcat(writeMessage, plotTitle);	
		
		logger.write(writeMessage);
	}
}



void GNUPlot::gnuPlot(fftw_complex *array, char const *plotTitle, Experiment* exp, plotType type, plotStyle style)
{
	if (exp->is_debug)
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
						for (int i = 0; i < exp->ncs_padded; i++) 
						{
							float magnitude = 10*log(sqrt(pow(array[i][0], 2) + pow(array[i][1], 2)));							
							fprintf(pipe_gp, "%i %f\n", i, magnitude);
						}
						break;
					}
					
					case IQ:
					{
						fputs("plot '-' title 'I Samples' with lines, '-' title 'Q Samples' with lines\n", pipe_gp);
						
						for (int i = 0; i < exp->ncs_padded; i++) 
						{
							fprintf(pipe_gp, "%i %f\n", i, array[i][0]);
						}
						fflush(pipe_gp);
						fprintf(pipe_gp, "e\n");						

						for (int i = 0; i < exp->ncs_padded; i++) 
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
				for (int i = 0; i < exp->ncs_padded; i++) 
				{
					if (i < (exp->ncs_padded/2 + 1)) 
					{
						fprintf(pipe_gp, "%i %f\n", i, (abs(sqrt(array[i + (exp->ncs_padded/2 - 1)][0]*array[i + (exp->ncs_padded/2 - 1)][0] +
																 array[i + (exp->ncs_padded/2 - 1)][1]*array[i + (exp->ncs_padded/2 - 1)][1]))));
					}
					else
					{
						fprintf(pipe_gp, "%i %f\n", i, (abs(sqrt(array[i - (exp->ncs_padded/2 + 1)][0]*array[i - (exp->ncs_padded/2 + 1)][0] + 
																 array[i - (exp->ncs_padded/2 + 1)][1]*array[i - (exp->ncs_padded/2 + 1)][1]))));
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

OpenCVPlot::OpenCVPlot(Experiment* exp)
{
	experiment = exp;
}


void OpenCVPlot::initOpenCV(void)
{	
	rtCMapSldr = experiment->cm_rti;
	rdCMapSldr = experiment->cm_doppler;
	thrsSldr = experiment->threshold;
	slowSldr = experiment->slow;
	histSldr = experiment->hist_equal;
	avrgSldr = experiment->n_plot_average;
	
	rtSize = cv::Size(500, 500);
	rdSize = cv::Size(250, 500);
	
	cv::namedWindow("RTI Plot");
	cv::moveWindow("RTI Plot", 0, 0);	
	rtImage = cv::Mat::ones(experiment->n_range_lines, experiment->ncs_padded, CV_64F);
	
	cv::namedWindow("Control Window", cv::WINDOW_NORMAL);	
	cv::moveWindow("Control Window", rtSize.width + rdSize.width, 0); 
	//cv::resizeWindow("Control Window", 300, 100);
	
	if (experiment->is_doppler)
	{
		cv::namedWindow("Doppler Plot");
		cv::moveWindow("Doppler Plot", rtSize.width, 0); 
		cv::createTrackbar( "Doppler Colour Map", "Control Window", &rdCMapSldr, cMapMax);
		rdImage = cv::Mat::ones(experiment->n_range_lines, experiment->ncs_doppler_cpi, CV_64F);
		
		avrgMax = (experiment->n_range_lines)/(experiment->update_rate);
		cv::createTrackbar( "Doppler Averaging", "Control Window", &avrgSldr, avrgMax);
		
		rdIndex = 0;			
		rdImageAvg = cv::Mat(rdSize.height, rdSize.width, CV_64F, cv::Scalar::all(0));	
	}	

	cv::createTrackbar( "Threshold Value", "Control Window", &thrsSldr, thrsMax);	
	cv::createTrackbar( "RTI Colour Map", "Control Window", &rtCMapSldr, cMapMax);
	cv::createTrackbar( "Slow Processing", "Control Window", &slowSldr, slowMax);
	cv::createTrackbar( "Histogram Equalisation", "Control Window", &histSldr, histMax);
}


void OpenCVPlot::addRTI(int rangeLine, double  *imageValues)
{
	cv::Mat matchedRow = cv::Mat(1, experiment->ncs_padded, CV_64F, imageValues);	
	cv::log(matchedRow, matchedRow);
	matchedRow.copyTo(rtImage(cv::Rect(0, rangeLine, matchedRow.cols, matchedRow.rows)));
	
	if (((rangeLine%(experiment->update_rate - 1) == 0) || rangeLine == (experiment->n_range_lines - 1)) && rangeLine != 0)
	{
		experiment->mutex.lock();
		plotRTI();
		experiment->mutex.unlock();
	}
}


void OpenCVPlot::addRD(int dopplerLine, double *imageValues)
{
	cv::Mat row = cv::Mat(1, experiment->ncs_doppler_cpi, CV_64F, imageValues);
	rdImage.push_back(row);
}


void OpenCVPlot::plotRTI(void)
{
	// uses bilinear interpolation to reduce number of pixels (decimation)
	cv::resize(rtImage, rtImageResize, rtSize);		
	
	cv::normalize(rtImageResize, rtImageResize, 0.0, 1.0, cv::NORM_MINMAX);

	rtImageResize.convertTo(rtImage8bit, CV_8U, 255);	

	if (histSldr)
	{
		cv::equalizeHist(rtImage8bit, rtImage8bit);
	}
	
	cv::threshold(rtImage8bit, rtImage8bit, thrsSldr, thrsMax, cv::THRESH_TOZERO);	
	cv::applyColorMap(rtImage8bit, rtImage8bit, rtCMapSldr);	
	cv::transpose(rtImage8bit, rtImage8bit);
	cv::flip(rtImage8bit, rtImage8bit, 0);		
	
	cv::imshow("RTI Plot", rtImage8bit);

	cv::waitKey(1 + slowSldr);	
}


void OpenCVPlot::plotRD(void)
{
	int summable_plots = 0;
	
	//add new dummy doppler plot to the vector
	rdVector.push_back(cv::Mat::ones(1, 1, CV_64F));
	
	cv::resize(rdImage, rdVector[rdIndex], rdSize);		
	rdImage.release();
	
	//init average as latest plot
	rdImageAvg = rdVector[rdIndex];
	
	//clear the average
	rdImageAvg = cv::Mat(rdSize.height, rdSize.width, CV_64F, cv::Scalar::all(0));
	
	//determine the number of plots available for averaging
	if (avrgSldr >= rdIndex)
	{
		summable_plots = rdIndex;
	}
	else
	{
		summable_plots = avrgSldr;
	}
	
	//sum all appropriate plots
	for (int i = (rdIndex - 1); i > (rdIndex - summable_plots); i--)
	{
		cv::add(rdImageAvg, rdVector[i], rdImageAvg);
	}
	
	cv::log(rdImageAvg, rdImageAvg);
	
	cv::normalize(rdImageAvg, rdImageAvg, 0.0, 1.0, cv::NORM_MINMAX);
	
	rdImageAvg.convertTo(rdImage8bit, CV_8U, 255);
	
	cv::threshold(rdImage8bit, rdImage8bit, thrsSldr, thrsMax, cv::THRESH_TOZERO);	
	
	if (histSldr)
	{
		cv::equalizeHist(rdImage8bit, rdImage8bit);
	}
	
	cv::applyColorMap(rdImage8bit, rdImage8bit, rdCMapSldr);
	
	//vertical flip through x-axis
	cv::flip(rdImage8bit, rdImage8bit, 0);
	
	cv::imshow("Doppler Plot", rdImage8bit);
	
	//increment the doppler plot index
	rdIndex++;
}

void OpenCVPlot::savePlots(void)
{
	std::string doppler_path = experiment->save_path + "/Range-Doppler.jpg";
	std::string water_path = experiment->save_path + "/Range-Time-Intensity.jpg";
	
	cv::imwrite(doppler_path.c_str(), rdImage8bit); 
	cv::imwrite(water_path.c_str(), rtImage8bit);	
	
	for (int i = 0; i < experiment->n_plot_average; i++)
	{
		rdVector[i].release();	
	}
	
	rtImage8bit.release();
}

