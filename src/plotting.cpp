#include "plotting.hpp"

double LOG10 = log(10);

void GNUPlot::plot(double *array, int ncs_padded, char const *plotTitle, std::string storDir)
{
	FILE *pipe_gp = popen("gnuplot", "w");	

	char writeTitle[100];
	strcpy(writeTitle, "set title '");
	strcat(writeTitle, plotTitle);
	strcat(writeTitle, "'\n");
	
	char writeFileName[200];
	strcpy(writeFileName, "set output '");
	strcat(writeFileName, storDir.c_str());
	strcat(writeFileName, "/");
	strcat(writeFileName, plotTitle);
	strcat(writeFileName, ".eps'\n");	

	fputs("set terminal postscript eps enhanced color font 'Helvetica,20' linewidth 0.5\n", pipe_gp);
	fputs(writeTitle, pipe_gp);	
	fputs(writeFileName, pipe_gp);
	
	fputs("plot '-' using 1:2 with lines notitle\n", pipe_gp);
	for (int i = 0; i < ncs_padded; i++) 
	{
		//float magnitude = 10*log(sqrt(pow(array[i][0], 2) + pow(array[i][1], 2)));							
		fprintf(pipe_gp, "%i %f\n", i, array[i]);
	}
	
	fputs("e\n", pipe_gp);
	pclose(pipe_gp);	
	
	char writeMessage[100];
	strcpy(writeMessage, "Generated Plot: ");
	strcat(writeMessage, plotTitle);	
	
	logger.write(writeMessage);
}



void GNUPlot::plot(fftw_complex *array, int ncs_padded, char const *plotTitle, plotType type, plotStyle style, std::string storDir)
{	
	FILE *pipe_gp = popen("gnuplot", "w");	
	
	char writeTitle[100];
	strcpy(writeTitle, "set title '");
	strcat(writeTitle, plotTitle);
	strcat(writeTitle, "'\n");
	
	char writeFileName[200];
	strcpy(writeFileName, "set output '");
	strcat(writeFileName, storDir.c_str());
	strcat(writeFileName, "/");
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
				case MAG: 
				{
					fputs("plot '-' using 1:2 with lines notitle\n", pipe_gp);
					for (int i = 0; i < ncs_padded; i++) 
					{
						fprintf(pipe_gp, "%i %f\n", i, mag(array[i]));
					}
					break;
				}
				
				case LOG: 
				{
					fputs("plot '-' using 1:2 with lines notitle\n", pipe_gp);
					for (int i = 0; i < ncs_padded; i++) 
					{
						fprintf(pipe_gp, "%i %f\n", i, log10(mag(array[i])));
					}
					break;
				}
				
				case IQ:
				{
					fputs("plot '-' title 'I Samples' with lines, '-' title 'Q Samples' with lines\n", pipe_gp);
					
					for (int i = 0; i < ncs_padded; i++) 
					{
						fprintf(pipe_gp, "%i %f\n", i, array[i][0]);
					}
					fflush(pipe_gp);
					fprintf(pipe_gp, "e\n");						

					for (int i = 0; i < ncs_padded; i++) 
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
			for (int i = 0; i < ncs_padded; i++) 
			{
				if (i < (ncs_padded/2 + 1)) 
				{
					fprintf(pipe_gp, "%i %f\n", i, (double_t)(abs(sqrt(array[i + (ncs_padded/2 - 1)][0]*array[i + (ncs_padded/2 - 1)][0] +
															 array[i + (ncs_padded/2 - 1)][1]*array[i + (ncs_padded/2 - 1)][1]))));
				}
				else
				{
					fprintf(pipe_gp, "%i %f\n", i, (double_t)(abs(sqrt(array[i - (ncs_padded/2 + 1)][0]*array[i - (ncs_padded/2 + 1)][0] + 
															 array[i - (ncs_padded/2 + 1)][1]*array[i - (ncs_padded/2 + 1)][1]))));
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

OpenCVPlot::OpenCVPlot(Experiment* exp)
{
	experiment = exp;
}


void OpenCVPlot::initOpenCV(void)
{	
	cmapSldr = experiment->cm;
	thrsSldr = experiment->threshold;
	slowSldr = experiment->slow;
	histSldr = experiment->hist_equal;
	
	rtSize = cv::Size(500, 500);
	rdSize = cv::Size(200, 500);
	spSize = cv::Size(200, 200);
	
	if (experiment->is_visualisation)
	{
		cv::namedWindow("Range-Time-Intensity", cv::WINDOW_AUTOSIZE);
		cv::moveWindow("Range-Time-Intensity", 0, 0);	
		
		cv::namedWindow("Control", cv::WINDOW_AUTOSIZE);
		cv::moveWindow("Control", (rtSize.width + 2), 0);		
		cv::createTrackbar( "Colour Map", "Control", &cmapSldr, cMapMax);
	}
	
	//init the range-time image with n_range_lines rows
	//the number of columns is equal to ncs_padded - ncs_reference off the front and ncs_reference/2 off the back of the image
	//all values are init to the blanking threshold 
	rtImage = cv::Mat(experiment->n_range_lines, experiment->ncs_padded - (experiment->ncs_blank + experiment->ncs_blank/2), CV_64F, cv::Scalar::all(experiment->blanking_threshold));	
	
	if (experiment->is_doppler)
	{
		rdImage = cv::Mat::ones(experiment->n_range_lines, experiment->ncs_doppler_cpi*experiment->doppler_padding_factor, CV_64F);
		
		avrgMax = (experiment->n_range_lines)/(experiment->update_rate);
		avrgSldr = avrgMax;	//set slider to max by default
		
		if (experiment->is_visualisation)
		{
			cv::namedWindow("Range-Doppler");
			cv::moveWindow("Range-Doppler", (rtSize.width + 2), 0); 
				
			cv::namedWindow("Spectrogram", cv::WINDOW_AUTOSIZE);
			cv::moveWindow("Spectrogram", (rtSize.width + rdSize.width + 2*2), 0);	
			
			cv::moveWindow("Control", (rtSize.width + rdSize.width + spSize.height + 3*2), 0); 		
		
			cv::createTrackbar( "RD Averaging", "Control", &avrgSldr, avrgMax);
		}
		
		rdIndex = 0;			
		rdImageAvg = cv::Mat(rdSize.height, rdSize.width, CV_64F, cv::Scalar::all(0));	
	}	

	if (experiment->is_visualisation)
	{
		cv::createTrackbar( "Threshold Value", "Control", &thrsSldr, thrsMax);	
		cv::createTrackbar( "Slow [ms]", "Control", &slowSldr, slowMax);
		cv::createTrackbar( "Histogram Eqn", "Control", &histSldr, histMax);
	}
}


void OpenCVPlot::addRTI(int rangeLine, double  *imageValues)
{
	//similarly to the init, the number of samples in the row equals
	//ncs_padded - ncs_reference off the front and ncs_reference/2 off the back of the image
	cv::Mat matchedRow = cv::Mat(1, experiment->ncs_padded - (experiment->ncs_blank + experiment->ncs_blank/2), CV_64F, &imageValues[experiment->ncs_blank]);	
	
	matchedRow.copyTo(rtImage(cv::Rect(0, rangeLine, matchedRow.cols, matchedRow.rows)));
	
	if (experiment->is_visualisation)
	{
		if (((rangeLine%(experiment->update_rate - 1) == 0) || (rangeLine == (experiment->n_range_lines - 1))) && rangeLine != 0) 
		{
			experiment->mutex.lock();
			plotRTI();
			experiment->mutex.unlock();
		}
	}
}


void OpenCVPlot::addRD(double *imageValues)
{
	cv::Mat row = cv::Mat(1, experiment->ncs_doppler_cpi*experiment->doppler_padding_factor, CV_64F, imageValues);	
	rdImage.push_back(row);
}

void OpenCVPlot::addSP(double *imageValues)
{
	cv::Mat row = cv::Mat(1, experiment->ncs_doppler_cpi*experiment->doppler_padding_factor, CV_64F, imageValues);
	spImage.push_back(row);
}


void OpenCVPlot::plotRTI(void)
{
	//use bilinear interpolation to reduce number of pixels (decimation)
	cv::resize(rtImage, rtImageResize, rtSize);		
	
	cv::normalize(rtImageResize, rtImageResize, 0.0, 1.0, cv::NORM_MINMAX);
	
	//gPlot.plot(matchedImageBuffer, experiment->ncs_padded, "Matched Image Buffer Pulse #1");	

	rtImageResize.convertTo(rtImage8bit, CV_8U, 255);	

	if (histSldr)
	{
		cv::equalizeHist(rtImage8bit, rtImage8bit);
	}
	
	cv::threshold(rtImage8bit, rtImage8bit, thrsSldr, thrsMax, cv::THRESH_TOZERO);
	
	cv::applyColorMap(rtImage8bit, rtImage8bit, cmapSldr);	
	cv::transpose(rtImage8bit, rtImage8bit);
	cv::flip(rtImage8bit, rtImage8bit, 0);		
	
	if (experiment->is_visualisation)
	{
		std::string rtCMapPath = cMapRoot + boost::to_string(cmapSldr) + ".jpg";
		cv::Mat rtCMap = cv::imread(rtCMapPath);
		cv::imshow("Control", rtCMap);
		
		cv::imshow("Range-Time-Intensity", rtImage8bit);
		cv::waitKey(1 + slowSldr);
	}		
}


void OpenCVPlot::plotSP(void)
{
	//use bilinear interpolation to reduce number of pixels (decimation)
	cv::resize(spImage, spImageResize, spSize);		
	
	cv::log(spImageResize, spImageResize);
	spImageResize = spImageResize/LOG10;
	
	cv::normalize(spImageResize, spImageResize, 0.0, 1.0, cv::NORM_MINMAX);
	
	spImageResize.convertTo(spImage8bit, CV_8U, 255);
	
	if (histSldr)
	{
		cv::equalizeHist(spImage8bit, spImage8bit);
	}
	
	cv::threshold(spImage8bit, spImage8bit, thrsSldr, thrsMax, cv::THRESH_TOZERO);	
	
	cv::applyColorMap(spImage8bit, spImage8bit, cmapSldr);
	
	//vertical flip through x-axis
	cv::flip(spImage8bit, spImage8bit, 0);
	
	cv::transpose(spImage8bit, spImage8bit);
	
	if (experiment->is_visualisation)
	{
		cv::imshow("Spectrogram", spImage8bit);
	}
}


void OpenCVPlot::plotRD(void)
{
	//legacy code from the stand alone version
	//cobalt version will always compute average over all Doppler plots.
	
	int summable_plots;
	
	//determine the number of plots available for averaging
	if (avrgSldr >= rdIndex)
	{
		summable_plots = rdIndex;
	}
	else
	{
		summable_plots = avrgSldr;
	}
	
	//add new dummy doppler plot to the vector
	rdVector.push_back(cv::Mat::ones(1, 1, CV_64F));
	
	//use bilinear interpolation to reduce number of pixels (decimation)
	cv::resize(rdImage, rdVector[rdIndex], rdSize);		
	
	//original range-Doppler plot now stored in vector
	rdImage.release();
	
	//clear the average
	rdImageAvg = cv::Mat(rdSize.height, rdSize.width, CV_64F, cv::Scalar::all(0));
	
	//init average as the current range-Doppler plot
	cv::add(rdImageAvg, rdVector[rdIndex], rdImageAvg);	
	
	//sum all appropriate plots
	for (int i = (rdIndex - 1); i > (rdIndex - summable_plots); i--)
	{
		cv::add(rdImageAvg, rdVector[i], rdImageAvg);
	}
	
	cv::log(rdImageAvg, rdImageAvg);
	rdImageAvg = rdImageAvg/LOG10;
	
	cv::normalize(rdImageAvg, rdImageAvg, 0.0, 1.0, cv::NORM_MINMAX);
	
	rdImageAvg.convertTo(rdImage8bit, CV_8U, 255);
	
	cv::threshold(rdImage8bit, rdImage8bit, thrsSldr, thrsMax, cv::THRESH_TOZERO);	
	
	if (histSldr)
	{
		cv::equalizeHist(rdImage8bit, rdImage8bit);
	}
	
	cv::applyColorMap(rdImage8bit, rdImage8bit, cmapSldr);	
	
	//vertical flip through x-axis
	cv::flip(rdImage8bit, rdImage8bit, 0);
	
	if (experiment->is_visualisation)
	{
		cv::imshow("Range-Doppler", rdImage8bit);
	}
	
	//increment the doppler plot index
	rdIndex++;
}

void OpenCVPlot::savePlots(void)
{
	
	rtPath = experiment->save_path + "/Range-Time-Intensity.jpg";	
	
	cv::imwrite(rtPath.c_str(), rtImage8bit);	
	
	rtImage8bit.release();
	
	if (experiment->is_doppler)
	{
		rdPath = experiment->save_path + "/Range-Doppler.jpg";
		cv::imwrite(rdPath.c_str(), rdImage8bit); 
		
		for (int i = 0; i < experiment->n_plot_average; i++)
		{
			rdVector[i].release();	
		}
		
		spPath = experiment->save_path + "/Spectrogram.jpg";
		cv::imwrite(spPath.c_str(), spImage8bit); 
	}
}

