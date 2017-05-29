#include "signal_processor.hpp"

SignalProcessor::SignalProcessor(Experiment* exp)
{	
	timer.start();
	dopplerDataStart = 0;
	experiment = exp;
}

void SignalProcessor::allocateMemory(void)
{
	binDataset 		= (int16_t*)malloc(experiment->n_range_lines*2*experiment->ncs_range_line*sizeof(int16_t));
	refDataset      = (int16_t*)malloc(experiment->ncs_reference*2*sizeof(int16_t));
	
	lineBuffer 		= (fftw_complex*)malloc(experiment->ncs_padded*sizeof(fftw_complex));
	refBuffer 		= (fftw_complex*)malloc(experiment->ncs_padded*sizeof(fftw_complex));
	
	resultBuffer 	= (fftw_complex*)malloc(experiment->ncs_padded*sizeof(fftw_complex));
	dopplerBuffer   = (fftw_complex*)malloc(experiment->ncs_doppler_cpi*sizeof(fftw_complex));

	dopplerData     = (fftw_complex*)malloc(experiment->ncs_padded*experiment->ncs_doppler_cpi*sizeof(fftw_complex));
	
	matchedImageBuffer  = (double*)malloc(experiment->ncs_padded*sizeof(double));
	dopplerImageBuffer  = (double*)malloc(experiment->ncs_doppler_cpi*sizeof(double));	
	
	logger.write("Memory Allocated", timer);		
}


void SignalProcessor::createPlans(void)
{	
	//init fftw for multi-threading
	fftw_init_threads();
	fftw_plan_with_nthreads((*experiment).n_threads);
	
	//fftw plans
	rangePlan   = fftw_plan_dft_1d((*experiment).ncs_padded		, lineBuffer   , lineBuffer   , FFTW_FORWARD , FFTW_MEASURE );
	refPlan     = fftw_plan_dft_1d((*experiment).ncs_padded		, refBuffer    , refBuffer    , FFTW_FORWARD , FFTW_ESTIMATE);
	resultPlan  = fftw_plan_dft_1d((*experiment).ncs_padded		, resultBuffer , lineBuffer   , FFTW_BACKWARD, FFTW_MEASURE );
	dopplerPlan = fftw_plan_dft_1d((*experiment).ncs_doppler_cpi, dopplerBuffer, dopplerBuffer, FFTW_FORWARD , FFTW_MEASURE );	
	
	logger.write("FFT Plans Created", timer);			
}


void SignalProcessor::destroyPlans(void)
{
	fftw_destroy_plan(rangePlan);	
	fftw_destroy_plan(refPlan);	
	fftw_destroy_plan(resultPlan);	
	fftw_destroy_plan(dopplerPlan);
}


void SignalProcessor::fftRefData(void)
{			
	fftw_execute(refPlan);	
	//gnu_plot.gnuPlot(refBuffer, "reference waveform frequency domain", FFT_SHIFT);
}


void SignalProcessor::fftRangeData(void)
{	
	fftw_execute(rangePlan);	
	//gnu_plot.gnuPlot(rangeBuffer, "return waveform frequency domain", FFT_SHIFT);
}


void SignalProcessor::fftDopplerData(void)
{
	fftw_execute(dopplerPlan);
}


void SignalProcessor::ifftMatchedData(void)
{
	fftw_execute(resultPlan);
	//gnu_plot.gnuPlot(rangeBuffer, "matched result time domain", NORMAL, AMPLITUDE);
}


void SignalProcessor::popDopplerData(int rangeLine)
{
	if (rangeLine%experiment->update_rate == 0)
		dopplerDataStart = rangeLine;

	if ((rangeLine + 1 - dopplerDataStart) <= experiment->ncs_doppler_cpi)
	{
		for (int j = 0; j < experiment->ncs_padded; j++)
		{
			dopplerData[j*experiment->ncs_doppler_cpi + (rangeLine - dopplerDataStart)][0] = lineBuffer[j][0];
			dopplerData[j*experiment->ncs_doppler_cpi + (rangeLine - dopplerDataStart)][1] = lineBuffer[j][1];
		}
	}
}


void SignalProcessor::processDoppler(int rangeLine, OpenCVPlot &plot)
{
	popDopplerData(rangeLine); 
	
	if ((rangeLine - dopplerDataStart + 1) == experiment->ncs_doppler_cpi)  //check that dopplerData is full
	{
		for (int i = 0; i < experiment->ncs_padded; i++)		
		{
			popDopplerBuffer(i);	
			fftDopplerData();
			addToDopplerPlot(i, plot);
		}
		plot.plotDoppler();
	}
}


void SignalProcessor::popDopplerBuffer(int dopplerLine)
{
	for (int j = 0; j < experiment->ncs_doppler_cpi; j++)
	{	
		dopplerBuffer[j][0] = dopplerData[dopplerLine*experiment->ncs_doppler_cpi + j][0]*dopplerWindow.getSample(j); 
		dopplerBuffer[j][1] = dopplerData[dopplerLine*experiment->ncs_doppler_cpi + j][1]*dopplerWindow.getSample(j);
	}	
}


void SignalProcessor::addToDopplerPlot(int dopplerLine, OpenCVPlot &plot)
{
	float maxResult = 0.0f;	
	float result = 0.0f;
	double processed = 0;

	//find max result
	for (int i = 0; i < experiment->ncs_doppler_cpi; i++)
	{
		result = (sqrt(dopplerBuffer[i][0]*dopplerBuffer[i][0] + dopplerBuffer[i][1]*dopplerBuffer[i][1]));

		if (result > maxResult)
			maxResult = result;
	}

	for (int i = 0; i < experiment->ncs_doppler_cpi; i++)
	{
		processed = (((sqrt(dopplerBuffer[i][0]*dopplerBuffer[i][0] + dopplerBuffer[i][1]*dopplerBuffer[i][1]))/maxResult)*255);

		//perform fft shift
		if (i < (experiment->ncs_doppler_cpi/2 + 1))		
			dopplerImageBuffer[i + (experiment->ncs_doppler_cpi/2 - 1)] = processed;
		else
			dopplerImageBuffer[i - (experiment->ncs_doppler_cpi/2 + 1)] = processed;
	}	
	
	plot.addToDopplerPlot(dopplerLine, dopplerImageBuffer);
}


void SignalProcessor::complxConjRef(void)
{
	for (int i = 0; i < (experiment->ncs_padded); i++)
		refBuffer[i][1] = -1*refBuffer[i][1];
		
	logger.write("Complex Conjugate Reference", timer);	
}


void SignalProcessor::complxMulti(void)
{
	for (int j = 0; j < (experiment->ncs_padded); j++)
	{			
		resultBuffer[j][0] = (lineBuffer[j][0]*refBuffer[j][0] - lineBuffer[j][1]*refBuffer[j][1]);
		resultBuffer[j][1] = (lineBuffer[j][0]*refBuffer[j][1] + lineBuffer[j][1]*refBuffer[j][0]);
	}
	//plot.gnuPlot(hilbertBuffer, "matched result frequency domain", FFT_SHIFT);
}


//process data extracted from the bin file into the complex lineBuffer line by line.
void SignalProcessor::popRangeBuffer(int rangeLine)
{
	int start = rangeLine*2*experiment->ncs_range_line;
		
	//populate complex range data and remove offset	
	for (int i = 0; i < experiment->ncs_padded; i++)
	{
		if (i < experiment->ncs_range_line)
		{	
			lineBuffer[i][0] = binDataset[i*2 + start    ]; //*rangeWindow.getSample(i);     //real component    
			lineBuffer[i][1] = binDataset[i*2 + start + 1]; //*rangeWindow.getSample(i);     //complex component
		}
		else
		{
			lineBuffer[i][0] = 0;
			lineBuffer[i][1] = 0; 
		}
	}	
	//plot.gnuPlot(rangeBuffer, "return waveform time domain", NORMAL, IQ);
}


void SignalProcessor::freeMemory(void)
{
	//fftw_free(fftRangeBuffer);
	fftw_free(refBuffer);
	fftw_free(lineBuffer);
	fftw_free(resultBuffer);
	fftw_free(dopplerBuffer);
	fftw_free(dopplerData);
	
	fftw_cleanup_threads();
	
	free(binDataset);
	free(refDataset);
	
	destroyPlans();

	logger.write("Memory Free \n", timer);
}


//pulse compressed range lines are normalized and converted to 8 bit for image generation. 
void SignalProcessor::addToWaterPlot(int rangeLine, OpenCVPlot &plot)
{
	for (int j = 0; j < experiment->ncs_padded; j++)
	{
		matchedImageBuffer[j] = mag(lineBuffer[j]);
	}	
	
	plot.addToWaterPlot(rangeLine, matchedImageBuffer);
	//plot.gnuPlot(matchedImageBuffer, "matched image buffer");
}


void SignalProcessor::loadBinaryDataset(void)
{
	//declare a file pointer
	FILE *binFile;	
	
	//assign pointer to file location
	binFile = fopen(experiment->dataset_filename, "rb");			
	
	//check that file exists in the location specified
	if (binFile != NULL)
	{
		//read from binary file into buffer
		fread(binDataset, sizeof(int16_t), experiment->n_range_lines*2*experiment->ncs_range_line, binFile);

		//close binary files
		fclose(binFile);	

		logger.write("Binary Dataset Loaded", timer);
	}
	//file does not exist in the specified location
	else
	{
		logger.write("Binary dataset could not be found in the specified location.");
		exit(0);
	}
}


void SignalProcessor::loadReferenceWaveform(void)
{
	//declare a file pointer
	FILE *refFile;	
	
	//assign pointer to file location
	refFile = fopen(experiment->reference_filename, "rb");
	
	//check that file exists in the location specified
	if (refFile != NULL)
	{
		//fseek(refFile, 2*COMPLEX_START_SAMPLE_IN_REFERENCE_WAVEFORM*sizeof(int16_t), SEEK_SET);
		//read from binary file into buffer
		fread(refDataset, sizeof(int16_t), experiment->ncs_reference*2, refFile);

		//close binary files
		fclose(refFile);	

		logger.write("Reference Dataset Loaded", timer);
	}
	//file does not exist in the specified location
	else
	{
		logger.write("Reference dataset could not be found in the specified location.");
		exit(0);
	} 
 

	//populate complex reference data and remove offset	
	for (int i = 0; i < experiment->ncs_padded; i++)
	{
		if (i < experiment->ncs_reference)
		{	
			refBuffer[i][0] = refDataset[i*2    ]*rangeWindow.getSample(i);     //real component    
			refBuffer[i][1] = refDataset[i*2 + 1]*rangeWindow.getSample(i);  	//complex component
		}
		else
		{
			refBuffer[i][0] = 0;
			refBuffer[i][1] = 0; 
		}
	}	

	logger.write("Reference Data Loaded", timer);	
	//gnu_plot.gnuPlot(refBuffer, "reference waveform time domain", NORMAL, IQ);		
}


void SignalProcessor::getExperimentParameters(void)
{
	//check that filepath is good
	std::ifstream check(EXP_FILE);

	if (!check.good()) 
	{
		printf("Please check location of experiment.ini and try again.\n");
		exit(EXIT_FAILURE);
	}
	
	CSimpleIniA ini;
	ini.LoadFile(EXP_FILE);	
	
	experiment->dataset_filename = (char *)ini.GetValue("dataset", "data_filename");
	experiment->ncs_range_line = atoi(ini.GetValue("dataset", "n_cmplx_samples_range_line"));
	experiment->n_range_lines = atoi(ini.GetValue("dataset", "n_range_lines"));
	
	experiment->reference_filename = (char *)ini.GetValue("dataset", "ref_filename");
	experiment->ncs_reference = atoi(ini.GetValue("dataset", "n_cmplx_samples_ref"));
	
	experiment->ncs_padded = atoi(ini.GetValue("dataset", "n_cmplx_samples_padded"));
	
	experiment->n_threads = atoi(ini.GetValue("processing", "n_threads"));	
	experiment->ncs_doppler_cpi = atoi(ini.GetValue("processing", "doppler_cpi"));		
		
	std::string doppler_flag = ini.GetValue("processing", "doppler_enabled");	
		
	if ((doppler_flag == "1") || (doppler_flag == "true") || (doppler_flag == "TRUE") || (doppler_flag == "True"))
		experiment->is_doppler = true;
	else
		experiment->is_doppler = false;
		
	std::string debug_flag = ini.GetValue("config", "debug_mode");	
		
	if ((debug_flag == "1") || (debug_flag == "true") || (debug_flag == "TRUE") || (debug_flag == "True"))
		experiment->is_debug = true;
	else
		experiment->is_debug = false;
		
	experiment->update_rate = atoi(ini.GetValue("visualisation", "update_rate"));
	experiment->cm_doppler = atoi(ini.GetValue("visualisation", "doppler_colour_map"));
	experiment->cm_rti = atoi(ini.GetValue("visualisation", "rti_colour_map"));
	experiment->hist_equal = atoi(ini.GetValue("visualisation", "histogram_equalization"));
	experiment->slow = atoi(ini.GetValue("visualisation", "slow"));
	experiment->threshold = atoi(ini.GetValue("visualisation", "threshold"));
	
	rangeWindow.init(HAMMING, experiment->ncs_reference);		
	dopplerWindow.init(HAMMING, experiment->ncs_doppler_cpi);	
}






