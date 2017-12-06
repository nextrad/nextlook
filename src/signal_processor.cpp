#include "signal_processor.hpp"

SignalProcessor::SignalProcessor(Experiment* exp)
{	
	timer.start();
	dopplerDataStart = 0;
	experiment = exp;
	
	experiment->dataset_filename = (char*)"-1";
	experiment->reference_filename = (char*)"-1";	
	experiment->n_range_lines = -1;		
	experiment->ncs_range_line = -1;	
	experiment->ncs_reference = -1;	
	experiment->ncs_padded = -1;			
	experiment->ncs_doppler_cpi = -1;		
	experiment->doppler_padding_factor = -1;		
	experiment->specro_range_bin = -1;	
	experiment->n_threads = -1;		
	experiment->pulse_blanking = -1;
	experiment->blanking_threshold = -1;
	experiment->dynamic_range = -1;
}

void SignalProcessor::allocateMemory(void)
{
	binDataset 		= (int16_t*)malloc(experiment->n_range_lines*2*experiment->ncs_range_line*sizeof(int16_t));
	refDataset      = (int16_t*)malloc(experiment->ncs_reference*2*sizeof(int16_t));
	
	lineBuffer 		= (fftw_complex*)malloc(experiment->n_threads*experiment->ncs_padded*sizeof(fftw_complex));
	refBuffer 		= (fftw_complex*)malloc(experiment->ncs_padded*sizeof(fftw_complex));
	
	resultBuffer 	= (fftw_complex*)malloc(experiment->n_threads*experiment->ncs_padded*sizeof(fftw_complex));
	dopplerBuffer   = (fftw_complex*)malloc(experiment->ncs_doppler_cpi*experiment->doppler_padding_factor*sizeof(fftw_complex));

	dopplerData     = (fftw_complex*)malloc(experiment->ncs_padded*experiment->ncs_doppler_cpi*sizeof(fftw_complex));
	
	matchedImageBuffer  = (double*)malloc(experiment->n_threads*experiment->ncs_padded*sizeof(double));
	dopplerImageBuffer  = (double*)malloc(experiment->ncs_doppler_cpi*experiment->doppler_padding_factor*sizeof(double));	
	spectroImageBuffer  = (double*)malloc(experiment->ncs_doppler_cpi*experiment->doppler_padding_factor*sizeof(double));
	
	rangePlan = (fftw_plan*)malloc(experiment->n_threads*sizeof(fftw_plan));
	resultPlan = (fftw_plan*)malloc(experiment->n_threads*sizeof(fftw_plan));	
	
	logger.write("Memory Allocated", timer);		
}


int SignalProcessor::getBlankedPeak(int thread_id)
{
	double max = 0;
	double value = 0;
	
	for (int j = 0; j < experiment->ncs_padded; j++)
	{
		value = 20*log10(mag(lineBuffer[j + thread_id*experiment->ncs_padded]));
		
		if (j == 0)
			max = value;
		else if (value > max)
			max = value;
	}
	
	return max;
}


void SignalProcessor::createPlans(int thread_id)
{
	rangePlan[thread_id] = fftw_plan_dft_1d(experiment->ncs_padded, &lineBuffer[thread_id*experiment->ncs_padded], &lineBuffer[thread_id*experiment->ncs_padded], FFTW_FORWARD, FFTW_MEASURE);
	resultPlan[thread_id] = fftw_plan_dft_1d(experiment->ncs_padded, &resultBuffer[thread_id*experiment->ncs_padded], &lineBuffer[thread_id*experiment->ncs_padded], FFTW_BACKWARD, FFTW_MEASURE);
	dopplerPlan = fftw_plan_dft_1d(experiment->ncs_doppler_cpi*experiment->doppler_padding_factor, dopplerBuffer, dopplerBuffer, FFTW_FORWARD , FFTW_MEASURE);
}


void SignalProcessor::destroyPlans(void)
{	
	fftw_destroy_plan(refPlan);			
	fftw_destroy_plan(dopplerPlan);
	
	for (int i = 0; i < experiment->n_threads; i++)
	{
		fftw_destroy_plan(rangePlan[i]);	
		fftw_destroy_plan(resultPlan[i]);
	}
}


void SignalProcessor::fftRefData(void)
{			
	refPlan = fftw_plan_dft_1d(experiment->ncs_padded, refBuffer, refBuffer, FFTW_FORWARD, FFTW_ESTIMATE);	
	fftw_execute(refPlan);	
	//gnu_plot.gnuPlot(refBuffer, "reference waveform frequency domain", FFT_SHIFT);
}


void SignalProcessor::fftRangeData(int thread_id)
{	
	fftw_execute(rangePlan[thread_id]);	
	//gnu_plot.gnuPlot(rangeBuffer, "return waveform frequency domain", FFT_SHIFT);
}


void SignalProcessor::fftDopplerData(void)
{
	fftw_execute(dopplerPlan);
}


void SignalProcessor::ifftMatchedData(int thread_id)
{
	fftw_execute(resultPlan[thread_id]);
	//gnu_plot.gnuPlot(rangeBuffer, "matched result time domain", NORMAL, AMPLITUDE);
}


void SignalProcessor::popDopplerData(int rangeLine)
{
	if (rangeLine%experiment->update_rate == 0)
		dopplerDataStart = rangeLine;

	if ((rangeLine - dopplerDataStart + 1) <= experiment->ncs_doppler_cpi)
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
		// experiment->pulse_blanking - experiment->ncs_reference corresponds to the resultant 
		// zero range after pulse compression
		for (int i = (experiment->pulse_blanking - experiment->ncs_reference); i < experiment->ncs_padded; i++)		
		{
			popDopplerBuffer(i);	
			fftDopplerData();
			addToDopplerPlot(plot);
			
			if (i == experiment->specro_range_bin)
			{
				plot.addSP(dopplerImageBuffer);
				plot.plotSP();
			}			
		}
		plot.plotRD();
	}
}


void SignalProcessor::popDopplerBuffer(int dopplerLine)
{
	for (int j = 0; j < experiment->ncs_doppler_cpi*experiment->doppler_padding_factor; j++)
	{	
		if (j < experiment->ncs_doppler_cpi)
		{
			float windowCoefficient = dopplerWindow.getSample(j);
			dopplerBuffer[j][0] = dopplerData[dopplerLine*experiment->ncs_doppler_cpi + j][0]*windowCoefficient; 
			dopplerBuffer[j][1] = dopplerData[dopplerLine*experiment->ncs_doppler_cpi + j][1]*windowCoefficient;
		}
		else
		{
			dopplerBuffer[j][0] = 0; 
			dopplerBuffer[j][1] = 0;
		}
	}	
}


void SignalProcessor::addToDopplerPlot(OpenCVPlot &plot)
{
	for (int i = 0; i < experiment->ncs_doppler_cpi*experiment->doppler_padding_factor; i++)
	{
		//perform fft shift
		if (i < ((experiment->ncs_doppler_cpi*experiment->doppler_padding_factor)/2 + 1))		
			dopplerImageBuffer[i + ((experiment->ncs_doppler_cpi*experiment->doppler_padding_factor)/2 - 1)] = mag(dopplerBuffer[i]);
		else
			dopplerImageBuffer[i - ((experiment->ncs_doppler_cpi*experiment->doppler_padding_factor)/2 + 1)] = mag(dopplerBuffer[i]);
	}	
	
	plot.addRD(dopplerImageBuffer);
}


void SignalProcessor::complxConjRef(void)
{
	for (int i = 0; i < (experiment->ncs_padded); i++)
		refBuffer[i][1] = -1*refBuffer[i][1];
		
	logger.write("Complex Conjugate Reference", timer);	
}


void SignalProcessor::complxMulti(int thread_id)
{
	for (int j = 0; j < (experiment->ncs_padded); j++)
	{			
		int k = j + thread_id*experiment->ncs_padded;
		resultBuffer[k][0] = (lineBuffer[k][0]*refBuffer[j][0] - lineBuffer[k][1]*refBuffer[j][1]);
		resultBuffer[k][1] = (lineBuffer[k][0]*refBuffer[j][1] + lineBuffer[k][1]*refBuffer[j][0]);
	}
	//plot.gnuPlot(hilbertBuffer, "matched result frequency domain", FFT_SHIFT);
}


//process data extracted from the bin file into the complex lineBuffer line by line.
void SignalProcessor::popRangeBuffer(int rangeLine, int thread_id)
{
	int start = rangeLine*2*experiment->ncs_range_line;
	
	//populate complex range data and remove offset	
	for (int i = 0; i < experiment->ncs_padded; i++)
	{
		if ((i > experiment->pulse_blanking) && (i < experiment->ncs_range_line))
		{	
			float windowCoefficient = rangeWindow.getSample(i);
			lineBuffer[i + thread_id*experiment->ncs_padded][0] = binDataset[i*2 + start    ]*windowCoefficient;     //real component    
			lineBuffer[i + thread_id*experiment->ncs_padded][1] = binDataset[i*2 + start + 1]*windowCoefficient;     //complex component
		}
		else
		{
			lineBuffer[i + thread_id*experiment->ncs_padded][0] = 0;
			lineBuffer[i + thread_id*experiment->ncs_padded][1] = 0; 
		}
	}
	
	if ((experiment->is_debug) && (rangeLine == 1))
	{
		gPlot.plot(lineBuffer, experiment->ncs_padded, "Time Domain Pulse #1", NORMAL, IQ, experiment->save_path);	
	}	
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
void SignalProcessor::addToWaterPlot(int rangeLine, OpenCVPlot &plot, int thread_id)
{
	for (int j = 0; j < experiment->ncs_padded; j++)
	{
		double pixel = 20*log10(mag(lineBuffer[j + thread_id*experiment->ncs_padded]));
		
		if (pixel < experiment->blanking_threshold)
			matchedImageBuffer[j + thread_id*experiment->ncs_padded] = experiment->blanking_threshold;
		else
			matchedImageBuffer[j + thread_id*experiment->ncs_padded] = pixel;		
	}	
	
	if ((experiment->is_debug) && (rangeLine == 1))
	{
		gPlot.plot(matchedImageBuffer, experiment->ncs_padded, "Pulse Compressed Pulse #1", experiment->save_path);	
	}
	
	plot.addRTI(rangeLine, &matchedImageBuffer[thread_id*experiment->ncs_padded]);
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
			float windowCoefficient = refWindow.getSample(i);
			refBuffer[i][0] = refDataset[i*2    ]*windowCoefficient;   //real component    
			refBuffer[i][1] = refDataset[i*2 + 1]*windowCoefficient;  	//complex component
		}
		else
		{
			refBuffer[i][0] = 0;
			refBuffer[i][1] = 0; 
		}
	}	

	logger.write("Reference Data Loaded", timer);	
	
	if (experiment->is_debug) 
	{
		gPlot.plot(refBuffer, experiment->ncs_reference, "Reference Waveform IQ", NORMAL, IQ, experiment->save_path);	
	}
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
	
	if (experiment->dataset_filename == "-1")
		experiment->dataset_filename 	= (char *)ini.GetValue("dataset", "data_filename");
		
	if (experiment->reference_filename == "-1")
		experiment->reference_filename 	= (char *)ini.GetValue("dataset", "ref_filename");
	
	if (experiment->n_range_lines == -1)	
		experiment->n_range_lines 		= atoi(ini.GetValue("dataset", "n_range_lines"));	
		
	if (experiment->ncs_range_line == -1)
		experiment->ncs_range_line 		= atoi(ini.GetValue("dataset", "n_cmplx_samples_range_line"));
	
	if (experiment->ncs_reference == -1)
		experiment->ncs_reference 		= atoi(ini.GetValue("dataset", "n_cmplx_samples_ref"));	
	
	if (experiment->ncs_padded == -1)
		experiment->ncs_padded 			= atoi(ini.GetValue("dataset", "n_cmplx_samples_padded"));	
		
	if (experiment->ncs_doppler_cpi == -1)
		experiment->ncs_doppler_cpi 	= atoi(ini.GetValue("processing", "doppler_cpi"));	
		
	if (experiment->doppler_padding_factor == -1)
		experiment->doppler_padding_factor = atoi(ini.GetValue("processing", "doppler_padding_factor"));
		
	if (experiment->specro_range_bin == -1)
		experiment->specro_range_bin = atoi(ini.GetValue("processing", "spectrogram_range_bin"));
	
	if (experiment->n_threads == -1)		
		experiment->n_threads 			= atoi(ini.GetValue("processing", "n_threads"));	
		
	if (experiment->pulse_blanking == -1)
		experiment->pulse_blanking 	= atoi(ini.GetValue("visualisation", "pulse_blanking"));	
		
	if (experiment->blanking_threshold == -1)
		experiment->blanking_threshold 	= atoi(ini.GetValue("visualisation", "plot_baseline"));
		
	if (experiment->dynamic_range == -1)
		experiment->dynamic_range 	= atoi(ini.GetValue("visualisation", "dynamic_range"));


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
		
	experiment->update_rate 	= atoi(ini.GetValue("visualisation", "update_rate"));
	experiment->cm 				= atoi(ini.GetValue("visualisation", "colour_map"));
	experiment->hist_equal 		= atoi(ini.GetValue("visualisation", "histogram_equalization"));
	experiment->slow 			= atoi(ini.GetValue("visualisation", "slow"));
	experiment->threshold 		= atoi(ini.GetValue("visualisation", "threshold"));
	experiment->n_plot_average 	= atoi(ini.GetValue("visualisation", "doppler_averaging"));	
	
	//extract and set windowing functions
	refWindow.init((WindowFunction)atoi(ini.GetValue("processing", "ref_window")), experiment->ncs_reference);	
	rangeWindow.init((WindowFunction)atoi(ini.GetValue("processing", "range_window")), (experiment->ncs_range_line));	
	dopplerWindow.init((WindowFunction)atoi(ini.GetValue("processing", "doppler_window")), experiment->ncs_doppler_cpi);	
	
	//calculate the number of range lines each thread is responsible for.
	experiment->n_range_lines_per_thread = experiment->n_range_lines/experiment->n_threads;
	
	//extract path from dataset location
	boost::filesystem::path experiment_filepath(experiment->dataset_filename);
	boost::filesystem::path filename(experiment_filepath.filename());	
	boost::filesystem::path path("../results");
	
	//make the results folder
	if(!boost::filesystem::exists(path))
	{
		std::string command = "mkdir " + path.string();
		system(command.c_str());
	}
	
	//make the dataset specific folder
	filename = path / filename.stem();
	
	if(!boost::filesystem::exists(filename.string()))
	{
		std::string command = "mkdir " + filename.string();
		system(command.c_str());
	}
	
	experiment->save_path = filename.string();	
	
	if (experiment->n_threads != 1)
	{
		experiment->is_doppler = false;
		printf("---> Doppler processing not available with multiple threads.\n");
	}
}


WindowFunction SignalProcessor::parseWindowOption(char* option)
{
	if (option == "HANNING")
		return HANNING;
	else if (option == "HAMMING")
		return HAMMING;
	else if (option == "UNIFORM")
		return UNIFORM;
	else if (option == "BLACKMAN")
		return BLACKMAN;
	else
	{
		printf("Unknown windowing function: '%s' requested, defaulting to UNIFORM.\n", option);
		return UNIFORM;
	}
}







