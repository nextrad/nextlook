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
		for (int i = 0; i < experiment->ncs_padded; i++)		
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
		if (i < experiment->ncs_range_line)
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
void SignalProcessor::addToWaterPlot(int rangeLine, OpenCVPlot &plot, int thread_id)
{
	for (int j = 0; j < experiment->ncs_padded; j++)
	{
		matchedImageBuffer[j + thread_id*experiment->ncs_padded] = mag(lineBuffer[j + thread_id*experiment->ncs_padded]);
	}	
	
	plot.addRTI(rangeLine, &matchedImageBuffer[thread_id*experiment->ncs_padded]);
	
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
	
	experiment->dataset_filename 	= (char *)ini.GetValue("dataset", "data_filename");
	experiment->ncs_range_line 		= atoi(ini.GetValue("dataset", "n_cmplx_samples_range_line"));
	experiment->n_range_lines 		= atoi(ini.GetValue("dataset", "n_range_lines"));	
	experiment->reference_filename 	= (char *)ini.GetValue("dataset", "ref_filename");
	experiment->ncs_reference 		= atoi(ini.GetValue("dataset", "n_cmplx_samples_ref"));	
	experiment->ncs_padded 			= atoi(ini.GetValue("dataset", "n_cmplx_samples_padded"));	
	experiment->n_threads 			= atoi(ini.GetValue("processing", "n_threads"));	
	experiment->ncs_doppler_cpi 	= atoi(ini.GetValue("processing", "doppler_cpi"));	
	experiment->doppler_padding_factor = atoi(ini.GetValue("processing", "doppler_padding_factor"));
	experiment->specro_range_bin = atoi(ini.GetValue("processing", "spectrogram_range_bin"));
	
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
	experiment->cm_doppler 		= atoi(ini.GetValue("visualisation", "doppler_colour_map"));
	experiment->cm_rti 			= atoi(ini.GetValue("visualisation", "rti_colour_map"));
	experiment->hist_equal 		= atoi(ini.GetValue("visualisation", "histogram_equalization"));
	experiment->slow 			= atoi(ini.GetValue("visualisation", "slow"));
	experiment->threshold 		= atoi(ini.GetValue("visualisation", "threshold"));
	experiment->n_plot_average 	= atoi(ini.GetValue("visualisation", "doppler_averaging"));	
	
	//extract and set windowing functions
	refWindow.init((WindowFunction)atoi(ini.GetValue("processing", "ref_window")), experiment->ncs_reference);	
	rangeWindow.init((WindowFunction)atoi(ini.GetValue("processing", "range_window")), experiment->ncs_range_line);	
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







