#include "signal_processor.hpp"

SignalProcessor::SignalProcessor(Experiment* exp)
{	
	timer.start();
	dopplerDataStart = 0;
	experiment = exp;
	
	experiment->dataset_filenames[0] = "-1";
	experiment->dataset_filenames[1] = "-1";
	experiment->dataset_filenames[2] = "-1";
	
	experiment->reference_filename = "-1";	
	
	experiment->output_filenames[0] = "-1";
	experiment->output_filenames[1] = "-1";
	experiment->output_filenames[2] = "-1";
	
	experiment->n_range_lines = -1;		
	experiment->ncs_range_line = -1;	
	experiment->ncs_reference = -1;	
	experiment->ncs_padded = -1;			
	experiment->ncs_doppler_cpi = -1;		
	experiment->doppler_padding_factor = -1;		
	experiment->specro_range_bin = -1;	
	experiment->n_threads = -1;		
	experiment->blanking_threshold = -1;
	experiment->dynamic_range = -1;
	experiment->node_id = -1;
	experiment->adc_channel = -1;
	experiment->is_move_file = false;	
	experiment->is_blanking = false;	
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
		// doppler plot chopped to have number of samples in pulse less in the front and half
		//of that from the back
		for (int i = experiment->ncs_blank; i < experiment->ncs_padded - experiment->ncs_blank/2; i++)		
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
void SignalProcessor::popRangeBuffer(int rangeLine, int thread_id, bool is_once_off)
{	
	int sample_index = rangeLine*2*experiment->ncs_range_line;	

	//populate complex range data and window
	for (int i = 0; i < experiment->ncs_padded; i++)
	{
		if ((i > experiment->ncs_blank) && (i < experiment->ncs_range_line))
		{	
			float windowCoefficient = rangeWindow.getSample(i);
			lineBuffer[i + thread_id*experiment->ncs_padded][0] = binDataset[i*2 + sample_index    ]*windowCoefficient;     //real component    
			lineBuffer[i + thread_id*experiment->ncs_padded][1] = binDataset[i*2 + sample_index + 1]*windowCoefficient;     //complex component
		}
		else
		{
			lineBuffer[i + thread_id*experiment->ncs_padded][0] = 0;
			lineBuffer[i + thread_id*experiment->ncs_padded][1] = 0; 
		}
	}
	
	if ((experiment->is_debug) && (rangeLine == 0))
	{
		gPlot.plot(lineBuffer, experiment->ncs_padded, "Time Domain Pulse #1", NORMAL, IQ, experiment->save_path);	
	}	
	
	if ((experiment->is_move_file) && (experiment->n_threads == 1))
	{				
		if (!is_once_off)
		{
			//temporary for moving one channel !!!!!!!!!!!!!!!!!!!!!!!!!!
			for (int channel = experiment->adc_channel; channel < experiment->adc_channel + 1; channel++)
			{
				if (rangeLine == 0)
				{
					std::cout << "Opening Output File: " << experiment->output_filenames[channel] << std::endl;
					outFile[channel] = fopen(experiment->output_filenames[channel].c_str(), "wb");
				}
				
				//check that file exists in the location specified
				if (outFile[channel] != NULL)
				{			
					//read from binary file into buffer
					fwrite(&binDataset[sample_index], sizeof(int16_t), experiment->ncs_range_line*2, outFile[channel]);

					if (rangeLine == experiment->n_range_lines - 1)
					{
						logger.write("Finished Copying Dataset.", timer);
						fclose(outFile[channel]);
						
						logger.write("Checking That Copy and Original are Identical.", timer);
						std::stringstream command;
						command << "diff " << experiment->dataset_filenames[channel] << " " << experiment->output_filenames[channel] << "\n";
						
						int ret = system(command.str().c_str());
						if (WEXITSTATUS(ret) == 0)
						{						
							//clear the stringstream
							command.str(std::string());
							
							command << "rm " << experiment->dataset_filenames[channel] << "\n";
							system(command.str().c_str());
							std::cout << "Deleting original file: " << command.str() << std::endl;
						}
						else
						{
							std::cout << "Copied Dataset is NOT Identical to the Original Dataset." << std::endl;
							std::cout << "Original Dataset will NOT be Deleted." << std::endl;
						}
					}
						
				}				
				else //file does not exist in the specified location
				{
					logger.write("Output File Location Cannot Be Written To.");
					exit(EXIT_FAILURE);
				}
			}	
		}
	}
	else
	{
		logger.write("File moving is not supported with multithreading, either disable file moving or specify one thread.");
		exit(EXIT_FAILURE);
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
	binFile = fopen(experiment->dataset_filenames[experiment->adc_channel].c_str(), "rb");			
	
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
		std::cout << "Binary dataset could not be found at: " << experiment->dataset_filenames[experiment->adc_channel] << std::endl;
		exit(EXIT_FAILURE);
	}
}


void SignalProcessor::loadReferenceWaveform(void)
{
	//declare a file pointer
	FILE *refFile;	
	
	//assign pointer to file location
	refFile = fopen(experiment->reference_filename.c_str(), "rb");
	
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
	
	if (experiment->doppler_padding_factor == -1)
		experiment->doppler_padding_factor = atoi(ini.GetValue("processing", "doppler_padding_factor"));
	
	if (experiment->n_threads == -1)		
		experiment->n_threads = atoi(ini.GetValue("processing", "n_threads"));	
		
	std::string move_flag = ini.GetValue("config", "is_move_file");	
		
	if ((move_flag == "1") || (move_flag == "true") || (move_flag == "TRUE") || (move_flag == "True"))
		experiment->is_move_file = true;
	else
		experiment->is_move_file = false;
		
	std::string blinking_flag = ini.GetValue("config", "is_blanking");	
		
	if ((blinking_flag == "1") || (blinking_flag == "true") || (blinking_flag == "TRUE") || (blinking_flag == "True"))
	{
		experiment->is_blanking = true;
		experiment->ncs_blank = experiment->ncs_reference;
	}
	else
	{
		experiment->is_blanking = false;
		experiment->ncs_blank = 0;
	}	

	std::string doppler_flag = ini.GetValue("config", "is_doppler");	
		
	if ((doppler_flag == "1") || (doppler_flag == "true") || (doppler_flag == "TRUE") || (doppler_flag == "True"))
		experiment->is_doppler = true;
	else
		experiment->is_doppler = false;
		
	std::string debug_flag = ini.GetValue("config", "is_debug");	
		
	if ((debug_flag == "1") || (debug_flag == "true") || (debug_flag == "TRUE") || (debug_flag == "True"))
		experiment->is_debug = true;
	else
		experiment->is_debug = false;		
	
	experiment->cm 				= atoi(ini.GetValue("visualisation", "colour_map"));
	experiment->hist_equal 		= atoi(ini.GetValue("visualisation", "histogram_equalization"));
	experiment->slow 			= atoi(ini.GetValue("visualisation", "slow"));
	experiment->threshold 		= atoi(ini.GetValue("visualisation", "threshold"));
	
	//extract and set windowing functions
	refWindow.init((WindowFunction)atoi(ini.GetValue("processing", "ref_window")), experiment->ncs_reference);	
	rangeWindow.init((WindowFunction)atoi(ini.GetValue("processing", "range_window")), (experiment->ncs_range_line));	
	dopplerWindow.init((WindowFunction)atoi(ini.GetValue("processing", "doppler_window")), experiment->ncs_doppler_cpi);	
	
	//calculate the number of range lines each thread is responsible for.
	experiment->n_range_lines_per_thread = experiment->n_range_lines/experiment->n_threads;
	
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


void SignalProcessor::readHeader(void)
{
	//check that filepath is good
	std::ifstream check(HDR_FILE);

	if (!check.good()) 
	{
		printf("Please check location of header file and try again.\n");
		exit(EXIT_FAILURE);
	}
	
	CSimpleIniA ini;
	ini.LoadFile(HDR_FILE);	
	
	experiment->n_range_lines = atoi(ini.GetValue("PulseParameters", "NUM_PRIS"));
	experiment->ncs_range_line = atoi(ini.GetValue("PulseParameters", "SAMPLES_PER_PRI"));
	//no zero pading is currently in place  
	experiment->ncs_padded = experiment->ncs_range_line;
	experiment->dynamic_range = atoi(ini.GetValue("Quicklook", "DYNAMIC_RANGE"));
	
	int waveform_index = atoi(ini.GetValue("PulseParameters", "WAVEFORM_INDEX"));
	int sampling_freq = 90e6;
	
	//get dataset filename
	experiment->adc_channel = atoi(ini.GetValue("Quicklook", "ADC_CHANNEL"));
	
	int amp_index = atoi(ini.GetValue("Quicklook", "AMPLIFIER"));
	
	std::stringstream ss_ref_filename;
	
	ss_ref_filename << REF_DIR;
	
	switch (experiment->adc_channel)
	{
		case 0:
			ss_ref_filename << "L_";
			break;
		case 1:
			ss_ref_filename << "X_";
			break;
		case 2:
			ss_ref_filename << "X_";
			break;			
	}
	
	switch (amp_index)
	{
		case 0:
			ss_ref_filename << "M_";
			break;
		case 1:
			ss_ref_filename << "H_";
			break;		
	}	
	
	ss_ref_filename << waveform_index << ".dat";
	
	experiment->reference_filename = ss_ref_filename.str();
	
	std::cout << experiment->reference_filename << std::endl;
	
	switch (waveform_index)
	{
		case 1:
			experiment->ncs_reference = 0.5e-6*sampling_freq;
			break;
		case 2:
			experiment->ncs_reference = 1.0e-6*sampling_freq;
			break;
		case 3:
			experiment->ncs_reference = 3.0e-6*sampling_freq;
			break;
		case 4:
			experiment->ncs_reference = 5.0e-6*sampling_freq;
			break;
		case 5:
			experiment->ncs_reference = 10.0e-6*sampling_freq;
			break;
		case 6:
			experiment->ncs_reference = 15.0e-6*sampling_freq;
			break;
		case 7:
			experiment->ncs_reference = 20.0e-6*sampling_freq;
			break;
	}
	
	experiment->dataset_filenames[0] = COBALT_ADC_DIR + (std::string)("adc0.dat");
	experiment->dataset_filenames[1] = COBALT_ADC_DIR + (std::string)("adc1.dat");
	experiment->dataset_filenames[2] = COBALT_ADC_DIR + (std::string)("adc2.dat");
	
	//get date and time
	experiment->year 	= atoi(ini.GetValue("Timing", "YEAR"));
	experiment->month 	= atoi(ini.GetValue("Timing", "MONTH"));
	experiment->day 	= atoi(ini.GetValue("Timing", "DAY"));
	experiment->hour 	= atoi(ini.GetValue("Timing", "HOUR"));
	experiment->minute 	= atoi(ini.GetValue("Timing", "MINUTE"));
	experiment->second 	= atoi(ini.GetValue("Timing", "SECOND"));
	
	//generate the output file name
	std::stringstream ss_output_path;
	std::stringstream ss_output_file;
	
	std::string folder_name;

	ss_output_path << experiment->year 	 << "_";
	ss_output_path << experiment->month	 << "_";
	ss_output_path << experiment->day 	 << "_";
	ss_output_path << experiment->hour 	 << "_";
	ss_output_path << experiment->minute << "_";
	ss_output_path << experiment->second << "_";		
	ss_output_path << "n" << experiment->node_id;
	
	//got folder name
	folder_name = ss_output_path.str();
	
	//clear stringstream
	ss_output_path.str(std::string());
	ss_output_path << EXT_STORAGE_DIR << folder_name;
	experiment->save_path = ss_output_path.str();	
	
	for (int i = 0; i < 3; i++)
	{
		ss_output_file.str(std::string());
		ss_output_file << experiment->save_path << "/" << folder_name << "_adc" << i << ".dat";
		experiment->output_filenames[i] = ss_output_file.str();
	}	
	
	std::string command = "mkdir " + experiment->save_path;
	system(command.c_str());	
	
	//temporary move command
	for (int i = 0; i < 3; i++)
	{
		if (i != experiment->adc_channel)
		{
			std::string command = "";
			command = "mv " + experiment->dataset_filenames[i] + " " + experiment->output_filenames[i] + "\n";
			system(command.c_str());
		}
	}

	experiment->specro_range_bin = atoi(ini.GetValue("Quicklook", "SPECTROGRAM_BIN"));
	experiment->ncs_doppler_cpi = atoi(ini.GetValue("Quicklook", "DOPPLER_FFT"));
	
	//hard coded update rate
	experiment->update_rate = experiment->ncs_doppler_cpi;
}







