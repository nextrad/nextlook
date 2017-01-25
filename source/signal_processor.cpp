#include "signal_processor.hpp"

int	dopplerDataStart = 0; 	
Window rangeWindow(HANNING, COMPLEX_SAMPLES_PER_RANGE_LINE);			
Window dopplerWindow(HANNING, RANGE_LINES_PER_DOPPLER_CPI);			

SignalProcessor::SignalProcessor(void)
{	
	timer.start();
}

void SignalProcessor::allocateMemory(void)
{
	binDataset 		= (int16_t*)malloc(NUMBER_OF_RANGE_LINES*VALUES_PER_RANGE_LINE*sizeof(int16_t));
	refDataset      = (int16_t*)malloc(COMPLEX_SAMPLES_IN_REFERENCE_WAVEFORM*2*sizeof(int16_t));
	
	rangeBuffer 	= (fftw_complex*)malloc(COMPLEX_SAMPLES_AFTER_ZERO_PADDING*sizeof(fftw_complex));
	refBuffer 		= (fftw_complex*)malloc(COMPLEX_SAMPLES_AFTER_ZERO_PADDING*sizeof(fftw_complex));
	
	resultBuffer 	= (fftw_complex*)malloc(COMPLEX_SAMPLES_AFTER_ZERO_PADDING*sizeof(fftw_complex));
	dopplerBuffer   = (fftw_complex*)malloc(RANGE_LINES_PER_DOPPLER_CPI*sizeof(fftw_complex));

	dopplerData     = (fftw_complex*)malloc(COMPLEX_SAMPLES_AFTER_ZERO_PADDING*RANGE_LINES_PER_DOPPLER_CPI*sizeof(fftw_complex));
	
	matchedImageBuffer  = (uint8_t*)malloc(COMPLEX_SAMPLES_AFTER_ZERO_PADDING*sizeof(uint8_t));
	dopplerImageBuffer  = (uint8_t*)malloc(RANGE_LINES_PER_DOPPLER_CPI*sizeof(uint8_t));
	
	logger.write("Memory Allocated", timer);
}


void SignalProcessor::createPlans(void)
{	
	//init fftw for multi-threading
	fftw_init_threads();
	fftw_plan_with_nthreads(THREADS);
	
	//fftw plans
	rangePlan   = fftw_plan_dft_1d(COMPLEX_SAMPLES_AFTER_ZERO_PADDING, rangeBuffer  , rangeBuffer  , FFTW_FORWARD , FFTW_MEASURE );
	refPlan     = fftw_plan_dft_1d(COMPLEX_SAMPLES_AFTER_ZERO_PADDING, refBuffer    , refBuffer    , FFTW_FORWARD , FFTW_ESTIMATE);
	resultPlan  = fftw_plan_dft_1d(COMPLEX_SAMPLES_AFTER_ZERO_PADDING, resultBuffer , rangeBuffer  , FFTW_BACKWARD, FFTW_MEASURE );
	dopplerPlan = fftw_plan_dft_1d(RANGE_LINES_PER_DOPPLER_CPI   	 , dopplerBuffer, dopplerBuffer, FFTW_FORWARD , FFTW_MEASURE );	
	
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
	plot.gnuPlot(refBuffer, "reference waveform frequency domain", FFT_SHIFT);
}


void SignalProcessor::fftRangeData(void)
{	
	fftw_execute(rangePlan);	
	plot.gnuPlot(rangeBuffer, "return waveform frequency domain", FFT_SHIFT);
}


void SignalProcessor::fftDopplerData(void)
{
	fftw_execute(dopplerPlan);
}


void SignalProcessor::ifftMatchedData(void)
{
	fftw_execute(resultPlan);
	plot.gnuPlot(rangeBuffer, "matched result time domain", NORMAL, AMPLITUDE);
}


void SignalProcessor::popDopplerData(int rangeLine)
{
	if (rangeLine%UPDATE_LINE == 0)
		dopplerDataStart = rangeLine;

	if ((rangeLine + 1 - dopplerDataStart) <= RANGE_LINES_PER_DOPPLER_CPI)
	{
		for (int j = 0; j < COMPLEX_SAMPLES_AFTER_ZERO_PADDING; j++)
		{
			dopplerData[j*RANGE_LINES_PER_DOPPLER_CPI + (rangeLine - dopplerDataStart)][0] = rangeBuffer[j][0];
			dopplerData[j*RANGE_LINES_PER_DOPPLER_CPI + (rangeLine - dopplerDataStart)][1] = rangeBuffer[j][1];
		}
	}
}


void SignalProcessor::processDoppler(int rangeLine)
{
	if ((rangeLine + 1 - dopplerDataStart) == RANGE_LINES_PER_DOPPLER_CPI)  //check that dopplerData is full
	{
		for (int i = 0; i < COMPLEX_SAMPLES_AFTER_ZERO_PADDING/2 + 1; i++)		
		{
			popDopplerBuffer(i);	
			fftDopplerData();
			postProcessDoppler();
			updateDoppler(dopplerImageBuffer);	
		}
		plotDoppler();
	}
}


void SignalProcessor::popDopplerBuffer(int dopplerLine)
{
	for (int j = 0; j < RANGE_LINES_PER_DOPPLER_CPI; j++)
	{	
		dopplerBuffer[j][0] = dopplerData[dopplerLine*RANGE_LINES_PER_DOPPLER_CPI + j][0]*dopplerWindow.getSample(j); 
		dopplerBuffer[j][1] = dopplerData[dopplerLine*RANGE_LINES_PER_DOPPLER_CPI + j][1]*dopplerWindow.getSample(j);
	}	
}


void SignalProcessor::postProcessDoppler(void)
{
	float maxResult = 0.0f;	
	float result = 0.0f;
	int processed = 0;

	//find max result
	for (int i = 0; i < RANGE_LINES_PER_DOPPLER_CPI; i++)
	{
		result = (sqrt(dopplerBuffer[i][0]*dopplerBuffer[i][0] + dopplerBuffer[i][1]*dopplerBuffer[i][1]));

		if (result > maxResult)
			maxResult = result;
	}

	for (int i = 0; i < RANGE_LINES_PER_DOPPLER_CPI; i++)
	{
		processed = (uint8_t)(((sqrt(dopplerBuffer[i][0]*dopplerBuffer[i][0] + dopplerBuffer[i][1]*dopplerBuffer[i][1]))/maxResult)*255);

		//perform fft shift
		if (i < (RANGE_LINES_PER_DOPPLER_CPI/2 + 1))		
			dopplerImageBuffer[i + (RANGE_LINES_PER_DOPPLER_CPI/2 - 1)] = processed;
		else
			dopplerImageBuffer[i - (RANGE_LINES_PER_DOPPLER_CPI/2 + 1)] = processed;
	}
}


void SignalProcessor::complxConjRef(void)
{
	for (int i = 0; i < (COMPLEX_SAMPLES_AFTER_ZERO_PADDING); i++)
		refBuffer[i][1] = -1*refBuffer[i][1];
		
	logger.write("Complex Conjugate Reference", timer);	
}


void SignalProcessor::complxMulti(void)
{
	for (int j = 0; j < (COMPLEX_SAMPLES_AFTER_ZERO_PADDING); j++)
	{			
		resultBuffer[j][0] = (rangeBuffer[j][0]*refBuffer[j][0] - rangeBuffer[j][1]*refBuffer[j][1]);
		resultBuffer[j][1] = (rangeBuffer[j][0]*refBuffer[j][1] + rangeBuffer[j][1]*refBuffer[j][0]);
	}
	//plot.gnuPlot(hilbertBuffer, "matched result frequency domain", FFT_SHIFT);
}


//process data extracted from the bin file into the complex rangeBuffer line by line.
void SignalProcessor::popRangeBuffer(int rangeLine)
{
	int start = rangeLine*COMPLEX_SAMPLES_PER_RANGE_LINE;
		
	//populate complex range data and remove offset	
	for (int i = 0; i < COMPLEX_SAMPLES_AFTER_ZERO_PADDING; i++)
	{
		if (i < COMPLEX_SAMPLES_PER_RANGE_LINE)
		{	
			rangeBuffer[i][0] = binDataset[i*2 + start    ]*rangeWindow.getSample(i);     //real component    
			rangeBuffer[i][1] = binDataset[i*2 + start + 1]*rangeWindow.getSample(i);     //complex component
		}
		else
		{
			rangeBuffer[i][0] = 0;
			rangeBuffer[i][1] = 0; 
		}
	}	
	plot.gnuPlot(rangeBuffer, "return waveform time domain", NORMAL, IQ);
}


void SignalProcessor::freeMemory(void)
{
	//fftw_free(fftRangeBuffer);
	fftw_free(refBuffer);
	fftw_free(rangeBuffer);
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
void SignalProcessor::postProcessMatched(int rangeLine)
{
	float maxResult = 0.0f;
	float magnitude = 0.0f;

	for (int j = 0; j < COMPLEX_SAMPLES_AFTER_ZERO_PADDING; j++)
	{
		magnitude = 10*log(sqrt(pow(rangeBuffer[j][0], 2) + pow(rangeBuffer[j][1], 2)));

		if (magnitude > maxResult)
		{
			maxResult = magnitude;
		}
	}

	for (int j = 0; j < COMPLEX_SAMPLES_AFTER_ZERO_PADDING; j++)
	{
		magnitude = 10*log(sqrt(pow(rangeBuffer[j][0], 2) + pow(rangeBuffer[j][1], 2)));
		matchedImageBuffer[j] = (uint8_t)((magnitude/maxResult)*255);
	}
	
	updateWaterfall(rangeLine, matchedImageBuffer);
}


void SignalProcessor::loadBinaryDataset(void)
{
	//declare a file pointer
	FILE *binFile;	
	
	//assign pointer to file location
	binFile = fopen(DATASET, "rb");
	
	//check that file exists in the location specified
	if (binFile != NULL)
	{
		//read from binary file into buffer
		fread(binDataset, sizeof(int16_t), NUMBER_OF_RANGE_LINES*VALUES_PER_RANGE_LINE, binFile);

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
	refFile = fopen(REFERENCE, "rb");
	
	//check that file exists in the location specified
	if (refFile != NULL)
	{
		fseek(refFile, 2*COMPLEX_START_SAMPLE_IN_REFERENCE_WAVEFORM*sizeof(int16_t), SEEK_SET);
		//read from binary file into buffer
		fread(refDataset, sizeof(int16_t), COMPLEX_SAMPLES_IN_REFERENCE_WAVEFORM*2, refFile);

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
	for (int i = 0; i < COMPLEX_SAMPLES_AFTER_ZERO_PADDING; i++)
	{
		if (i < COMPLEX_SAMPLES_IN_REFERENCE_WAVEFORM)
		{	
			refBuffer[i][0] = refDataset[i*2];     //real component    
			refBuffer[i][1] = refDataset[i*2 + 1]; //complex component
		}
		else
		{
			refBuffer[i][0] = 0;
			refBuffer[i][1] = 0; 
		}
	}	

	logger.write("Reference Data Loaded", timer);	
	plot.gnuPlot(refBuffer, "reference waveform time domain", NORMAL, IQ);	
}







