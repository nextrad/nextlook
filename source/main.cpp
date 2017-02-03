#include "includes.hpp"
#include "parameters.hpp"
#include "timer.hpp"
#include "window.hpp"
#include "logger.hpp"
#include "plotting.hpp"
#include "signal_processor.hpp"

void loopThroughDataset(void);
void initTerminal(void); 

SignalProcessor signalProcessor;

int main(int argc, char *argv[])
{
	initTerminal();	
	initOpenCV();
	
	signalProcessor.allocateMemory();
	signalProcessor.createPlans();
	signalProcessor.loadBinaryDataset();	
	
	signalProcessor.loadReferenceWaveform();	
	
	signalProcessor.fftRefData();		
	signalProcessor.complxConjRef();
	
	loopThroughDataset();
	
	signalProcessor.freeMemory();
	cv::waitKey(0);

	return 0;
}


void loopThroughDataset(void)
{
	for (int i = 0; i < NUMBER_OF_RANGE_LINES; i++)
	{
		signalProcessor.popRangeBuffer(i);
		signalProcessor.fftRangeData();		
		signalProcessor.complxMulti();			
		signalProcessor.ifftMatchedData();								
		signalProcessor.postProcessMatched(i);

		signalProcessor.popDopplerData(i); 		
		signalProcessor.processDoppler(i); 
	}	
}


void initTerminal(void)
{
	system("clear\n");
	printf("NeXtRAD Quicklook Processor\n");
	printf("---------------------------\n");
}


