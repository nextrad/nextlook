#include "source/includes.hpp"
#include "source/parameters.hpp"
#include "source/utilities/timer.hpp"
#include "source/utilities/window.hpp"
#include "source/utilities/logger.hpp"
#include "source/plotting.hpp"
#include "source/signal_processor.hpp"

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


