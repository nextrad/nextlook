#include "timer.hpp"
#include "window.hpp"
#include "logger.hpp"
#include "plotting.hpp"
#include "signal_processor.hpp"

void loopThroughDataset(void);
void initTerminal(void); 

Experiment experiment;
SignalProcessor signalProcessor(&experiment);
OpenCVPlot opencvPlot(&experiment);
GNUPlot gnuPlot(&experiment);

int main(int argc, char *argv[])
{
	initTerminal();	
	
	signalProcessor.getExperimentParameters();		
	
	signalProcessor.allocateMemory();
	signalProcessor.createPlans();
	signalProcessor.loadBinaryDataset();		
	signalProcessor.loadReferenceWaveform();	
	
	signalProcessor.fftRefData();		
	signalProcessor.complxConjRef();
	
	opencvPlot.initOpenCV();
	loopThroughDataset();
	
	signalProcessor.freeMemory();
	cv::waitKey(0);

	return 0;
}


void loopThroughDataset(void)
{
	for (int i = 0; i < experiment.n_range_lines; i++)
	{
		signalProcessor.popRangeBuffer(i);
		signalProcessor.fftRangeData();		
		signalProcessor.complxMulti();			
		signalProcessor.ifftMatchedData();								
		signalProcessor.addToWaterPlot(i, opencvPlot);
		
		signalProcessor.processDoppler(i, opencvPlot); 
	}	
}


void initTerminal(void)
{
	system("clear\n");
	printf("NeXtRAD Quicklook Processor\n");
	printf("---------------------------\n");
}


