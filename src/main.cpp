#include "timer.hpp"
#include "window.hpp"
#include "logger.hpp"
#include "plotting.hpp"
#include "signal_processor.hpp"

void perThread(int id);
void initTerminal(void); 

Experiment experiment;
SignalProcessor signalProcessor(&experiment);
OpenCVPlot opencvPlot(&experiment);
GNUPlot gnuPlot(&experiment);

int main(int argc, char *argv[])
{
	boost::thread_group threadGroup;	
	
	initTerminal();	
	
	signalProcessor.getExperimentParameters();			
	
	boost::thread threads[experiment.n_threads];
	fftw_make_planner_thread_safe();
	
	signalProcessor.allocateMemory();		
	
	signalProcessor.loadBinaryDataset();		
	signalProcessor.loadReferenceWaveform();	
	
	signalProcessor.fftRefData();		
	signalProcessor.complxConjRef();
	
	opencvPlot.initOpenCV();
	
	for (int i = 0; i < experiment.n_threads; i++)
	{
		threadGroup.create_thread(boost::bind(&perThread, i));
	}
	
	threadGroup.join_all();	
	
	opencvPlot.plotWaterfall();
	
	opencvPlot.savePlots();
	
	signalProcessor.freeMemory();
	cv::waitKey(0);

	return 0;
}


void perThread(int id)
{
	int start_index = id*experiment.n_range_lines_per_thread;
	
	signalProcessor.createPlans(id);
	
	for (int i = start_index; i < start_index + experiment.n_range_lines_per_thread; i++)
	{
		signalProcessor.popRangeBuffer(i, id);
		signalProcessor.fftRangeData(id);		
		signalProcessor.complxMulti(id);			
		signalProcessor.ifftMatchedData(id);								
		signalProcessor.addToWaterPlot(i, opencvPlot, id);
		
		if ((experiment.is_doppler) && (experiment.n_threads == 1))
		{
			signalProcessor.processDoppler(i, opencvPlot); 
		}
	}	
}


void initTerminal(void)
{
	system("clear\n");
	printf("NeXtlook\n");
	printf("--------\n");
}


