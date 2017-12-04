#include "timer.hpp"
#include "window.hpp"
#include "logger.hpp"
#include "plotting.hpp"
#include "signal_processor.hpp"

void help(void);
void perThread(int id);
void initTerminal(void); 
void parse_options(int argc, char *argv[]);

Experiment experiment;
SignalProcessor signalProcessor(&experiment);
OpenCVPlot opencvPlot(&experiment);
GNUPlot gnuPlot(&experiment);

int main(int argc, char *argv[])
{
	boost::thread_group threadGroup;	
	
	initTerminal();	
	
	parse_options(argc, argv);
	
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
	
	opencvPlot.plotRTI();
	
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
	printf("NeXtLook\n");
	printf("--------\n");
}


void parse_options(int argc, char *argv[])
{
	int opt;

    while ((opt = getopt(argc, argv, "hd:r:p:s:k:z:")) != -1 )
    {
        switch (opt)
        {
            case 'h':
				help();
				break;
            case 'd':
                experiment.dataset_filename = optarg;
                break;
			case 'r':
                experiment.reference_filename = optarg;
                break;
            case 'p':
                experiment.n_range_lines = atoi(optarg);
                break; 
            case 's':
                experiment.ncs_range_line = atoi(optarg);
                break; 
            case 'z':
                experiment.ncs_padded = atoi(optarg);
                break;
            case 'k':
                experiment.ncs_reference = atoi(optarg);
                break;                             
            case '?':
				printf("Unknown command line option.\n");
				exit(EXIT_FAILURE);
        }        
    }
}


void help(void)
{
	printf(" -h: display this help screen\n");
	printf(" -d: dataset_filename\n");
	printf(" -r: reference_filename\n");
	printf(" -p: n_range_lines (pulses)\n");
	printf(" -s: ncs_range_line (complex range bins)\n");
	printf(" -z: ncs_padded (complex range bins padded)\n");
	printf(" -k: ncs_reference (complex samples in reference)\n");
	exit(EXIT_SUCCESS);	
}

