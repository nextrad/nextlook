#include "timer.hpp"
#include "window.hpp"
#include "logger.hpp"
#include "plotting.hpp"
#include "signal_processor.hpp"

#define CNC_IP 			"192.168.1.100"
#define CNC_SAVE_PATH 	"/home/nextrad/Documents/NeXtLook/"

void help(void);
void perThread(int id);
void calcThreshold(int pulse_number);
void initTerminal(void); 
void parse_options(int argc, char *argv[]);
void pushImages(void);

Experiment experiment;
SignalProcessor signalProcessor(&experiment);
OpenCVPlot opencvPlot(&experiment);

int main(int argc, char *argv[])
{
	boost::thread_group threadGroup;	
	
	initTerminal();	
	
	//get cmd line options
	parse_options(argc, argv);
	
	//get parameters from the header file
	signalProcessor.readHeader();
	
	//get parameters from the experiment file
	//only if they haven't been set using cmd line or the header file
	signalProcessor.getExperimentParameters();			
	
	boost::thread threads[experiment.n_threads];
	fftw_make_planner_thread_safe();
	
	signalProcessor.allocateMemory();		
	
	signalProcessor.loadBinaryDataset();		
	signalProcessor.loadReferenceWaveform();	
	
	signalProcessor.fftRefData();		
	signalProcessor.complxConjRef();
	
	calcThreshold(1);	
	opencvPlot.initOpenCV();
	
	for (int i = 0; i < experiment.n_threads; i++)
	{
		threadGroup.create_thread(boost::bind(&perThread, i));
	}
	
	threadGroup.join_all();	
	
	opencvPlot.plotRTI();
	opencvPlot.plotSP();
	
	opencvPlot.savePlots();
	
	pushImages();
	
	signalProcessor.freeMemory();
	
	if (experiment.is_visualisation)
	{
		std::cout << "Press any key to exit..." << std::endl;	
		cv::waitKey(0);
	}

	return 0;
}


void calcThreshold(int pulse_number)
{	
	//pulse compress the first range line and get the peak value after
	//blanking. The threshold = blanked_peak - dynamic_range
	signalProcessor.createPlans(0);	
	signalProcessor.popRangeBuffer(pulse_number, 0, true);		
	signalProcessor.fftRangeData(0);		
	signalProcessor.complxMulti(0);			
	signalProcessor.ifftMatchedData(0);	
	//signalProcessor.destroyPlans();
	
	experiment.blanking_threshold = signalProcessor.getBlankedPeak(0) - experiment.dynamic_range;				
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

    while ((opt = getopt(argc, argv, "hn:vb")) != -1 )
    {
        switch (opt)
        {
            case 'h':
				help();
				break;
            case 'n':
                experiment.node_id = atoi(optarg);
                if ((experiment.node_id != 0) && (experiment.node_id != 1) && (experiment.node_id != 2))
                {
					std::cout << "Node ID can only be 0, 1 or 2." << std::endl;
					exit(EXIT_FAILURE);
				}
                break;   
            case 'v':
                experiment.is_visualisation = true;
                break; 
            case 'b':
                experiment.is_blanking = true;
                break;                         
            case '?':
				printf("Unknown command line option.\n");
				exit(EXIT_FAILURE);
        }        
    }
    
	if (experiment.node_id == -1)
	{
		std::cout << "Node ID has not set." << std::endl << std::endl;
		help();
		exit(EXIT_FAILURE);
	}
}


void help(void)
{
	printf("Usage: ./nextlook [options]\n");
	printf(" -h: display this help screen\n");
	printf(" -n: specify node ID (0, 1 or 2)\n");
	printf(" -v: display plots during processing (requires -X during ssh)\n");	
	printf(" -b: enable range blanking\n");	
	exit(EXIT_SUCCESS);	
}


void pushImages(void)
{
	std::stringstream ss_command;
	
	ss_command.str(std::string());
	ss_command << "scp " << opencvPlot.get_rt_path() << " nextrad@" << CNC_IP << ":" << CNC_SAVE_PATH << "RTI_N" << experiment.node_id << ".jpg\n";
	std::cout << ss_command.str();
	system(ss_command.str().c_str());
	
	ss_command.str(std::string());
	ss_command << "scp " << opencvPlot.get_rd_path() << " nextrad@" << CNC_IP << ":" << CNC_SAVE_PATH << "RD_N" << experiment.node_id << ".jpg\n";
	std::cout << ss_command.str();
	system(ss_command.str().c_str());
	
	ss_command.str(std::string());
	ss_command << "scp " << opencvPlot.get_sp_path() << " nextrad@" << CNC_IP << ":" << CNC_SAVE_PATH << "SP_N" << experiment.node_id << ".jpg\n";
	std::cout << ss_command.str();
	system(ss_command.str().c_str());
}

