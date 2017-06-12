#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

#include <boost/thread.hpp>
#include <string.h>

typedef struct
{
	bool is_doppler;	
	bool is_debug;		
	
	int n_threads;
	int ncs_reference;		
	int ncs_range_line;		
	int ncs_padded; 		
	int n_range_lines; 		
	int n_range_lines_per_thread;
	int ncs_doppler_cpi; 	
	int n_plot_average;
	
	int update_rate;	
	int cm_doppler;
	int cm_rti;
	int hist_equal;
	int slow;
	int threshold;

	char* dataset_filename; 	
	char* reference_filename; 	
	std::string save_path;
	
	boost::mutex mutex;
	
} Experiment;

#endif
