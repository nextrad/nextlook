#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

typedef struct
{
	bool is_doppler;	
	bool is_debug;		
	
	int n_threads;
	int ncs_reference;		
	int ncs_range_line;		
	int ncs_padded; 		
	int n_range_lines; 		
	int ncs_doppler_cpi; 	
	
	int update_rate;	
	int cm_doppler;
	int cm_rti;
	int hist_equal;
	int slow;
	int threshold;

	char* dataset_filename; 	
	char* reference_filename; 	
	
} Experiment;

#endif
