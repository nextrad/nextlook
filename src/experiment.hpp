#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

#include <boost/thread.hpp>
#include <string.h>

typedef struct
{
	bool is_doppler;	
	bool is_debug;		
	bool is_move_file;
	bool is_blanking;
	bool is_visualisation;
	
	int n_threads;
	int ncs_reference;		
	int ncs_range_line;	
	int ncs_range_line_image;	
	int ncs_padded; 		
	int n_range_lines; 		
	int n_range_lines_per_thread;
	int ncs_doppler_cpi; 	
	int n_plot_average;
	int dynamic_range;
	int ncs_blank;
	
	int doppler_padding_factor;
	int specro_range_bin;
	int blanking_threshold;
	int update_rate;	
	int cm;
	int hist_equal;
	int slow;
	int threshold;
	
	int year;
	int month;
	int day;
	int hour;
	int minute;
	int second;
	
	int node_id;
	int adc_channel;

	std::string dataset_filenames[3]; 	
	std::string reference_filename; 	
	std::string output_filenames[3];
	std::string save_path;
	
	boost::mutex mutex;
	
} Experiment;

#endif
