#ifndef SIGNAL_PROCESSOR_HPP
#define SIGNAL_PROCESSOR_HPP

#include <fftw3.h>

#include "parameters.hpp"
#include "plotting.hpp"
#include "timer.hpp"
#include "window.hpp"
#include "logger.hpp" 

typedef struct
{
	bool is_doppler;	
	bool is_debug;		
	int n_threads;
	int ncs_reference;		
	int ncs_range_line;		
	int ncs_padded; 		
	int n_range_lines; 		
	int ns_doppler_cpi; 	
	int update_line;

	char* dataset_filename; 	
	char* reference_filename; 	

} Experiment;

class SignalProcessor
{
	private:
		//buffer pointers
		int16_t 		*binDataset; 
		int16_t 		*refDataset;
		fftw_complex	*refBuffer;
		fftw_complex 	*rangeBuffer;
		fftw_complex 	*resultBuffer;
		fftw_complex 	*dopplerBuffer;
		fftw_complex 	*dopplerData;
		
		uint8_t  		*matchedImageBuffer;
		uint8_t  		*dopplerImageBuffer;
		
		Timer timer;
		Logger logger;
		Plot plot;
		Experiment experiment;
		
		//fftw plans
		fftw_plan rangePlan;
		fftw_plan refPlan;
		fftw_plan resultPlan;
		fftw_plan dopplerPlan;		
		
	public:
		//constructor
		SignalProcessor(void);
		
		//memory management
		void allocateMemory(void);
		void freeMemory(void);
		
		//data extraction
		void loadBinaryDataset(void);
		void loadReferenceWaveform(void);
		
		//plan management
		void createPlans(void);
		void destroyPlans(void);
		
		//transforms
		void fftDopplerData(void);
		void fftRefData(void);
		void fftRangeData(void);
		void ifftMatchedData(void);

		void popDopplerData(int rangeLine);
		void popDopplerBuffer(int dopplerLine);		
		void processDoppler(int rangeLine);
		void postProcessDoppler(void);
		
		void complxConjRef(void);
		void complxMulti(void);		
		void postProcessMatched(int rangeLine);
		void popRangeBuffer(int rangeLine);	
		
		void getExperimentParameters(void);
};


#endif


