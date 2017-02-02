cmake_minimum_required(VERSION 2.8)
project( quicklook )

find_package( OpenCV REQUIRED )
find_library(FFTW_LIBRARY NAMES fftw3)


add_executable( quicklook 
				main.cpp 
				source/plotting.cpp
				source/plotting.hpp
				source/signal_processor.cpp
				source/signal_processor.hpp
				source/utilities/timer.cpp
				source/utilities/timer.hpp
				source/utilities/window.cpp
				source/utilities/window.hpp
				source/utilities/logger.cpp
				source/utilities/logger.hpp
				source/parameters.hpp
				source/includes.hpp
				)


target_link_libraries( 	quicklook 
						${OpenCV_LIBS} 
						${FFTW_LIBRARY}
						-lfftw3_threads
						)
