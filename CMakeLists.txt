cmake_minimum_required(VERSION 2.8)
project(quicklook)

find_package(OpenCV REQUIRED)
find_library(FFTW_LIBRARY NAMES fftw3)

file(GLOB CPPFILES "source/*.cpp" "source/*.h")
add_executable(quicklook ${CPPFILES})

target_link_libraries(quicklook ${OpenCV_LIBS} ${FFTW_LIBRARY} -lfftw3_threads)
