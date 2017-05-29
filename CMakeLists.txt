cmake_minimum_required(VERSION 2.8)
project(nextlook)

find_package(OpenCV REQUIRED)
find_library(FFTW_LIBRARY NAMES fftw3)
find_package(Boost COMPONENTS thread date_time program_options filesystem system REQUIRED)

file(GLOB CPPFILES "src/*.cpp" "src/*.h")
add_executable(nextlook ${CPPFILES})

target_link_libraries(nextlook ${OpenCV_LIBS} ${FFTW_LIBRARY} ${Boost_LIBRARIES} -lfftw3_threads -pthread)
