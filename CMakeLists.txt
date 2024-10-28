# Minimum required CMake version and project name
cmake_minimum_required(VERSION 3.10)
project(coor2graph LANGUAGES Fortran)

# Configure basic compilation options and default configuration
set(CMAKE_Fortran_FLAGS "-Wall")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")

# Choose between static or dynamic linkig
option(BUILD_STATIC "Enable static linking" OFF)
if(BUILD_STATIC)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -static")
endif()

# Define output directories for binaries and intermediate files
set(BIN_DIR "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Add source files (modules and executable)
file(GLOB SRC_FILES "${CMAKE_SOURCE_DIR}/src/*.f90")

# Create the executable from the source files
add_executable(coor2graph ${SRC_FILES})

# Display configuration messages
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
if(BUILD_STATIC)
  message(STATUS "Static linking enabled")
else()
	message(STATUS "Dynamic linking enabled")
endif()

# Ensure the output directory for the executable is created
install(TARGETS coor2graph DESTINATION ${BIN_DIR})

