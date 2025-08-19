file(GLOB_RECURSE SRCS *.cpp)

add_library(aby3-Feddb-Core ${SRCS})
include_directories("/usr/include/jsoncpp")  
include_directories("/usr/include/mpich")

#######################################
# json
if(NOT TARGET jsoncpp)

    find_path(JSON_INCLUDE_DIR 
      NAMES "json/json.h"
      REQUIRED)
    if(NOT EXISTS ${JSON_INCLUDE_DIR})
        message(FATAL_ERROR "failed to find json")
    endif()
    message(STATUS "In aby3-Feddb-Core, JSON_INCLUDE_DIR=${JSON_INCLUDE_DIR}")
    add_library(json INTERFACE IMPORTED)
    target_include_directories(json INTERFACE 
                    $<BUILD_INTERFACE:${JSON_INCLUDE_DIR}>
                    $<INSTALL_INTERFACE:>)
endif()

#######################################
# mpi
if(NOT TARGET mpi)
    find_path(MPI_INCLUDE_DIR 
      NAMES "mpi.h" 
      REQUIRED)
    if(NOT EXISTS ${MPI_INCLUDE_DIR})
        message(FATAL_ERROR "failed to find mpi.h")
    endif()
    message(STATUS "In aby3-Feddb-Core, MPI_INCLUDE_DIR=${MPI_INCLUDE_DIR}\n\n")
    add_library(mpi INTERFACE IMPORTED)
    target_include_directories(mpi INTERFACE 
                    $<BUILD_INTERFACE:${MPI_INCLUDE_DIR}>
                    $<INSTALL_INTERFACE:>)
endif()

target_link_libraries(aby3-Feddb-Core aby3 jsoncpp mpi)


# file(GLOB_RECURSE SRCS *.cpp)
 

# add_library(aby3-Feddb-Core ${SRCS})


# target_link_libraries(aby3-Feddb-Core aby3)
