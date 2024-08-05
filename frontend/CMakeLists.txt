
#project(fronend)


 

#############################################
#               Build fronend.exe           #
#############################################

file(GLOB_RECURSE SRC_FRONTEND ${CMAKE_SOURCE_DIR}/frontend/*.cpp)
include_directories(${CMAKE_SOURCE_DIR}/frontend/)

add_executable(frontend ${SRC_FRONTEND})

target_link_libraries(frontend com-psi)
target_link_libraries(frontend aby3-ML)
target_link_libraries(frontend com-psi_Tests)
target_link_libraries(frontend aby3_Tests)
target_link_libraries(frontend oc::tests_cryptoTools)
target_link_libraries(frontend aby3-RTR)
target_link_libraries(frontend aby3-Basic)
target_link_libraries(frontend aby3-GraphQuery)
target_link_libraries(frontend aby3-Benchmark)

find_package(MPI REQUIRED)  
target_link_libraries(frontend MPI::MPI_CXX)

