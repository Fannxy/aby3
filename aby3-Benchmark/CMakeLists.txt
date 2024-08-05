file(GLOB_RECURSE SRCS *.cpp)

add_library(aby3-Benchmark ${SRCS})

target_link_libraries(aby3-Benchmark aby3 aby3-Basic aby3-RTR)

option(MPI_APP "Enable MPI support" OFF)

if(MPI_APP)
    target_compile_definitions(aby3-Benchmark PRIVATE "MPI_APP=ON")
endif()