#TARGET_LINK_LIBRARIES (tsqr )
#LINK_DIRECTORIES ()

# set source list

SET (CMAKE_LIBRARY_PATH $CMAKE_LIBRARY_PATH:"/usr/apps.sp3/isv/intel/xe2013.1.046/composer_xe_2013_sp1.2.144/mkl/lib/intel64")

#INCLUDE_DIRECTORIES(/opt/local/include/openmpi-devel-gcc45)
#INCLUDE_DIRECTORIES(/usr/include/mpi)
file(GLOB_RECURSE SOURCES "*.c")
#file(GLOB_RECURSE HEADERS "src/*.h")
ADD_EXECUTABLE(tsqr ${SOURCES})

# Require MPI for this project:
find_package(MPI REQUIRED)

# For openMPI+Intel+default MPI on tsubame2.5 (default MPI)
include_directories(${MPI_INCLUDE_PATH} "/usr/apps.sp3/free/papi/4.2.1/include")

SET(ICC_COVERAGE_COMPILE_FLAGS "-mkl=sequential")
SET(ICC_COVERAGE_LINK_FLAGS    "-lmkl_lapack95_lp64")

#SET(MPI_INCLUDE_PATH "/usr/apps.sp3/mpi/openmpi/1.6.5/p14.6/include")

SET( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${ICC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${ICC_COVERAGE_LINK_FLAGS} -lpapi -L/usr/apps.sp3/free/papi/4.2.1/lib" )

#set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS})
target_link_libraries(tsqr ${MPI_LIBRARIES})
