cmake_minimum_required (VERSION 3.9)
project (mp2)
find_package( BLAS REQUIRED )
find_package( LAPACK REQUIRED)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executableTR(mp2 main.cc
      		    MOTRANS/ops_cis.cc
		          IO/ops_io.cc
		          MATH/ops_mat.cc
		          SCF/ops_rhf.cc
		          SCF/scf.cc )
target_link_libraries(mp2 blas lapack cblas)
