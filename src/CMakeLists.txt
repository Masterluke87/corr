cmake_minimum_required (VERSION 3.9)
project (mp2)

set (CMAKE_CXX_STANDARD 11)

find_package( BLAS REQUIRED )
find_package( CBLAS REQUIRED )
find_package( LAPACK REQUIRED )
FIND_PACKAGE( OpenMP REQUIRED )
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

message("libint include: ${LIBINT_INCLUDE_DIR}")
if(EXISTS "${LIBINT_INCLUDE_DIR}/libint2.hpp")
  message("LIBINT HEADER FOUND")
else()
  message(SEND_ERROR "Can't find libint2.hpp in ${LIBINT_INCLUDE_DIR}")
endif()

if(EXISTS "${LIBINT_LIBRARY_DIR}/libint2.a")
  message("LIBINT LIBRARY FOUND")
else()
  message(SEND_ERROR "Can't find libint2.a in ${LIBINT_LIBRARY_DIR}")
endif()

include_directories(${LIBINT_INCLUDE_DIR})

add_executable(mp2 main.cc
      		    POSTHF/motrans.cc
              POSTHF/loc_mp2.cc
              POSTHF/can_mp2.cc
		          IO/ops_io.cc
		          MATH/ops_mat.cc
		          SCF/ops_rhf.cc
		          SCF/scf.cc )

target_link_libraries(mp2 ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${CBLAS_LIBRARIES} ${LIBINT_LIBRARY_DIR}/libint2.a)
