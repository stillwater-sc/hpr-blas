########################################################################
# Config.cmake
#
# Authors: Matthias Moller, Nauman Ahmed, Theodore Omtzigt
########################################################################


########################################################################
# Configure the index data type
########################################################################
if(NOT HPR_BLAS_INDEX_TYPE)
  set(HPR_BLAS_INDEX_TYPE "unsigned long int" CACHE STRING
   "Index data type" FORCE)
endif()
set_property(CACHE HPR_BLAS_INDEX_TYPE PROPERTY STRINGS
  "unsigned short"
  "unsigned int"
  "unsigned long"
  "unsigned long long"
  )

########################################################################
# Configure the backend type
########################################################################
if(NOT HPR_BLAS_BACKEND_TYPE)
  set(HPR_BLAS_BACKEND_TYPE "cpu" CACHE STRING
   "Backend type" FORCE)
endif()
set_property(CACHE HPR_BLAS_BACKEND_TYPE PROPERTY STRINGS
  "cpu"
  "cuda"
  "xsmm"
  )


# Include half precision library if needed
if (HAVE_HALF)
  include_directories("${PROJECT_SOURCE_DIR}/external/half")
endif()

# Include universal library if needed
if (HAVE_POSIT OR HAVE_FIXPNT)
  include_directories("${PROJECT_SOURCE_DIR}/external/universal/include")
  set(CMAKE_CXX_STANDARD 17)
else()
  set(CMAKE_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

########################################################################
# CUDA
########################################################################
option(HPR_BLAS_WITH_CUDA "Enable CUDA support" OFF)

if(HPR_BLAS_WITH_CUDA)
  find_package(CUDA QUIET REQUIRED)
  if(CUDA_FOUND)
    enable_language("CUDA")

    set(HPR_BLAS_NVCC_ARCHS_SUPPORTED "")
    if (NOT CUDA_VERSION VERSION_LESS 7.5)
      list(APPEND HPR_BLAS_NVCC_ARCHS_SUPPORTED 53)
    endif()
    if (NOT CUDA_VERSION VERSION_LESS 8.0)
      list(APPEND HPR_BLAS_NVCC_ARCHS_SUPPORTED 60 61)
    endif()
    if (NOT CUDA_VERSION VERSION_LESS 9.0)
      list(APPEND HPR_BLAS_NVCC_ARCHS_SUPPORTED 70)
    endif()
    if (NOT CUDA_VERSION VERSION_LESS 9.2)
      list(APPEND HPR_BLAS_NVCC_ARCHS_SUPPORTED 72)
    endif()
    if (NOT CUDA_VERSION VERSION_LESS 10.0)
      list(APPEND HPR_BLAS_NVCC_ARCHS_SUPPORTED 75)
    endif()
    if (NOT CUDA_VERSION VERSION_LESS 11.0)
      list(APPEND HPR_BLAS_NVCC_ARCHS_SUPPORTED 80)
    endif()

    set(QNS_NVCC_ARCHS ${HPR_BLAS_NVCC_ARCHS_SUPPORTED} CACHE STRING "The SM architectures requested.")
    set(QNS_NVCC_ARCHS_ENABLED ${HPR_BLAS_NVCC_ARCHS} CACHE STRING "The SM architectures to build code for.")

    set(NVCC_FLAGS)
    foreach(ARCH ${HPR_BLAS_NVCC_ARCHS_ENABLED})
      set(HAVE_CUDA_SM_${ARCH} ON)
      list(APPEND NVCC_FLAGS -gencode=arch=compute_${ARCH},code=sm_${ARCH})
    endforeach()

    set(HAVE_CUDA ON)
    set(CUDA_CXX_STANDARD          ${CMAKE_CXX_STANDARD})
    set(CUDA_CXX_STANDARD_REQUIRED ${CMAKE_CXX_STANDARD_REQUIRED})
    set(CUDA_CXX_EXTENSIONS        ${CMAKE_CXX_EXTENSIONS})
    set(CUDA_PROPAGATE_HOST_FLAGS  ON)
    include_directories("${PROJECT_SOURCE_DIR}/external/cutlass/include")
    include_directories("${PROJECT_SOURCE_DIR}/external/cutlass/tools/util/include")
  endif()
endif()

########################################################################
# MPI
########################################################################
option(HPR_BLAS_WITH_MPI "Enable MPI support" OFF)

if(HPR_BLAS_WITH_MPI)
  find_package(MPI QUIET REQUIRED)
  if(MPI_CXX_FOUND)
    set(HAVE_MPI ON)
  endif()
endif()

########################################################################
# OpenMP
########################################################################
option(HPR_BLAS_WITH_OPENMP "Enable OpenMP support" ON)

if(HPR_BLAS_WITH_OPENMP)
  # Apple explicitly disabled OpenMP support in their compilers that
  # are shipped with XCode but there is an easy workaround as
  # described at https://mac.r-project.org/openmp/
  if (CMAKE_C_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_C_COMPILER_ID STREQUAL "Clang" AND ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" OR
      CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    find_path(OpenMP_C_INCLUDE_DIR
      NAMES "omp.h" PATHS /usr/local /opt /opt/local /opt/homebrew PATH_SUFFICES include)
    find_path(OpenMP_CXX_INCLUDE_DIR
      NAMES "omp.h" PATHS /usr/local /opt /opt/local /opt/homebrew PATH_SUFFICES include)
    find_library(OpenMP_libomp_LIBRARY
      NAMES "omp" PATHS /usr/local /opt /opt/local /opt/homebrew PATH_SUFFICES lib)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xclang -fopenmp -I${OpenMP_C_INCLUDE_DIR}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp -I${OpenMP_CXX_INCLUDE_DIR}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_libomp_LIBRARY}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_libomp_LIBRARY}")
  else() 
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_libomp_LIBRARY}")
  endif()
  set(HAVE_OPENMP ON)
endif()

########################################################################
# LibXSMM
########################################################################
option(HPR_BLAS_WITH_XSMM "Enable LibXSMM support" OFF)

if(HPR_BLAS_WITH_XSMM)
  include_directories("${PROJECT_SOURCE_DIR}/external/libxsmm/include")
  set(HAVE_XSMM ON)
endif()

########################################################################
# Summary
########################################################################
message("")
message("---------- Configuration ----------")
message("Build type.........................: ${CMAKE_BUILD_TYPE}")
message("Build shared libraries.............: ${BUILD_SHARED_LIBS}")
message("Build directory....................: ${PROJECT_BINARY_DIR}")
message("Source directory...................: ${PROJECT_SOURCE_DIR}")
message("Install directory..................: ${CMAKE_INSTALL_PREFIX}")

message("")
message("AR command.........................: ${CMAKE_AR}")
message("RANLIB command.....................: ${CMAKE_RANLIB}")

if(CMAKE_C_COMPILER)
  message("")
  message("C compiler.........................: ${CMAKE_C_COMPILER}")
  message("C compiler flags ..................: ${CMAKE_C_FLAGS}")
  message("C compiler flags (debug)...........: ${CMAKE_C_FLAGS_DEBUG}")
  message("C compiler flags (release).........: ${CMAKE_C_FLAGS_RELEASE}")
  message("C compiler flags (release+debug)...: ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
endif()

if(CMAKE_CXX_COMPILER)
  message("")
  message("CXX compiler.......................: ${CMAKE_CXX_COMPILER}")
  message("CXX standard.......................: ${CMAKE_CXX_STANDARD}")
  message("CXX compiler flags ................: ${CMAKE_CXX_FLAGS}")
  message("CXX compiler flags (debug).........: ${CMAKE_CXX_FLAGS_DEBUG}")
  message("CXX compiler flags (release).......: ${CMAKE_CXX_FLAGS_RELEASE}")
  message("CXX compiler flags (release+debug).: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif()

message("")
message("EXE linker flags...................: ${CMAKE_EXE_LINKER_FLAGS}")
message("EXE linker flags (debug)...........: ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
message("EXE linker flags (release).........: ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
message("EXE linker flags (release+debug)...: ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}")

message("")
message("------------- Options -------------")
message("HPR_BLAS_BACKEND_TYPE..............: ${HPR_BLAS_BACKEND_TYPE}")

message("")
message("HPR_BLAS_NVCC_ARCHS_SUPPORTED......: ${HPR_BLAS_NVCC_ARCHS_SUPPORTED}")
message("HPR_BLAS_WITH_CUDA.................: ${HPR_BLAS_WITH_CUDA}")
message("HPR_BLAS_WITH_MPI..................: ${HPR_BLAS_WITH_MPI}")
message("HPR_BLAS_WITH_OPENMP...............: ${HPR_BLAS_WITH_OPENMP}")
message("HPR_BLAS_WITH_XSMM.................: ${HPR_BLAS_WITH_XSMM}")

message("")
message("HAVE_CUDA..........................: " ${HAVE_CUDA})
message("HAVE_FIXPNT........................: " ${HAVE_FIXPNT})
message("HAVE_HALF..........................: " ${HAVE_HALF})
message("HAVE_MPI...........................: " ${HAVE_MPI})
message("HAVE_OPENMP........................: " ${HAVE_OPENMP})
message("HAVE_POSIT.........................: " ${HAVE_POSIT})
message("HAVE_XSMM..........................: " ${HAVE_XSMM})
