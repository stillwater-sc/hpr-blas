############################
#This configuration file defines some cmake variables:
#HPR_BLAS_INCLUDE_DIRS: list of include directories for the universal library
#HPR_BLAS_LIBRARIES: libraries needed for interfaces like umfpack and arprec, see below
#HPR_BLAS_CXX_DEFINITIONS: definitions to enable the requested interfaces
#HPR_BLAS_VERSION: version (current: 1)
#HPR_BLAS_MINOR_VERSION: minor version 
#
#supported components:
#

unset(HPR_BLAS_LIBRARIES )
unset(HPR_BLAS_CXX_DEFINITIONS )
unset(HPR_BLAS_INCLUDE_DIRS )

if (MSVC)
    add_definitions(/wd4522) # multiple assignment ops for single type, to be investigated further if avoidable
endif()

if (USE_ASSERTS)
  list(APPEND HPR_BLAS_CXX_DEFINITIONS "-DHPR_BLAS_ASSERT_FOR_THROW")
endif()

if(EXISTS ${HPR_BLAS_DIR}/blas/blas.hpp)
	list(APPEND HPR_BLAS_INCLUDE_DIRS "${HPR_BLAS_DIR}/blas")
else()
	list(APPEND HPR_BLAS_INCLUDE_DIRS "${HPR_BLAS_DIR}/../../include")
endif(EXISTS ${HPR_BLAS_DIR}/blas/blas.hpp)

macro(hpr_blas_check_cxx_compiler_flag FLAG RESULT)
  # counts entirely on compiler's return code, maybe better to combine it with check_cxx_compiler_flag
  file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx" "int main() { return 0;}\n")
  try_compile(${RESULT}
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx
    COMPILE_DEFINITIONS ${FLAG})  
endmacro()

