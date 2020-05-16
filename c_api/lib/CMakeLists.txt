set(lib_name hprblas_c_api)
set(folder Libraries)

add_library(${lib_name} STATIC hprblas_c_api.cpp)
set_target_properties(${lib_name} PROPERTIES FOLDER ${folder})

if(C_API_LIB_PIC)
  if(CMAKE_COMPILER_IS_GNUCXX OR MINGW OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
endif()

install(TARGETS ${lib_name} DESTINATION lib)
install(FILES hprblas_c_api.h DESTINATION include)

# TODO: for adaptive precision algorithms there are hundreds of valuable posit configurations that the C++ version supports
# we will need a stub generator to capture them all in a productive manner
# OR: we simply require that you use the C++ library as C and Fortran environments would have difficulty creating
# adaptive precision algorithms anyways.
