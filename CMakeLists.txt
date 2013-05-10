cmake_minimum_required(VERSION 2.8)
project(hammad_mpm)
include(SBELUtils.cmake)

# Additional packages
#enable_cuda_support()
enable_openmp_support()
find_package("GLUT")
	find_package("OpenGL")
set (LINK_LIBS
		${FREEGLUT_LIB}
		${GLUT_LIBRARY}
		${OPENGL_LIBRARIES})
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")

# Add executables
ADD_EXECUTABLE(mpm main.cpp model.cpp)		
TARGET_LINK_LIBRARIES(mpm	${LINK_LIBS})
