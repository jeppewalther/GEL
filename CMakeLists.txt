# add stuff to use MFC in this executable
ADD_DEFINITIONS(-D_AFXDLL)
SET(CMAKE_MFC_FLAG 6)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(GEL)

FILE(GLOB CGLA ./src/CGLA/*.cpp ./src/CGLA/*.c ./src/CGLA/*.h)
FILE(GLOB GLGraphics ./src/GLGraphics/*.cpp ./src/GLGraphics/*.c ./src/GLGraphics/*.h)
FILE(GLOB GL ./src/GL/*.cpp ./src/GL/*.c ./src/GL/*.h)
FILE(GLOB HMesh ./src/HMesh/*.cpp ./src/HMesh/*.c ./src/HMesh/*.h)
FILE(GLOB Geometry ./src/Geometry/*.cpp ./src/Geometry/*.c ./src/Geometry/*.h)
FILE(GLOB Util ./src/Util/*.cpp ./src/Util/*.c ./src/Util/*.h)
FILE(GLOB LinAlg ./src/LinAlg/*.cpp ./src/LinAlg/*.c ./src/LinAlg/*.h)

source_group(CGLA FILES ${CGLA})
source_group(GLGraphics FILES ${GLGraphics})
source_group(GL FILES ${GL})
source_group(HMesh FILES ${HMesh})
source_group(Geometry FILES ${Geometry})
source_group(Util FILES ${Util})
source_group(LinAlg FILES ${LinAlg})

ADD_LIBRARY(GEL ${CGLA} ${HMesh} ${GL} ${GLGraphics} ${Geometry} ${Util} ${LinAlg}) 

IF(APPLE)
SET_TARGET_PROPERTIES(GEL PROPERTIES XCODE_ATTRIBUTE_PRODUCT_NAME GEL)
FIND_LIBRARY(COCOA_LIBRARY Cocoa)
FIND_LIBRARY(OpenGL_LIBRARY OpenGL)
FIND_LIBRARY(GLUT_LIBRARY GLUT)
FIND_LIBRARY(vecLib_LIBRARY vecLib)
SET(EXTRA_LIBS ${COCOA_LIBRARY} ${OpenGL_LIBRARY} ${GLUT_LIBRARY} ${vecLib_LIBRARY})
TARGET_LINK_LIBRARIES(GEL ${EXTRA_LIBS})
SET_TARGET_PROPERTIES(GEL PROPERTIES LINKER_LANGUAGE CXX)
SET_TARGET_PROPERTIES(GEL PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD c++0x)
SET_TARGET_PROPERTIES(GEL PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY libc++)
ELSE(APPLE)

set(DEPENDENCY_DIRS "lib")
set(INCLUDE_DIRS "include")

# Make sure we have the libraries
find_library(LIB_GLUT glut32 ${DEPENDENCY_DIRS})
find_library(LIB_LAPACK clapack ${DEPENDENCY_DIRS})

# Include directories for dependencies
find_path(INCLUDE_GLUT GL/glut.h ${INCLUDE_DIRS})
include_directories(${INCLUDE_GLUT})
TARGET_LINK_LIBRARIES(GEL ${LIB_LAPACK} ${LIB_GLUT})
ENDIF(APPLE)

# Set flags for CLang/C++11
if (("${CMAKE_GENERATOR}" MATCHES "Unix Makefiles") AND APPLE)
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
endif ()

# Remove all the deprecated stuff from VS2005
IF(WIN32)
  OPTION(DISABLE_VS2005DEPRECATED "Disable the deprecated warnings in VS2005" ON)
  IF(DISABLE_VS2005DEPRECATED)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
  ENDIF(DISABLE_VS2005DEPRECATED)
ENDIF(WIN32)

# Append a d to all debug libraries
SET(CMAKE_DEBUG_POSTFIX "d")

OPTION(BUILD_APPLICATIONS "Build Applications" ON)
OPTION(BUILD_TESTING "Build test applications" ON)

INCLUDE_DIRECTORIES(src)

IF(BUILD_APPLICATIONS)
  SUBDIRS(
    apps/MeshDistance
    apps/MeshEdit
    apps/OBJViewer
    apps/Pathtracer
    apps/Raytracer
    )
ENDIF(BUILD_APPLICATIONS)

IF(BUILD_TESTING)
  SUBDIRS(
    test/CGLA-covariance
    test/CGLA-mat
    test/CGLA-ogl
    test/CGLA-simple
    test/CGLA-vec
    test/Geometry-kdtree
    test/GLGraphics-console
    test/LinAlg-leastsq
    test/LinAlg-eigensolver
    test/Util-resman
    )
ENDIF(BUILD_TESTING)

# Try to locate test data
FIND_PATH(GEL_DATA_DIR head.mtl
  ${GEL_SOURCE_DIR}/data
  )

# Now define it so it can be used as a path in test programs
IF(GEL_DATA_DIR)
  ADD_DEFINITIONS(-DGEL_DATA_DIR=${GEL_DATA_DIR})
ENDIF(GEL_DATA_DIR)

IF(NOT LIBRARY_OUTPUT_PATH)
  SET(LIBRARY_OUTPUT_PATH ${GEL_BINARY_DIR}/lib CACHE INTERNAL "Single output directory for building all libraries.")
ENDIF(NOT LIBRARY_OUTPUT_PATH)
IF(NOT EXECUTABLE_OUTPUT_PATH)
  SET(EXECUTABLE_OUTPUT_PATH ${GEL_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
ENDIF(NOT EXECUTABLE_OUTPUT_PATH)

# Configure the file describing how to use GEL
CONFIGURE_FILE(
  ${GEL_SOURCE_DIR}/GELConfig.cmake.in
  ${GEL_BINARY_DIR}/GELConfig.cmake
  @ONLY IMMEDIATE
  )

# Inspired by the VXL install routines
CONFIGURE_FILE( ${GEL_SOURCE_DIR}/GELUse.cmake
                ${GEL_BINARY_DIR}/GELUSE.cmake COPYONLY )
