#---------------------------------------------------
# OpenMM CustomAnisotropicNonbonded Plugin
# modified for CAN tjanicki 5.12.17
# Please note that all subsequent CMake files are adapted from
# 	OpenMM use
#----------------------------------------------------

#find_package(PythonInterp)
#mark_as_advanced(CLEAR PYTHON_EXECUTABLE)

CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
SET(CMAKE_CXX_STANDARD 11)
#SET(CMAKE_BUILD_TYPE DEBUG)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

# We need to know where OpenMM is installed so we can access the headers and libraries.
SET(OPENMM_DIR "/usr/local/openmm" CACHE PATH "Where OpenMM is installed")
SET(OPENMM_SOURCE "/usr/local/openmm" CACHE PATH "Where OpenMM source files are located")
INCLUDE_DIRECTORIES("${OPENMM_DIR}/include")
LINK_DIRECTORIES("${OPENMM_DIR}/lib" "${OPENMM_DIR}/lib/plugins")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "{CMAKE_SOURCE_DIR}/cmake_modules")

# set flags for linking on mac
IF(APPLE)
    SET (CMAKE_INSTALL_NAME_DIR "@rpath")
    SET(EXTRA_COMPILE_FLAGS "-msse2 -stdlib=libc++")
ENDIF(APPLE)

# Select where to install
IF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
    IF(WIN32)
        SET(CMAKE_INSTALL_PREFIX "$ENV{ProgramFiles}/OpenMM" CACHE PATH "Where to install the plugin" FORCE)
    ELSE(WIN32)
        SET(CMAKE_INSTALL_PREFIX "/usr/local/openmm" CACHE PATH "Where to install the plugin" FORCE)
    ENDIF(WIN32)
ENDIF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(CAN_PLUGIN_SOURCE_SUBDIRS openmmapi serialization)

# Set the library name
SET(CAN_LIBRARY_NAME CustomAnisotropicNonbondedPlugin)
SET(SHARED_CAN_TARGET ${CAN_LIBRARY_NAME})

# These are all the places to search for header files which are to be part of the API.
SET(API_INCLUDE_DIRS "openmmapi/include" "openmmapi/include/internal")

# Locate header files.
SET(API_INCLUDE_FILES)
FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)
    SET(API_INCLUDE_FILES ${API_INCLUDE_FILES} ${fullpaths})
ENDFOREACH(dir)

# Collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)
FOREACH(subdir ${CAN_PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

FILE(GLOB src_files ${OPENMM_SOURCE}/libraries/asmjit/*/*.cpp)
FILE(GLOB incl_files ${OPENMM_SOURCE}/libraries/asmjit/*.h)
SET(SOURCE_FILES ${SOURCE_FILES} ${src_files})
SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
INCLUDE_DIRECTORIES(BEFORE "${OPENMM_SOURCE}/libraries/asmjit")
SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DLEPTON_USE_JIT")

# Create the library.

ADD_LIBRARY(${SHARED_CAN_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_CAN_TARGET}
    PROPERTIES COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DCAN_BUILDING_SHARED_LIBRARY -DLEPTON_BUILDING_SHARED_LIBRARY"
    LINK_FLAGS "${EXTRA_LINK_FLAGS}")

#link to libdl
FIND_LIBRARY(DL_LIBRARY dl)
IF(DL_LIBRARY)
  TARGET_LINK_LIBRARIES(${SHARED_CAN_TARGET} ${DL_LIBRARY} ${PTHREADS_LIB})
ENDIF(DL_LIBRARY)

INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${SHARED_CAN_TARGET})
TARGET_LINK_LIBRARIES(${SHARED_CAN_TARGET} OpenMM)

# install headers
FILE(GLOB API_ONLY_INCLUDE_FILES "openmmapi/include/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES} DESTINATION include)
FILE(GLOB API_ONLY_INCLUDE_FILES_INTERNAL "openmmapi/include/internal/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES_INTERNAL} DESTINATION include/internal)


# Enable testing

ENABLE_TESTING()
ADD_SUBDIRECTORY(serialization/tests)

# Build the implementations for different platforms

ADD_SUBDIRECTORY(platforms/reference)

FIND_PACKAGE(CUDA QUIET)
IF(CUDA_FOUND)
    SET(CAN_BUILD_CUDA_LIB ON CACHE BOOL "Build implementation for CUDA")
ELSE(CUDA_FOUND)
    SET(CAN_BUILD_CUDA_LIB OFF CACHE BOOL "Build implementation for CUDA")
ENDIF(CUDA_FOUND)
IF(CAN_BUILD_CUDA_LIB)
    ADD_SUBDIRECTORY(platforms/cuda)
ENDIF(CAN_BUILD_CUDA_LIB)

# Build the Python API
#

FIND_PROGRAM(PYTHON_EXECUTABLE python)
FIND_PROGRAM(SWIG_EXECUTABLE swig)
IF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
     SET(CAN_BUILD_PYTHON_WRAPPERS ON CACHE BOOL "Build wrappers for Python")
ELSE(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(CAN_BUILD_PYTHON_WRAPPERS OFF CACHE BOOL "Build wrappers for Python")
ENDIF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
IF(CAN_BUILD_PYTHON_WRAPPERS)
    ADD_SUBDIRECTORY(python)
ENDIF(CAN_BUILD_PYTHON_WRAPPERS)

