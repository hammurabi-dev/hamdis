CMAKE_MINIMUM_REQUIRED(VERSION 3.2.3)

#-------------- customized zone --------------#
SET(CMAKE_CXX_COMPILER "g++")
OPTION(ENABLE_TESTING "Enable testing for this project" ON)
OPTION(ENABLE_TIMING "Enable timing for this porect" ON)
OPTION(ON_DOCKER "Build on docker image" ON)
OPTION(BUILD_SHARED_LIB "Build shared library" ON)
OPTION(ENABLE_REPORT "Enable verbose report" ON)
#-------------- instruction ------------------#
# ENABLE_TESTING by default ON,
# ENABLE_TIMING is also not necessary for users,
# ON_DOCKER by defauly ON, we highly recommend Docker image for non-HPC tasks,
# but if install manually, switch it off and modify the path hints listed below,
# BUILD_SHARED_LIB by default ON,
# you will be overwhelmed by ENABLE_REPORT, switch it off for non-testing tasks,
# 
# you have to specify your local paths of external libraries just below here,
# in some special cases you have to modify FIND_PATH/FIND_LIBRARY functions,
#
# in some special cases, you may want to manually fix LFLAGS or CFLAGS,
#
# if you add new modules/derived classes beyond original code,
# please manually add source file paths to SET(SRC_FILES ...) function,
#
# @ enthusiastic or paranoid developers:
# we use Google Test for tests,
# Google Test package is assembled INTO testing modules manually,
# you can either install GoogleTest and cp src dir into install path,
# or just download GoogleTest and specify root dir to GTEST_HINTS,
#--------------------------------------------#
IF(ON_DOCKER)
	SET(GTEST_HINTS /usr/local)
	SET(INSTALL_ROOT_DIR /usr/local/hamdis)
ELSE()
	MESSAGE("apply customized library paths")
	SET(GTEST_HINTS $HOME/local/googletest)
ENDIF()
SET(INSTALL_ROOT_DIR $HOME/local/hamdis)
SET(CMAKE_EXE_LINKER_FLAGS -lcurl)
#---------------------------------------------#

PROJECT(hamdis CXX)

# SETup cxx standard
SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_STANDARD_REQUIRED YES)

# compile flags under defferent options
IF(ENABLE_TESTING AND ENABLE_TIMING)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated")
ELSEIF(ENABLE_TIMING)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -DNDEBUG")
ELSEIF(ENABLE_TESTING)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -DNTIMING")
ELSE()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -DNTIMING -DNDEBUG")
ENDIF()

IF (ENABLE_REPORT)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")
ENDIF()

# openmp and thread support
# if FindOpenMP fails, try add -fopenmp to CMAKE_CXX_FLAGS above
# the same solution applies to -pthread
include(FindOpenMP)
IF(OPENMP_FOUND)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ELSE()
	MESSAGE(FATAL_ERROR "openmp unsupported")
ENDIF()
IF(ENABLE_TESTING)
	FIND_PACKAGE(Threads REQUIRED)
	LIST(APPEND ALL_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
ENDIF()

# we assemble include and external libs together
SET(ALL_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
SET(ALL_LIBRARIES)

# find GoogleTest
# GoogleTest package must be assembled INTO testing modules manually
# you can either install GoogleTest and cp src dir into install path
# or just download GoogleTest
IF(ENABLE_TESTING)
	FIND_PATH(GTEST_INC_DIR 
		NAMES gtest/gtest.h
  		HINTS
		${GTEST_HINTS}/include
		${GTEST_HINTS}/googletest/include
	)
	IF(NOT GTEST_INC_DIR)
  		MESSAGE(FATAL_ERROR ${GTEST_INC_DIR})
	ENDIF()
	FIND_PATH(GTEST_SOURCE_DIR gtest-all.cc
  		${GTEST_HINTS}/src
  		${GTEST_HINTS}/googletest/src
	)
	IF(NOT GTEST_SOURCE_DIR)
  		MESSAGE(FATAL_ERROR ${GTEST_SOURCE_DIR})
	ENDIF()
	SET(GTEST_LIB_SOURCES ${GTEST_SOURCE_DIR}/gtest-all.cc)
	SET(GTEST_MAIN_SOURCES ${GTEST_SOURCE_DIR}/gtest_main.cc)
	SET(GTEST_INCLUDE_DIRS ${GTEST_INC_DIR} ${GTEST_SOURCE_DIR} ${GTEST_SOURCE_DIR}/..)
ENDIF()

# assemble lib
INCLUDE_DIRECTORIES(${ALL_INCLUDE_DIR})
ADD_LIBRARY(hamdis INTERFACE)

# build testing cases
IF(ENABLE_TESTING)
	ENABLE_TESTING()
	ADD_SUBDIRECTORY(tests)
ENDIF()
