#############################
## Global Configuration    ##
#############################

# set Memory check command
set(CTEST_MEMORYCHECK_COMMAND valgrind)

#get the proper Snowpack library
IF(BUILD_SHARED_LIBS)
	SET(LIBSNOWPACK_LIBRARY ${PROJECT_NAME})
ELSE(BUILD_SHARED_LIBS)
	IF(BUILD_STATIC_LIBS)
		SET(LIBSNOWPACK_LIBRARY "${PROJECT_NAME}_STATIC")
	ELSE(BUILD_STATIC_LIBS)
		MESSAGE(SEND_ERROR "Not building Snowpack, the test won't be able to build")
	ENDIF(BUILD_STATIC_LIBS)
ENDIF(BUILD_SHARED_LIBS)
INCLUDE_DIRECTORIES(../)

FIND_PACKAGE(MeteoIO REQUIRED)
INCLUDE_DIRECTORIES(${METEOIO_INCLUDE_DIR})

IF(APPLE)
	#this is necessary for GUI exceptions
	SET(EXTRA_LINKS "-framework CoreServices")
ENDIF(APPLE)
IF(UNIX)
	SET(EXTRA_LINKS "dl;pthread")
ENDIF(UNIX)

#in order to link a new executable: 
#TARGET_LINK_LIBRARIES(${BINARY} ${LIBSNOWPACK_LIBRARY} ${METEOIO_LIBRARY} ${EXTRA_LINKS})

###################
## Add Tests     ##
###################
ADD_SUBDIRECTORY(res1exp)
ADD_SUBDIRECTORY(res5exp)
ADD_SUBDIRECTORY(basics)
ADD_SUBDIRECTORY(mass_and_energy_balance)
