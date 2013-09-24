INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}")

IF(COMMAND cmake_policy)
	CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

IF(MSVC)
	SET(getopt_src getopt.c getopt_long.c)
ENDIF(MSVC)

SET(snowpack-app_sources
	Main.cc
	${getopt_src}
)

IF(POPC)
        SET(SN_LINKS "snowpackpopc")
ELSE(POPC)
        SET(SN_LINKS "snowpack")
ENDIF(POPC)
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/lib)

FIND_PACKAGE(MeteoIO REQUIRED)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${METEOIO_INCLUDE_DIR})
LINK_DIRECTORIES(${METEOIO_LIBRARIES}) #this, so that if the users types in paths, it works
LIST(APPEND SN_LINKS ${METEOIO_LIBRARIES})

IF(APPLE)
	#this is necessary for GUI exceptions
	LIST(APPEND SN_LINKS "-framework CoreServices")
ENDIF(APPLE)

#special features # currently set in the main CMakeLists.txt
#SET(DEBUG_ARITHM ON CACHE BOOL "Force-crash the application when doing an arithmetic exception")

IF(PLUGIN_IMISIO)
	FIND_PACKAGE(OCCI REQUIRED)
	INCLUDE_DIRECTORIES(SYSTEM ${OCCI_INCLUDE_DIR})
	LINK_DIRECTORIES(${OCCI_LIBRARIES})
	LIST(APPEND SN_LINKS ${OCCI_LIBRARIES})
ENDIF(PLUGIN_IMISIO)

IF(ENABLE_LAPACK)
	FIND_PACKAGE(Lapack REQUIRED)
	LINK_DIRECTORIES(${Lapack_LIBRARIES})
	LINK_DIRECTORIES(${Blas_LIBRARIES})
	LIST(APPEND SN_LINKS "lapack;blas")
ENDIF(ENABLE_LAPACK)

INCLUDE("${CMAKE_MODULE_PATH}/BuildVersion.cmake")
BuildVersion()

IF(UNIX)
	LIST(APPEND SN_LINKS "dl;pthread")
ENDIF(UNIX)

#Prepare executable
SET(BINARY "snowpack.app")
ADD_EXECUTABLE(${BINARY} ${snowpack-app_sources})

TARGET_LINK_LIBRARIES(${BINARY} ${SN_LINKS})

SET_TARGET_PROPERTIES(${BINARY} PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
	CLEAN_DIRECT_OUTPUT 1
	OUTPUT_NAME "snowpack"
)

INSTALL(TARGETS ${BINARY}
	RUNTIME DESTINATION bin
	COMPONENT exe
)

###########################################################
## Bundle MeteoIO
###########################################################
IF(WIN32)
	IF(CPACK_BUNDLE_ALL) #we consider that if !bundle, the user knows how to handle is library path
		FILE(GLOB METEOIO_DLL "${METEOIO_ROOT}/lib/libmeteo*.dll")
		FILE(GLOB SNOWPACK_DLL "${PROJECT_SOURCE_DIR}/lib/libsnow*.dll")
		INSTALL(FILES ${METEOIO_DLL} ${SNOWPACK_DLL} DESTINATION bin COMPONENT exe)
	ENDIF(CPACK_BUNDLE_ALL)
ENDIF(WIN32)
