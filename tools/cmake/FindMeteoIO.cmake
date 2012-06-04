INCLUDE(LibFindMacros)

# Finally the library itself
IF(WIN32)
	GET_FILENAME_COMPONENT(METEOIO_ROOT1 "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\Uninstall\\MeteoIO;UninstallString]" PATH CACHE INTERNAL)
	GET_FILENAME_COMPONENT(METEOIO_ROOT2 "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Wow6432Node\\Microsoft\\Windows\\CurrentVersion\\Uninstall\\MeteoIO;UninstallString]" PATH CACHE INTERNAL)
	GET_FILENAME_COMPONENT(METEOIO_ROOT3 "[HKEY_LOCAL_MACHINE\\SOFTWARE\\WSL Institute for Snow and Avalanche Research\\MeteoIO]" ABSOLUTE CACHE INTERNAL)
	#GET_FILENAME_COMPONENT(METEOIO_ROOT4 "C:/Progra~1/MeteoIO*" ABSOLUTE CACHE INTERNAL)
	SET(SEARCH_PATH
		ENV LIB
		${METEOIO_ROOT1}/lib
		${METEOIO_ROOT2}/lib
		${METEOIO_ROOT3}/lib
		"C:/Progra~1/MeteoIO*/lib" )

	IF(MSVC)
		FIND_LIBRARY(METEOIO_LIBRARY
			NAMES libmeteoio.lib
			PATHS ${SEARCH_PATH}
			DOC "Location of the libmeteoio, like c:/Program Files/MeteoIO-2.0.0/lib/libmeteoio.lib"
			)
	ELSE(MSVC)
		FIND_LIBRARY(METEOIO_LIBRARY
			NAMES libmeteoio.dll.a
			PATHS ${SEARCH_PATH}
			DOC "Location of the libmeteoio, like c:/Program Files/MeteoIO-2.0.0/lib/libmeteoio.dll.a"
			)
	ENDIF(MSVC)
ELSE(WIN32)
	IF(APPLE)
		FIND_LIBRARY(METEOIO_LIBRARY
		NAMES meteoio meteoiopopc
		PATHS
			"/Applications/MeteoIO/lib"
			ENV LD_LIBRARY_PATH
			ENV DYLD_FALLBACK_LIBRARY_PATH
			"~/usr/lib"
			"/usr/local/lib"
			"/usr/lib"
			"/opt/lib"
		DOC "Location of the libmeteoio, like /usr/lib/libmeteoio.dylib"
		)
	ELSE(APPLE)
		FIND_LIBRARY(METEOIO_LIBRARY
		NAMES meteoio meteoiopopc
		PATHS
			ENV LD_LIBRARY_PATH
			"~/usr/lib"
			"/usr/local/lib"
			"/usr/lib"
			"/opt/lib"
		DOC "Location of the libmeteoio, like /usr/lib/libmeteoio.so"
		)
	ENDIF(APPLE)
ENDIF(WIN32)

#build METEOIO_ROOT so we can provide a hint for searching for the header file
IF("${METEOIO_LIBRARY}" MATCHES "^(.+)lib[\\/]libmeteoio\\.(.+)$")
   SET(METEOIO_ROOT "${CMAKE_MATCH_1}")
ENDIF("${METEOIO_LIBRARY}" MATCHES "^(.+)lib[\\/]libmeteoio\\.(.+)$")

# locate main header file
FIND_PATH(METEOIO_INCLUDE_DIR
  NAMES meteoio/MeteoIO.h
  #HINTS ${METEOIO_ROOT}/include
  PATHS
	"${METEOIO_ROOT}/include"
	"~/usr/include"
	"/usr/local/include"
	"/usr/include"
	"/opt/include"
  DOC "Location of the meteoio headers, like /usr/include"
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET(METEOIO_PROCESS_INCLUDES METEOIO_INCLUDE_DIR)
SET(METEOIO_PROCESS_LIBS METEOIO_LIBRARY)
libfind_process(METEOIO)
