INCLUDE(LibFindMacros)

# Finally the library itself
IF(WIN32)
	GET_FILENAME_COMPONENT(LIBSNOWPACK_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\WSL Institute for Snow and Avalanche Research\\Snowpack]" ABSOLUTE CACHE)

	FIND_LIBRARY(LIBSNOWPACK_LIBRARY
		NAMES libsnowpack.lib
		PATHS
			ENV LIB
			${LIBSNOWPACK_ROOT}/lib
			"C:/Program Files/Snowpack/lib"
		DOC "Location of the libsnowpack, like c:/Program Files/Snowpack-2.0.0/lib/libsnowpack.lib"
		)
ELSE(WIN32)
	FIND_LIBRARY(LIBSNOWPACK_LIBRARY
	NAMES snowpack
	PATHS
		ENV LD_LIBRARY_PATH
		"~/usr/lib"
		"/usr/local/lib"
		"/usr/lib"
		"/opt/lib"
	DOC "Location of the libsnowpack, like /usr/lib"
	)
END(WIN32)

#build LIBSNOWPACK_ROOT so we can provide a hint for searching for the header file
IF("${LIBSNOWPACK_LIBRARY}" MATCHES "^(.+)lib[\\/]libsnowpack\\.(.+)$")
   SET(LIBSNOWPACK_ROOT "${CMAKE_MATCH_1}")
ENDIF("${LIBSNOWPACK_LIBRARY}" MATCHES "^(.+)lib[\\/]libsnowpack\\.(.+)$")

# locate main header file
FIND_PATH(LIBSNOWPACK_INCLUDE_DIR
  NAMES snowpack/libsnowpack.h
  #HINTS ${LIBSNOWPACK_ROOT}/include
  PATHS
	"${LIBSNOWPACK_ROOT}/include"
	"~/usr/include"
	"/usr/local/include"
	"/usr/include"
	"/opt/include"
  DOC "Location of the libsnowpack headers, like /usr/include"
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET(LIBSNOWPACK_PROCESS_INCLUDES LIBSNOWPACK_INCLUDE_DIR)
SET(LIBSNOWPACK_PROCESS_LIBS LIBSNOWPACK_LIBRARY)
libfind_process(LIBSNOWPACK)
