include(LibFindMacros)

# Finally the library itself
find_library(METEOIO_LIBRARY
  NAMES meteoio
  PATHS "/usr/lib" "/usr/local/lib" "~/usr/lib" "/opt/lib" ENV LD_LIBRARY_PATH
  DOC "Location of the libmeteoio, like /usr/lib"
 )

#build METEOIO_ROOT so we can provide a hint for searching for the header file
if ("${METEOIO_LIBRARY}" MATCHES "^(.+)lib[\\/]libmeteoio.(.+)$")
   set(METEOIO_ROOT "${CMAKE_MATCH_1}")
endif ("${METEOIO_LIBRARY}" MATCHES "^(.+)lib[\\/]libmeteoio.(.+)$")

# locate main header file
find_path(METEOIO_INCLUDE_DIR
  NAMES meteoio/MeteoIO.h
  HINTS ${METEOIO_ROOT}/include
  PATHS "/usr/include" "/usr/local/include" "~/usr/include" "/opt/include"
  DOC "Location of the meteoio headers, like /usr/include"
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(METEOIO_PROCESS_INCLUDES METEOIO_INCLUDE_DIR)
set(METEOIO_PROCESS_LIBS METEOIO_LIBRARY)
libfind_process(METEOIO)
