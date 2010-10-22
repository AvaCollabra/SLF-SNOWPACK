include(LibFindMacros)

# Glib-related libraries also use a separate config header, which is in lib dir
find_path(METEOIO_INCLUDE_DIR
  NAMES meteoio/MeteoIO.h
  PATHS "/usr/include" "/usr/local/include" "~/usr/include" "/opt/include"
)

# Finally the library itself
find_library(METEOIO_LIBRARY
  NAMES meteoio
  PATHS "/usr/lib" "/usr/local/lib" "~/usr/lib" "/opt/lib"
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(METEOIO_PROCESS_INCLUDES METEOIO_INCLUDE_DIR)
set(METEOIO_PROCESS_LIBS METEOIO_LIBRARY)
libfind_process(METEOIO)
