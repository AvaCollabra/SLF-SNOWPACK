include(LibFindMacros)

# Finally the library itself
find_library(LIBSNOWPACK_LIBRARY
  NAMES snowpack
  PATHS "/usr/lib" "/usr/local/lib" "~/usr/lib" "/opt/lib" ENV LD_LIBRARY_PATH
  DOC "Location of the libsnowpack, like /usr/lib"
)

#build LIBSNOWPACK_ROOT so we can provide a hint for searching for the header file
if ("${LIBSNOWPACK_LIBRARY}" MATCHES "^(.+)lib[\\/]libsnowpack.(.+)$")
   set(LIBSNOWPACK_ROOT "${CMAKE_MATCH_1}")
endif ("${LIBSNOWPACK_LIBRARY}" MATCHES "^(.+)lib[\\/]libsnowpack.(.+)$")

# locate main header file
find_path(LIBSNOWPACK_INCLUDE_DIR
  NAMES snowpack/libsnowpack.h
  HINTS ${LIBSNOWPACK_ROOT}/include
  PATHS "/usr/include" "/usr/local/include" "~/usr/include" "/opt/include"
  DOC "Location of the libsnowpack headers, like /usr/include"
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(LIBSNOWPACK_PROCESS_INCLUDES LIBSNOWPACK_INCLUDE_DIR)
set(LIBSNOWPACK_PROCESS_LIBS LIBSNOWPACK_LIBRARY)
libfind_process(LIBSNOWPACK)
