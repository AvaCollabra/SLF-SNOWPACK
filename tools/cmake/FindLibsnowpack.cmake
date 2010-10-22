include(LibFindMacros)

# Glib-related libraries also use a separate config header, which is in lib dir
find_path(LIBSNOWPACK_INCLUDE_DIR
  NAMES snowpack/libsnowpack.h
  PATHS "/usr/include" "/usr/local/include" "~/usr/include" "/opt/include"
)

# Finally the library itself
find_library(LIBSNOWPACK_LIBRARY
  NAMES snowpack
  PATHS "/usr/lib" "/usr/local/lib" "~/usr/lib" "/opt/lib"
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(LIBSNOWPACK_PROCESS_INCLUDES LIBSNOWPACK_INCLUDE_DIR)
set(LIBSNOWPACK_PROCESS_LIBS LIBSNOWPACK_LIBRARY)
libfind_process(LIBSNOWPACK)
