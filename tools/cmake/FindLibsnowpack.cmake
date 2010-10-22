FIND_PATH(LIBSNOWPACK_INCLUDE_DIR snowpack/libsnowpack.h "/usr/include" "/usr/local/include" "~/usr/include" "/opt/include")
FIND_LIBRARY(LIBSNOWPACK_SHARED libsnowpack PATHS "/usr/lib/" "/usr/local/lib/" "~/usr/lib/" "/opt/lib")

IF (LIBSNOWPACK_INCLUDE_DIR AND LIBSNOWPACK_SHARED)
   SET(LIBSNOWPACK_FOUND TRUE)
ENDIF (LIBSNOWPACK_INCLUDE_DIR AND LIBSNOWPACK_SHARED)

IF (LIBSNOWPACK_FOUND)
   IF (NOT Libsnowpack_FIND_QUIETLY)
      MESSAGE(STATUS "Found libsnowpackO: ${LIBSNOWPACK_SHARED}")
   ENDIF (NOT Libsnowpack_FIND_QUIETLY)
ELSE (LIBSNOWPACK_FOUND)
   IF (Libsnowpack_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find libsnowpack")
   ENDIF (Libsnowpack_FIND_REQUIRED)
ENDIF (LIBSNOWPACK_FOUND)
