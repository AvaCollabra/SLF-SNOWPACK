cmake_minimum_required(VERSION 2.4) 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tools/cmake/")

# Version of library
SET(VERSION_MAJOR 1 CACHE INTERNAL "")
SET(VERSION_MINOR 0 CACHE INTERNAL "")
SET(VERSION_PATCH 0 CACHE INTERNAL "")

PROJECT(snowpack)
ADD_SUBDIRECTORY(snowpack)

# require proper c++
set(DEBUG_ARITHM ON CACHE BOOL "Use the debug flags -DDEBUG_ARITHM to catch arithmetic exceptions")
set(WARNINGS "-Wall")
#set(EXTRA_WARNINGS "-Wextra -ansi -pedantic")

IF(DEBUG_ARITHM)
  LIST(APPEND CFLAGS "-DDEBUG_ARITHM")
ENDIF(DEBUG_ARITHM)

###########################################################
#platform specific
IF(UNIX)
  IF(APPLE)
    set(PLATFORM "-DAPPLE")
  ELSE(APPLE)
    set(PLATFORM "-DLINUX")
  ENDIF(APPLE)
ELSE(UNIX)
  IF(WIN32)
    set(PLATFORM "-DWIN32 -U__STRICT_ANSI__")
  ENDIF(WIN32)
ENDIF(UNIX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "debug" CACHE STRING
      "Choose the type of build, options are: debug release "
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(PLUGIN_IMISIO OFF CACHE BOOL "Compilation IMISDBIO ON or OFF")
IF(PLUGIN_IMISIO)
  SET(MYDEFINES "-DIMISDBIO")
ENDIF(PLUGIN_IMISIO)

###########################################################
#targets providing sets of compiler options
if(NOT DEST)
  set(DEST "safe" CACHE STRING "Choose safe or optimized" FORCE)
endif(NOT DEST)

if (DEST STREQUAL "safe")
  set(OPTIM "-g -O3 -DNDEBUG")
  set(ARCH  "")
  set(DEBUG "-g -O0 -D__DEBUG")
endif(DEST STREQUAL "safe")

IF(DEST STREQUAL "optimized")
  SET(OPTIM "-g -O3 -DNDEBUG")
  SET(ARCH  "-mtune=native")
  SET(DEBUG "-g -O0 -D__DEBUG")
ENDIF(DEST STREQUAL "optimized")

set(CMAKE_CXX_FLAGS_RELEASE "${ARCH} ${OPTIM} ${CFLAGS} ${PLATFORM} ${MYDEFINES}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "${DEBUG} ${ARCH} ${CFLAGS} ${PLATFORM} ${WARNINGS} ${EXTRA_WARNINGS} ${MYDEFINES}" CACHE STRING "" FORCE)

message(STATUS "DEBUG CXXFLAGS  : ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "RELEASE CXXFLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

###########################################################
set(LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

###########################################################
#for the install target
FILE(GLOB hfiles "${CMAKE_CURRENT_SOURCE_DIR}/snowpack/*.h")
INSTALL(FILES ${hfiles} DESTINATION include/snowpack)

FILE(GLOB solibs "${PROJECT_SOURCE_DIR}/lib/libsnowpack*.so*")
FILE(GLOB alibs "${PROJECT_SOURCE_DIR}/lib/*.a*")

INSTALL(FILES ${solibs} ${alibs} DESTINATION lib)

add_custom_target(distclean make clean
  COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake meteoio/cmake_install.cmake meteoio/plugins/cmake_install.cmake meteoio/plugins/gsn/cmake_install.cmake
  COMMAND cmake -E remove_directory CMakeFiles
  COMMAND cmake -E remove_directory snowpack/CMakeFiles
)
###########################################################
#for the uninstall target
CONFIGURE_FILE(
  "${PROJECT_SOURCE_DIR}/tools/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

###########################################################
#creation of packages NOT FULLY WORKING YET
set(CPACK_PACKAGE_DESCRIPTION "Snowpack Library")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "..")
set(CPACK_PACKAGE_NAME "libsnowpack")
set(CPACK_PACKAGE_CONTACT "Mathias Bavay ")
set(CPACK_PACKAGE_VENDOR "WSL Institute for Snow and Avalanche Research")
set(CPACK_DEBIAN_PACKAGE_NAME "libsnowpack")
set(CPACK_DEBIAN_PACKAGE_SECTION "devel")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
set(CPACK_RPM_PACKAGE_LICENSE "LPGLv3")
set(CPACK_RPM_PACKAGE_REQUIRES "libstdc++6")
set(CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
set(CPACK_RPM_PACKAGE_GROUP "Development/Libraries")

set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

set(CPACK_GENERATOR "DEB;RPM;TGZ;")
#set(CPACK_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_ARCH "i386")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}${PACKAGE_EXTRA}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}-${CPACK_ARCH}")

include(CPack)
