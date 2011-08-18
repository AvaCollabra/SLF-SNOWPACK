cmake_minimum_required(VERSION 2.4)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tools/cmake/")

# Version of library
SET(VERSION_MAJOR 2 CACHE INTERNAL "")
SET(VERSION_MINOR 0 CACHE INTERNAL "")
SET(VERSION_PATCH 0 CACHE INTERNAL "")

SET(LIBPREFIX "lib") #make sure ALL plugins are libXXX for ALL plateforms
SET(CMAKE_IMPORT_LIBRARY_PREFIX "lib") #for Visual C++ import libraries

PROJECT(snowpack)
ADD_SUBDIRECTORY(snowpack)

###########################################################
#plugins
set(PLUGIN_IMISIO OFF CACHE BOOL "Compilation IMISDBIO ON or OFF - only relevant for SLF")

###########################################################
#platform specific. First, compiler specific flags, then system specific flags
#SET(DEBUG_ARITHM ON CACHE BOOL "Force crash the application when doing an arithmetic exception")
OPTION( BUILD_SHARED_LIBS "Set to OFF to build only static libraries" ON )
INCLUDE("${CMAKE_MODULE_PATH}/BuildVersion.cmake")
BuildVersion()

IF(CMAKE_COMPILER_IS_GNUCXX)
	IF(PLUGIN_IMISIO)
		LIST(APPEND CFLAGS " -D${IMISDBIO}")
	ENDIF(PLUGIN_IMISIO)
	SET(WARNINGS "-Wall")
	SET(PROFILING "-pg -fprofile-arcs")
	SET(EXTRA_WARNINGS "-Wextra -ansi -pedantic")
	IF(WIN32)
		LIST(APPEND CFLAGS " -D_USE_MATH_DEFINES") #USE_MATH_DEFINES needed for Win32
	ENDIF(WIN32)
	SET(OPTIM "-g -O3 -DNDEBUG -DNOSAFECHECKS")
	EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	IF (GCC_VERSION VERSION_GREATER 4.2 OR GCC_VERSION VERSION_EQUAL 4.2)
		SET(ARCH_OPTIM  "-march=native -mtune=native")
	ENDIF()
	SET(ARCH_SAFE  "")
	SET(DEBUG "-g3 -O0 -D__DEBUG")
	SET(_VERSION "-D_VERSION=\\\"${_versionString}\\\"")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

IF(MSVC)
	IF(PLUGIN_IMISIO)
		LIST(APPEND CFLAGS " -/D${IMISDBIO}")
	ENDIF(PLUGIN_IMISIO)
	#SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited configs"  FORCE)
	SET(WARNINGS "/W4")
	SET(EXTRA_WARNINGS "/Wp64") #/Wall
	SET(OPTIM "/O2 /Vmb /DNDEBUG /MD /DNOSAFECHECKS")
	SET(ARCH_OPTIM "/arch:SSE2")
	SET(ARCH_SAFE "")
	SET(DEBUG "/Z7 /Od /D__DEBUG /MDd")
	LIST(APPEND CFLAGS " /D_CRT_SECURE_NO_WARNINGS /D_USE_MATH_DEFINES") #Za: strict ansi
	SET(_VERSION "/D_VERSION=\\\"${_versionString}\\\"")
ENDIF(MSVC)

IF(UNIX)
	IF(APPLE)
		SET(PLATFORM "-DAPPLE")
		SET(DLL_EXT "dylib")
		SET(STAT_EXT "a")
	ELSE(APPLE)
		SET(PLATFORM "-DLINUX")
		SET(DLL_EXT "so")
		SET(STAT_EXT "a")
	ENDIF(APPLE)
ELSE(UNIX)
	IF(WIN32)
		SET(PLATFORM "-D_WIN32 -U__STRICT_ANSI__")
		SET(DLL_EXT "dll")
		SET(STAT_EXT "a")
	ENDIF(WIN32)
ENDIF(UNIX)

IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE "release" CACHE STRING "Choose the type of build, options are: debug release " FORCE)
	SET(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE})
ENDIF(NOT CMAKE_BUILD_TYPE)

###########################################################
#targets providing sets of compiler options
if(NOT DEST)
  set(DEST "safe" CACHE STRING "Choose safe or optimized" FORCE)
endif(NOT DEST)

IF (DEST STREQUAL "safe")
	SET(ARCH  "${ARCH_SAFE}")
ENDIF(DEST STREQUAL "safe")

IF(DEST STREQUAL "optimized")
	SET(ARCH  "${ARCH_OPTIM}")
ENDIF(DEST STREQUAL "optimized")

###########################################################
#finally, SET compile flags
SET(CMAKE_CXX_FLAGS_RELEASE "${_VERSION} ${ARCH} ${OPTIM} ${CFLAGS} ${PLATFORM}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_DEBUG "${DEBUG} ${_VERSION} ${ARCH} ${CFLAGS} ${PLATFORM} ${WARNINGS} ${EXTRA_WARNINGS}" CACHE STRING "" FORCE)

MESSAGE(STATUS "DEBUG CXXFLAGS  : ${CMAKE_CXX_FLAGS_DEBUG}")
MESSAGE(STATUS "RELEASE CXXFLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

###########################################################
#for the install target
FILE(GLOB hfiles "snowpack/*.h")
INSTALL(FILES ${hfiles} DESTINATION include/snowpack COMPONENT headers)

IF(WIN32) #handle the special case of WIN32 needing import libraries
	FILE(GLOB_RECURSE importlibs "lib/libsnowpack*.lib*")
ENDIF(WIN32)
FILE(GLOB_RECURSE solibs "lib/libsnowpack*.${DLL_EXT}*")
FILE(GLOB_RECURSE alibs "lib/*.${STAT_EXT}*")

INSTALL(FILES ${solibs} ${importlibs} ${alibs} DESTINATION lib COMPONENT libraries)

ADD_CUSTOM_TARGET(distclean make clean
	COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake snowpack/cmake_install.cmake
	COMMAND cmake -E remove_directory CMakeFiles
	COMMAND cmake -E remove_directory snowpack/CMakeFiles)

###########################################################
#for the uninstall target
CONFIGURE_FILE(
	"${PROJECT_SOURCE_DIR}/tools/cmake/cmake_uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

#Build the documentation if doxygen is found
include(FindDoxygen)
FIND_PACKAGE(Doxygen)
if (DOXYGEN_FOUND)
	ADD_CUSTOM_TARGET(documentation "doxygen" ${PROJECT_SOURCE_DIR}/snowpack/config.dox)
	get_target_property(DOC_TARGET doc TYPE)
	if(NOT DOC_TARGET)
		add_custom_target(doc)
	endif(NOT DOC_TARGET)
	add_dependencies(doc documentation)
endif (DOXYGEN_FOUND)
#Install the documentation if the user desires it
SET(INSTALL_DOC OFF CACHE BOOL "Install documentation ON or OFF")
IF(INSTALL_DOC)
	IF(WIN32)
		SET(DOC_PATH "doc")
	ELSE(WIN32)
		SET(DOC_PATH "share/doc/${CMAKE_PROJECT_NAME}")
	ENDIF(WIN32)
	SET(DOC_FILES License.txt snowpack/icons/icon.ico)
	IF (DOXYGEN_FOUND)
		INSTALL(DIRECTORY doc/html DESTINATION ${DOC_PATH} COMPONENT documentation)
		LIST(APPEND DOC_FILES doc/Snowpack_doc.html)
	ELSE (DOXYGEN_FOUND)
		MESSAGE("Doxygen has not been found and is needed to build online documentation. Please either install doxygen or download a documentation-only package!")
	ENDIF (DOXYGEN_FOUND)

	INSTALL(FILES ${DOC_FILES} DESTINATION ${DOC_PATH} COMPONENT documentation)
	INSTALL(DIRECTORY doc/examples DESTINATION ${DOC_PATH} COMPONENT examples PATTERN ".svn" EXCLUDE)
ENDIF(INSTALL_DOC)

###########################################################
## CPACK section
###########################################################
set(CPACK_PACKAGE_DESCRIPTION "Snowpack Library")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "This is the SNOWPACK model ([Bartelt 2002], [Lehning 2002a], [Lehning 2002b]) for physical modeling of snow pack as a library. It can therefore be integrated into other models.")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/Welcome.txt")
set(CPACK_PACKAGE_NAME "libsnowpack")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Snowpack")
set(CPACK_PACKAGE_CONTACT "Mathias Bavay ")
set(CPACK_PACKAGE_VENDOR "WSL Institute for Snow and Avalanche Research")
SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/snowpack/icons\\\\icon.ico")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/License.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/Readme.txt")

set(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")
set(VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
#SET(CPACK_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
SET(CPACK_ARCH "i386")
SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}${PACKAGE_EXTRA}-${VERSION}-${CPACK_ARCH}")

#installation by components for formats that support it
#SET(CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE 1)
SET(CPACK_COMPONENTS_ALL libraries headers documentation examples)
SET(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "Snowpack library")
SET(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "The SNOWPACK library itself")
SET(CPACK_COMPONENT_LIBRARIES_GROUP Runtime)
SET(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "Header files")
SET(CPACK_COMPONENT_HEADERS_DESCRIPTION "Necessary C++ header files for developing applications using SNOWPACK")
SET(CPACK_COMPONENT_HEADERS_GROUP Development)
SET(CPACK_COMPONENT_HEADERS_DEPENDS libraries)
SET(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "Documentation")
SET(CPACK_COMPONENT_DOCUMENTATION_DESCRIPTION "Documentation of the SNOWPACK usage and API both for the end user and the developer using SNOWPACK in his/her application")
SET(CPACK_COMPONENT_DOCUMENTATION_GROUP Documentation)
SET(CPACK_COMPONENT_EXAMPLES_DISPLAY_NAME "Examples")
SET(CPACK_COMPONENT_EXAMPLES_DESCRIPTION "Various simulations with their example data")
SET(CPACK_COMPONENT_EXAMPLES_GROUP Documentation)
SET(CPACK_COMPONENT_EXAMPLES_DEPENDS libraries headers runtime)

SET(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")
SET(CPACK_COMPONENT_GROUP_RUNTIME_DESCRIPTION "What you need to run an application based on SNOWPACK")
SET(CPACK_COMPONENT_DEVELOPMENT_GROUP "Development")
SET(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "What you need to develop an application based on SNOWPACK")

IF(UNIX)
	IF(APPLE)
		set(CPACK_GENERATOR "PackageMaker;")
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/Applications/${CPACK_PACKAGE_NAME}-${VERSION}")
	ELSE(APPLE)
		set(CPACK_GENERATOR "DEB;RPM;TGZ")
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
		set(CPACK_DEBIAN_PACKAGE_NAME "libsnowpack")
		set(CPACK_DEBIAN_PACKAGE_SECTION "devel")
		set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6, meteoio") #objdump -p lib/libmeteoio.so | grep NEEDED | tr -s ' ' | cut -d' ' -f3 | tr '\n' ','
		set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386") #dpkg --print-architecture
		set(CPACK_RPM_PACKAGE_NAME "libsnowpack")
		set(CPACK_RPM_PACKAGE_LICENSE "LPGLv3")
		set(CPACK_RPM_PACKAGE_REQUIRES "libstdc++6, meteoio")
		set(CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
		set(CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
	ENDIF(APPLE)
ELSE(UNIX)
	IF(WIN32)
		SET(CPACK_GENERATOR "NSIS")
		SET(CPACK_NSIS_MUI_ICON "${CPACK_PACKAGE_ICON}")
		SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}-${VERSION}")
		SET(CPACK_NSIS_MENU_LINKS
			"${DOC_PATH}/Snowpack_doc.html" "SNOWPACK documentation"
			"${DOC_PATH}/examples" "SNOWPACK examples"
			"https://slfsmm.indefero.net/p/snowpack/" "Home page")
		#info for the add/uninstall software administration
		SET(CPACK_NSIS_INSTALLED_ICON_NAME "${DOC_PATH}\\\\icon.ico")
		SET(CPACK_NSIS_URL_INFO_ABOUT "http://www.slf.ch")
		SET(CPACK_NSIS_HELP_LINK "https://slfsmm.indefero.net/p/snowpack")
		SET(CPACK_NSIS_CONTACT "Mathias.Bavay@wsl.ch")
	ENDIF(WIN32)
ENDIF(UNIX)


set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${VERSION}" CACHE INTERNAL "tarball basename")
set(CPACK_SOURCE_GENERATOR TGZ)
# The following components are regex's to match anywhere (unless anchored)
# in absolute path + filename to find files or directories to be excluded
# from source tarball.
set(CPACK_SOURCE_IGNORE_FILES
	"~$"
	"\\\\.svn/"
	"/CMakeFiles/"
	"/CMakeCache*"
	"/cmake_install\\\\.cmake$"
	"/cmake_uninstall\\\\.cmake$"
	"^_CPack_.*/"
	"/CPack*"
	"\\\\.o$"
	"/m4/"
)

include(CPack)
