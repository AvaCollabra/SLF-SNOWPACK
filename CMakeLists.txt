cmake_minimum_required(VERSION 2.4)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tools/cmake/")

# Version of library
SET(VERSION_MAJOR 3 CACHE INTERNAL "")
SET(VERSION_MINOR 0 CACHE INTERNAL "")
SET(VERSION_PATCH 0 CACHE INTERNAL "")

SET(LIBPREFIX "lib") #make sure ALL plugins are libXXX for ALL plateforms
SET(CMAKE_IMPORT_LIBRARY_PREFIX "lib") #for Visual C++ import libraries

PROJECT(snowpack)
ADD_SUBDIRECTORY(snowpack)

IF(SNOWPACK_STANDALONE)
	#it seems that the only reliable way is to test-compile...
	INCLUDE (CheckCXXSourceCompiles)
	CHECK_CXX_SOURCE_COMPILES("
	#include <fenv.h>
	int main() {
	feenableexcept(FE_DIVBYZERO);
	return 0;
	}" HAVE_FEENABLE)

	IF(HAVE_FEENABLE)
		SET(DEBUG_ARITHM ON CACHE BOOL "Force-crash the application when doing an arithmetic exception")
	ENDIF(HAVE_FEENABLE)
	ADD_SUBDIRECTORY(applications/snowpack)
#ELSE(SNOWPACK_STANDALONE)
#	SET(DEBUG_ARITHM OFF CACHE BOOL "Force-crash the application when doing an arithmetic exception")
ENDIF(SNOWPACK_STANDALONE)

###########################################################
#plugins and applications
SET(PLUGIN_IMISIO OFF CACHE BOOL "Compilation IMISDBIO ON or OFF - only relevant for SLF")
SET(SNOWPACK_STANDALONE ON CACHE BOOL "Compilation SNOWPACK standalone ON or OFF")

###########################################################
#platform specific. First, compiler specific flags, then system specific flags
OPTION( BUILD_SHARED_LIBS "Set to ON to build dynamic libraries" ON )
INCLUDE("${CMAKE_MODULE_PATH}/BuildVersion.cmake")
BuildVersion()

IF(MSVC)
	IF(PLUGIN_IMISIO)
		SET(EXTRA "${EXTRA} -/DIMISDBIO")
	ENDIF(PLUGIN_IMISIO)
	#SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited configs"  FORCE)
	SET(WARNINGS "/W4")
	#SET(EXTRA_WARNINGS "/Wp64") #/Wall
	SET(OPTIM "/O2 /DNDEBUG /MD /DNOSAFECHECKS")
	SET(ARCH_OPTIM "/arch:SSE2")
	SET(ARCH_SAFE "")
	SET(DEBUG "/Z7 /Od /D__DEBUG /MDd")
	IF(DEBUG_ARITHM)
		SET(EXTRA "${EXTRA} /EHa")
	ENDIF(DEBUG_ARITHM)
	LIST(APPEND CFLAGS " /D_CRT_SECURE_NO_WARNINGS /D_USE_MATH_DEFINES /EHsc") #Za: strict ansi
	SET(_VERSION "/D_VERSION=\\\"${_versionString}\\\"")
ELSE(MSVC)
	#we consider that all other compilers support "-" options and silently ignore what they don't know
	IF(PLUGIN_IMISIO)
		SET(EXTRA "${EXTRA} -DIMISDBIO")
	ENDIF(PLUGIN_IMISIO)
	SET(WARNINGS "-Wall")
	SET(PROFILING "-pg -fprofile-arcs")
	SET(EXTRA_WARNINGS "-Wextra -pedantic")
	IF(CMAKE_COMPILER_IS_GNUCXX)
		SET(EXTRA_WARNINGS "${EXTRA_WARNINGS} -ansi")
		IF(WIN32) #for gcc on windows
			SET(CMAKE_SHARED_LINKER_FLAGS "--enable-auto-import")
		ENDIF(WIN32)
	ENDIF(CMAKE_COMPILER_IS_GNUCXX)
	IF(WIN32)
		LIST(APPEND CFLAGS " -D_USE_MATH_DEFINES") #USE_MATH_DEFINES needed for Win32
	ENDIF(WIN32)
	SET(OPTIM "-g -O3 -DNDEBUG -DNOSAFECHECKS")
	EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	IF(CMAKE_COMPILER_IS_GNUCXX AND (GCC_VERSION VERSION_GREATER 4.2 OR GCC_VERSION VERSION_EQUAL 4.2) )
		SET(ARCH_OPTIM  "-march=native -mtune=native")
	ENDIF()
	SET(ARCH_SAFE  "")
	SET(DEBUG "-g3 -O0 -D__DEBUG")
	IF(DEBUG_ARITHM)
		SET(EXTRA "${EXTRA} -DDEBUG_ARITHM")
	ENDIF(DEBUG_ARITHM)
	SET(_VERSION "-D_VERSION=\\\"${_versionString}\\\"")
ENDIF(MSVC)

IF(UNIX)
	IF(APPLE)
		SET(PLATFORM "-DAPPLE")
		SET(DLL_EXT "dylib")
		SET(STAT_EXT "a")
	ELSE(APPLE)
		SET(PLATFORM "-DLINUX")
		SET(DLL_EXT "so*")
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
SET(CMAKE_CXX_FLAGS_RELEASE "${_VERSION} ${EXTRA} ${ARCH} ${OPTIM} ${CFLAGS} ${PLATFORM}" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_DEBUG "${DEBUG} ${_VERSION} ${EXTRA} ${ARCH} ${CFLAGS} ${PLATFORM} ${WARNINGS} ${EXTRA_WARNINGS}" CACHE STRING "" FORCE)

MESSAGE(STATUS "DEBUG CXXFLAGS  : ${CMAKE_CXX_FLAGS_DEBUG}")
MESSAGE(STATUS "RELEASE CXXFLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

###########################################################
#for the install target
FILE(GLOB hfiles "snowpack/*.h")

IF(WIN32) #handle the special case of WIN32 needing import libraries
	FILE(GLOB_RECURSE importlibs "lib/libsnow*.lib*")
ENDIF(WIN32)
FILE(GLOB_RECURSE solibs "lib/libsnow*.${DLL_EXT}")
FILE(GLOB_RECURSE alibs "lib/*.${STAT_EXT}")

INSTALL(FILES ${hfiles} DESTINATION include/snowpack COMPONENT headers)
INSTALL(FILES ${solibs} ${importlibs} ${alibs} DESTINATION lib COMPONENT libraries)

# DAVID: WHY twice ??? command should work even if file is not present !
IF(SNOWPACK_STANDALONE)
	ADD_CUSTOM_TARGET(distclean make clean
		COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake snowpack/cmake_install.cmake applications/snowpack/cmake_install.cmake tests/*.cmake
		COMMAND cmake -E remove_directory CMakeFiles
		COMMAND cmake -E remove_directory snowpack/CMakeFiles
		COMMAND cmake -E remove_directory applications/snowpack/CMakeFiles
		COMMAND cmake -E remove_directory tests/CMakeFiles
		COMMAND cmake -E remove tests/res1exp/*.cmake tests/res1exp/output/*
		COMMAND cmake -E remove_directory tests/res1exp/CMakeFiles
		COMMAND cmake -E remove_directory tests/res1exp/Testing
		)
ELSE(SNOWPACK_STANDALONE)
	ADD_CUSTOM_TARGET(distclean make clean
		COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake snowpack/cmake_install.cmake  tests/*.cmake
		COMMAND cmake -E remove_directory CMakeFiles
		COMMAND cmake -E remove_directory snowpack/CMakeFiles
		COMMAND cmake -E remove_directory tests/CMakeFiles
		COMMAND cmake -E remove tests/res1exp/*.cmake tests/res1exp/output/*
		COMMAND cmake -E remove_directory tests/res1exp/CMakeFiles
		COMMAND cmake -E remove_directory tests/res1exp/Testing
		)
ENDIF(SNOWPACK_STANDALONE)

###########################################################
#for the uninstall target
CONFIGURE_FILE(
	"${PROJECT_SOURCE_DIR}/tools/cmake/cmake_uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

#Build the documentation if doxygen is found
INCLUDE(FindDoxygen)
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND)
	ADD_CUSTOM_TARGET(documentation "doxygen" ${PROJECT_SOURCE_DIR}/snowpack/config.dox)
	GET_TARGET_PROPERTY(DOC_TARGET doc TYPE)
	IF(NOT DOC_TARGET)
		ADD_CUSTOM_TARGET(doc)
	ENDIF(NOT DOC_TARGET)
	ADD_DEPENDENCIES(doc documentation)
ENDIF (DOXYGEN_FOUND)
#Install the documentation if the user desires it
SET(INSTALL_DOC OFF CACHE BOOL "Install documentation ON or OFF")
IF(INSTALL_DOC)
	IF(WIN32)
		SET(DOC_PATH "doc")
	ELSE(WIN32)
		SET(DOC_PATH "share/doc/${CMAKE_PROJECT_NAME}")
	ENDIF(WIN32)
	SET(DOC_FILES License.txt Readme.txt snowpack/icons/icon.ico)
	IF (DOXYGEN_FOUND)
		INSTALL(DIRECTORY doc/html DESTINATION ${DOC_PATH} COMPONENT documentation)
		LIST(APPEND DOC_FILES doc/Snowpack_doc.html)
	ELSE (DOXYGEN_FOUND)
		MESSAGE("Doxygen has not been found and is needed to build online documentation. Please either install doxygen or download a documentation-only package!")
	ENDIF (DOXYGEN_FOUND)

	INSTALL(FILES ${DOC_FILES} DESTINATION ${DOC_PATH} COMPONENT documentation)
	INSTALL(DIRECTORY doc/examples USE_SOURCE_PERMISSIONS DESTINATION ${DOC_PATH} COMPONENT examples PATTERN ".svn" EXCLUDE)
ENDIF(INSTALL_DOC)

###########################################################
## CPACK section
###########################################################
SET(CPACK_BUNDLE_ALL OFF CACHE BOOL "Pack all necessary libs and tools for a standalone release")
IF(CPACK_BUNDLE_ALL)
	#bundle MeteoIO
	FIND_PACKAGE(MeteoIO REQUIRED)
	FILE(GLOB meteoio_libs "${METEOIO_ROOT}/lib/libmeteo*.*") #this contains dynamic, static and import libs
	INSTALL(DIRECTORY ${METEOIO_ROOT}/include/meteoio DESTINATION include COMPONENT headers)
	INSTALL(FILES ${meteoio_libs} DESTINATION lib COMPONENT libraries)
	INSTALL(DIRECTORY ${METEOIO_ROOT}/lib/meteoio/plugins DESTINATION lib/meteoio COMPONENT libraries)
	IF(INSTALL_DOC)
		IF(WIN32)
			SET(METEOIO_INSTALLED_DOC "${METEOIO_ROOT}/doc")
			SET(METEOIO_DOC ${DOC_PATH}/meteoio)
		ELSE(WIN32)
			SET(METEOIO_INSTALLED_DOC "${METEOIO_ROOT}/share/doc/meteoio")
			SET(METEOIO_DOC ${DOC_PATH}/../meteoio)
		ENDIF(WIN32)
		INSTALL(DIRECTORY "${METEOIO_INSTALLED_DOC}/examples" DESTINATION ${METEOIO_DOC} COMPONENT examples)
		INSTALL(DIRECTORY "${METEOIO_INSTALLED_DOC}/html" DESTINATION ${METEOIO_DOC} COMPONENT documentation)
		FILE(GLOB meteoio_docs "${METEOIO_INSTALLED_DOC}/*.*")
		INSTALL(FILES ${meteoio_docs} DESTINATION ${METEOIO_DOC} COMPONENT documentation)
	ENDIF(INSTALL_DOC)

	#bundle inishell
	FILE(GLOB inishell "${CMAKE_SOURCE_DIR}/tools/inish*.jar")
	INSTALL(FILES "${inishell}" DESTINATION bin COMPONENT exe)
ENDIF(CPACK_BUNDLE_ALL)

SET(CPACK_PACKAGE_DESCRIPTION "Snowpack Library")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "This is the SNOWPACK model ([Bartelt 2002], [Lehning 2002a], [Lehning 2002b]) for physical modeling of snow pack as a library. It can therefore be integrated into other models.")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_PACKAGE_NAME "Snowpack")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Snowpack")
SET(CPACK_PACKAGE_CONTACT "Mathias Bavay ")
SET(CPACK_PACKAGE_VENDOR "WSL Institute for Snow and Avalanche Research")
SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/snowpack/icons\\\\icon.ico")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/License.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/Readme.txt")

SET(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")
SET(VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
#SET(CPACK_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
SET(CPACK_ARCH "i386")
SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}${PACKAGE_EXTRA}-${VERSION}-${CPACK_ARCH}")

#installation by components for formats that support it
#SET(CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE 1)
SET(CPACK_COMPONENTS_ALL libraries headers documentation examples exe)
SET(CPACK_COMPONENT_EXE_DISPLAY_NAME "Snowpack application")
SET(CPACK_COMPONENT_EXE_DESCRIPTION "The SNOWPACK standalone model")
SET(CPACK_COMPONENT_EXE_GROUP Runtime)
SET(CPACK_COMPONENT_EXE_DEPENDS libraries)
SET(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "Snowpack library")
SET(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "The SNOWPACK library itself")
SET(CPACK_COMPONENT_LIBRARIES_GROUP Runtime)
SET(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "Header files")
SET(CPACK_COMPONENT_HEADERS_DESCRIPTION "Necessary C++ header files for developing applications using SNOWPACK")
SET(CPACK_COMPONENT_HEADERS_GROUP Development)
SET(CPACK_COMPONENT_HEADERS_DEPENDS libraries)
IF(INSTALL_DOC)
	SET(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "Documentation")
	SET(CPACK_COMPONENT_DOCUMENTATION_DESCRIPTION "Documentation of the SNOWPACK usage and API both for the end user and the developer using SNOWPACK in his/her application")
	SET(CPACK_COMPONENT_DOCUMENTATION_GROUP Docs)
	SET(CPACK_COMPONENT_EXAMPLES_DISPLAY_NAME "Examples")
	SET(CPACK_COMPONENT_EXAMPLES_DESCRIPTION "Various simulations with their example data")
	SET(CPACK_COMPONENT_EXAMPLES_GROUP Docs)
	#SET(CPACK_COMPONENT_EXAMPLES_DEPENDS libraries headers runtime)
	SET(CPACK_COMPONENT_DOCS_GROUP "Documentation")
	SET(CPACK_COMPONENT_GROUP_DOCS_DESCRIPTION "How to use Snowpack as a standalone model or as a library called by your own model")
ENDIF(INSTALL_DOC)

SET(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")
SET(CPACK_COMPONENT_GROUP_RUNTIME_DESCRIPTION "What you need to run an application based on SNOWPACK")
SET(CPACK_COMPONENT_DEVELOPMENT_GROUP "Development")
SET(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "What you need to develop an application based on SNOWPACK")


IF(UNIX)
	IF(APPLE)
		SET(CPACK_GENERATOR "PackageMaker;")
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/Applications/${CPACK_PACKAGE_NAME}-${VERSION}")
		SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/snowpack/icons/icon.icns")
		SET(CPACK_BUNDLE_ICON "${CMAKE_SOURCE_DIR}/snowpack/icons/icon.icns")
		SET_SOURCE_FILES_PROPERTIES(${CPACK_BUNDLE_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
		INSTALL(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink bin/snowpack /usr/bin/snowpack)")
	ELSE(APPLE)
		SET(CPACK_GENERATOR "DEB;RPM;TGZ")
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
		SET(CPACK_DEBIAN_PACKAGE_NAME "snowpack")
		SET(CPACK_DEBIAN_PACKAGE_SECTION "science")
		SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6, meteoio") #objdump -p lib/libmeteoio.so | grep NEEDED | tr -s ' ' | cut -d' ' -f3 | tr '\n' ','
		SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386") #dpkg --print-architecture
		SET(CPACK_RPM_PACKAGE_NAME "snowpack")
		SET(CPACK_RPM_PACKAGE_LICENSE "LPGLv3")
		SET(CPACK_RPM_PACKAGE_REQUIRES "libstdc++6, meteoio")
		SET(CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
		SET(CPACK_RPM_PACKAGE_GROUP "Applications/Science")
	ENDIF(APPLE)
ELSE(UNIX)
	IF(WIN32)
		SET(CPACK_GENERATOR "NSIS")
		SET(CPACK_NSIS_MUI_ICON "${CPACK_PACKAGE_ICON}")
		SET(CPACK_PACKAGE_ICON "${CPACK_PACKAGE_ICON}")
		SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}-${VERSION}")
		SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME}")
		#SET(CPACK_NSIS_EXTRA_INSTALL_COMMANDS "${EnvVarUpdate} $0 \"LIB\"  \"A\" \"HKLM\" \"C:\MyLib\"")
		SET(CPACK_NSIS_MODIFY_PATH ON)
		IF(INSTALL_DOC)
			IF(CPACK_BUNDLE_ALL)
				SET(CPACK_NSIS_MENU_LINKS
				"bin/inishell.jar" "inishell"
				"${DOC_PATH}/Snowpack_doc.html" "SNOWPACK documentation"
				"${DOC_PATH}/examples" "SNOWPACK examples"
				"https://slfsmm.indefero.net/p/snowpack/" "Home page"
				"${METEOIO_DOC}/MeteoIO_doc.html" "MeteoIO documentation"
				"${METEOIO_DOC}/examples" "MeteoIO examples"
				"https://slfsmm.indefero.net/p/meteoio/" "MeteoIO Home page")
			ELSE(CPACK_BUNDLE_ALL)
				SET(CPACK_NSIS_MENU_LINKS
				"${DOC_PATH}/Snowpack_doc.html" "SNOWPACK documentation"
				"${DOC_PATH}/examples" "SNOWPACK examples"
				"https://slfsmm.indefero.net/p/snowpack/" "Home page")
			ENDIF(CPACK_BUNDLE_ALL)
		ELSEIF(INSTALL_DOC)
			IF(CPACK_BUNDLE_ALL)
				SET(CPACK_NSIS_MENU_LINKS
				"bin/inishell.jar" "inishell"
				"https://slfsmm.indefero.net/p/snowpack/" "Home page"
				"https://slfsmm.indefero.net/p/meteoio/" "MeteoIO Home page")
			ELSE(CPACK_BUNDLE_ALL)
				SET(CPACK_NSIS_MENU_LINKS "https://slfsmm.indefero.net/p/snowpack/" "Home page")
			ENDIF(CPACK_BUNDLE_ALL)
		ENDIF(INSTALL_DOC)

		#info for the add/uninstall software administration
		SET(CPACK_NSIS_INSTALLED_ICON_NAME "doc\\\\icon.ico")
		SET(CPACK_NSIS_URL_INFO_ABOUT "http://www.slf.ch")
		SET(CPACK_NSIS_HELP_LINK "https://slfsmm.indefero.net/p/snowpack")
		SET(CPACK_NSIS_CONTACT "Mathias.Bavay@wsl.ch")
	ENDIF(WIN32)
ENDIF(UNIX)


SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${VERSION}" CACHE INTERNAL "tarball basename")
SET(CPACK_SOURCE_GENERATOR TGZ)
# The following components are regex's to match anywhere (unless anchored)
# in absolute path + filename to find files or directories to be excluded
# from source tarball.
SET(CPACK_SOURCE_IGNORE_FILES
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

INCLUDE(CPack)

###########################################################
## CTest section
###########################################################
OPTION(BUILD_TESTING "Build the QA tests" OFF)
IF(BUILD_TESTING)
    INCLUDE(CTest)
    ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTING)
