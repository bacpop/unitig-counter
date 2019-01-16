cmake_minimum_required(VERSION 3.6)

# Handle superbuild first
option (USE_SUPERBUILD "Whether or not a superbuild should be invoked" ON)

if (USE_SUPERBUILD)
  project (SUPERBUILD NONE)
  # execute the superbuild (this script will be invoked again without the
  # USE_SUPERBUILD option this time)
  include (boostSuperBuild/SuperBuild.cmake)
  return()                      # stop processing this file further
else()
  project(DBGWAS)
endif()

#
# <usual CMakeLists.txt stuff goes here...>
#
################################################################################
# The version number.
################################################################################
SET (gatb-tool_VERSION_MAJOR 0)
SET (gatb-tool_VERSION_MINOR 5)
SET (gatb-tool_VERSION_PATCH 2)

#write the version number
file(WRITE ${PROJECT_SOURCE_DIR}/tools/DBGWAS/src/version.h "\#define VERSION \"${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH}\"")


IF (DEFINED MAJOR)
    SET (gatb-tool_VERSION_MAJOR ${MAJOR})
ENDIF()
IF (DEFINED MINOR)
    SET (gatb-tool_VERSION_MINOR ${MINOR})
ENDIF()
IF (DEFINED PATCH)
    SET (gatb-tool_VERSION_PATCH ${PATCH})
ENDIF()

set (gatb-tool-version ${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH})

################################################################################
# Define cmake modules directory
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/gatb/gatb-core)
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

################################################################################
# SUPPORTED KMER SIZES
################################################################################

# One can uncomment this line and set the wanted values
set (KSIZE_LIST "32   64   96  128")

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1)
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
#  TOOLS
################################################################################
#ADD_SUBDIRECTORY(tools)
MACRO(SUBDIRLIST result curdir)
    FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
    SET (dirlist "")
    FOREACH(child ${children})
      IF(IS_DIRECTORY ${curdir}/${child})
          LIST(APPEND dirlist ${child})
      ENDIF()
    ENDFOREACH()
    SET(${result} ${dirlist})
ENDMACRO()

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

# We add the gatb-core include directory
include_directories (${gatb-core-includes})

# We add the pstream include directory
set (PSTREAMS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/pstreams/)
include_directories (${PSTREAMS_SOURCE_DIR})

# We add the whereami include directory
set (WHEREAMI_SOURCE_DIR ${PROJECT_SOURCE_DIR}/whereami/src)
include_directories (${WHEREAMI_SOURCE_DIR})

set(BOOST_ROOT ${CMAKE_BINARY_DIR}/../Dependencies/Source/boost_1_62_0)
set(Boost_USE_STATIC_LIBS   ON) #This makes boost plug the libraries we specify as static libraries (not dynamic ones - the user does not even need to have boost installed in their pc)
FIND_PACKAGE(Boost 1.62.0 COMPONENTS filesystem system regex REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})
link_directories(${Boost_LIBRARY_DIR})

# We get the list of subdirectories names
SUBDIRLIST(PROGRAMS ${PROJECT_SOURCE_DIR}/tools)

# We create a target for each
FOREACH (program ${PROGRAMS})
  set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/tools/${program}/src)
  include_directories (${PROGRAM_SOURCE_DIR})
  file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/*.cpp)
  file (GLOB_RECURSE  WhereamiFiles  ${WHEREAMI_SOURCE_DIR}/*.c)
  add_executable(${program} ${ProjectFiles} ${WhereamiFiles})
  target_link_libraries(${program} ${gatb-core-libraries} ${Boost_LIBRARIES} -static-libgcc -static-libstdc++)
ENDFOREACH (program)

################################################################################
#  INSTALLATION
################################################################################
install (TARGETS ${PROGRAMS} DESTINATION bin COMPONENT precompiled)

################################################################################
# Packaging
################################################################################
# We set the version number
SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "gatb-tool ${PROJECT_NAME}")
SET (CPACK_PACKAGE_VENDOR               "bioMerieux; LBBE; Erable")
SET (CPACK_PACKAGE_VERSION_MAJOR        "${gatb-tool_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${gatb-tool_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${gatb-tool_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${gatb-tool-version}")  

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# Packaging the source ; we ignore unwanted files 
SET (CPACK_SOURCE_IGNORE_FILES
        "^${PROJECT_SOURCE_DIR}/build/"
        "^${PROJECT_SOURCE_DIR}/builddebug/"
        "^${PROJECT_SOURCE_DIR}/buildlinux/"
        "^${PROJECT_SOURCE_DIR}/TODO"
        "^${GATB_CORE_HOME}/.project"
        "^${GATB_CORE_HOME}/.gitignore"
        "^${GATB_CORE_HOME}/doc"
        "^${GATB_CORE_HOME}/DELIVERY.md"
)

# Packaging the binary ; we want to include some additional files
INSTALL (FILES   ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE                     DESTINATION . COMPONENT precompiled)
INSTALL (FILES   ${CMAKE_CURRENT_SOURCE_DIR}/THIRDPARTIES                     DESTINATION . COMPONENT precompiled)
INSTALL (FILES   ${CMAKE_CURRENT_SOURCE_DIR}/README.md                     DESTINATION . COMPONENT precompiled)

set(CPACK_ARCHIVE_COMPONENT_INSTALL ON)
set(CPACK_COMPONENTS_ALL precompiled)
include (CPack)
