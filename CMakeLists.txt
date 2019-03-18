cmake_minimum_required(VERSION 3.6)

project(UNITIGS)
set (CMAKE_CXX_STANDARD 11)
#
# <usual CMakeLists.txt stuff goes here...>
#
################################################################################
# The version number.
################################################################################
SET (gatb-tool_VERSION_MAJOR 1)
SET (gatb-tool_VERSION_MINOR 0)
SET (gatb-tool_VERSION_PATCH 2)

#write the version number
file(WRITE ${PROJECT_SOURCE_DIR}/src/version.h "\#define VERSION \"${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH}\"")
file(WRITE ${PROJECT_SOURCE_DIR}/unitig-graph/version.h "\#define VERSION \"${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH}\"")


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
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/gatb-core/gatb-core)
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

# this stops the boost in gatb-core/thirdparty being used
if(EXISTS ${GATB_CORE_HOME}/thirdparty/boost)
    file(RENAME ${GATB_CORE_HOME}/thirdparty/boost ${GATB_CORE_HOME}/thirdparty/boost-old)
endif()
execute_process(COMMAND git apply ${PROJECT_SOURCE_DIR}/0001-Fix-for-boost-include-dir.patch
                WORKING_DIRECTORY ${GATB_CORE_HOME})

# GATB CORE
include (GatbCore)

################################################################################
#  TOOLS
################################################################################

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

# We add the pstream include directory
set (PSTREAMS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/pstreams/)
include_directories (${PSTREAMS_SOURCE_DIR})

set(BOOST_ROOT /usr/bin/include)
set(Boost_USE_STATIC_LIBS   ON)
FIND_PACKAGE(Boost 1.69.0 COMPONENTS filesystem system regex REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# We add the gatb-core include directory
include_directories (${gatb-core-includes})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})
link_directories(${Boost_LIBRARY_DIR})

# cdbg-ops target
set(PROGRAM cdbg-ops)
set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/unitig-graph)
include_directories (${PROGRAM_SOURCE_DIR})
add_library(${PROGRAM}obj OBJECT ${PROGRAM_SOURCE_DIR}/node_dists.cpp)
add_library(${PROGRAM} STATIC $<TARGET_OBJECTS:${PROGRAM}obj>)
add_executable(${PROGRAM} $<TARGET_OBJECTS:${PROGRAM}obj> ${PROGRAM_SOURCE_DIR}/graph_ops.cpp)

# unitig-counter target
set(PROGRAM unitig-counter)
set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
include_directories (${PROGRAM_SOURCE_DIR})
file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/*.cpp)
add_executable(${PROGRAM} ${ProjectFiles})
target_link_libraries(${PROGRAM} ${gatb-core-libraries} ${Boost_LIBRARIES} -static-libgcc -static-libstdc++)

################################################################################
#  INSTALLATION
################################################################################
install (TARGETS ${PROGRAMS} DESTINATION bin COMPONENT precompiled)

################################################################################
# Packaging
################################################################################
# We set the version number
SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "gatb-tool ${PROJECT_NAME}")
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
