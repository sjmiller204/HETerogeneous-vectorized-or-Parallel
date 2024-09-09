# HETP/CMakeLists.txt

cmake_minimum_required(VERSION 3.5)
project(HETP
  VERSION 1.0
  LANGUAGES Fortran
)

#-----------------------------------------------------------------------------
# Set CMake policies.  For more information, see:
#
# https://cmake.org/cmake/help/latest/policy/CMP0054.html
# https://cmake.org/cmake/help/latest/policy/CMP0057.html
# https://cmake.org/cmake/help/latest/policy/CMP0074.html
# https://cmake.org/cmake/help/latest/policy/CMP0079.html
#-----------------------------------------------------------------------------
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0057 NEW)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
if(POLICY CMP0079)
  cmake_policy(SET CMP0079 NEW)
endif()

#-----------------------------------------------------------------------------
# Add CMakeScripts/ to the module path
#-----------------------------------------------------------------------------
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/CMakeScripts)
include(HETP-Helpers)

#-----------------------------------------------------------------------------
# Print header with the CMake project version and the HETP repo version
#-----------------------------------------------------------------------------
message("=================================================================")
message("HETP ${PROJECT_VERSION}")
message("=================================================================")

#-----------------------------------------------------------------------------
# Declare the HETPBuildProperties
#-----------------------------------------------------------------------------
add_library(HETPBuildProperties INTERFACE)

#---------------------------------------------------------------------
# Default compiler options
#
# HETP_Fortran_FLAGS_<COMPILER_ID> : common flags (all build types) 
# HETP_Fortran_FLAGS_<CONFIG>_<COMPILER_ID> : types-specific flags
#
# Valid COMPILER_ID: Intel and GNU
# Valid CONFIG: RELEASE, DEBUG, RELWITHDEBINFO
#---------------------------------------------------------------------
set(HETP_DETECTED_FORTRAN_COMPILER_ID ${CMAKE_Fortran_COMPILER_ID}
    CACHE INTERNAL "Logging the COMPILER_ID to CMakeCache.txt"
)
set(HETP_DETECTED_FORTRAN_COMPILER_VERSION ${CMAKE_Fortran_COMPILER_VERSION}
    CACHE INTERNAL "Logging the compiler version to CMakeCache.txt"
)

set(HETP_Fortran_FLAGS_Intel
   -cpp -w -auto -noalign "SHELL:-convert big_endian" "SHELL:-fp-model source" -mcmodel=medium
   -shared-intel -traceback -DLINUX_IFORT
   CACHE STRING "HETP compiler flags for all build types with Intel compilers"
)
set(HETP_Fortran_FLAGS_RELEASE_Intel
   -O2
   CACHE STRING "HETP compiler flags for build type release with Intel compilers"
)
set(HETP_Fortran_FLAGS_RELWITHDEBINFO_Intel
   -O2
   CACHE STRING "HETP compiler flags for build type relwithdebinfo with Intel compilers"
)
set(HETP_Fortran_FLAGS_DEBUG_Intel
   -g -O0 "SHELL:-check arg_temp_created" "SHELL:-debug all" -fpe0 -ftrapuv -check,bounds -DDEBUG
   CACHE STRING "HETP compiler flags for build type debug with Intel compilers"
)

set(HETP_Fortran_FLAGS_GNU
   -cpp -w -std=legacy -fautomatic -fno-align-commons -fconvert=big-endian
   -fno-range-check -mcmodel=medium -fbacktrace -g -DLINUX_GFORTRAN
   -ffree-line-length-none
   CACHE STRING "HETP compiler flags for all build types with GNU compilers"
)
set(HETP_Fortran_FLAGS_RELEASE_GNU
   -O3 -funroll-loops
   CACHE STRING "HETP compiler flags for build type release with GNU compilers"
)
set(HETP_Fortran_FLAGS_RELWITHDEBINFO_GNU
   -O3 -funroll-loops
   CACHE STRING "HETP compiler flags for build type relwithdebinfo with GNU compilers"
)
set(HETP_Fortran_FLAGS_DEBUG_GNU
   -g -gdwarf-2 -gstrict-dwarf -O0 -Wall -Wextra -Wconversion -Warray-temporaries 
   -fcheck=array-temps -ffpe-trap=invalid,zero,overflow -finit-real=snan -fcheck=bounds -fcheck=pointer
   CACHE STRING "HETP compiler flags for build type debug with GNU compilers"
)

set(HETP_SUPPORTED_COMPILER_IDS "Intel" "GNU")
if(NOT CMAKE_Fortran_COMPILER_ID IN_LIST HETP_SUPPORTED_COMPILER_IDS)
   message(FATAL_ERROR "HETP does not support ${CMAKE_Fortran_COMPILER_ID} compilers")
endif()

#---------------------------------------------------------------------
# Assign comiler options to build properties
#---------------------------------------------------------------------
target_compile_options(HETPBuildProperties
   INTERFACE
   $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},Intel>:
      ${HETP_Fortran_FLAGS_Intel}
      $<$<CONFIG:Debug>:${HETP_Fortran_FLAGS_DEBUG_Intel}>
      $<$<CONFIG:RelWithDebInfo>:${HETP_Fortran_FLAGS_RELWITHDEBINFO_Intel}>
      $<$<CONFIG:Release>:${HETP_Fortran_FLAGS_RELEASE_Intel}>
   >
   $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},GNU>:
      ${HETP_Fortran_FLAGS_GNU}
      $<$<CONFIG:Debug>:${HETP_Fortran_FLAGS_DEBUG_GNU}>
      $<$<CONFIG:RelWithDebInfo>:${HETP_Fortran_FLAGS_RELWITHDEBINFO_GNU}>
      $<$<CONFIG:Release>:${HETP_Fortran_FLAGS_RELEASE_GNU}>
   >
)

#-----------------------------------------------------------------------------
# Put all of HETP's mod files in build subdir called mod
#-----------------------------------------------------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY ${HETP_BINARY_DIR}/mod)
target_include_directories(HETPBuildProperties
    INTERFACE ${HETP_BINARY_DIR}/mod
)

#-----------------------------------------------------------------------------
# For HETP testing
#
# This conditional block configures the HETP build for HETP standalone test.
# It sets HETP_EXE_TARGETS and it configures the HETPBuildProperties.
#-----------------------------------------------------------------------------
if(NOT HETP_EXTERNAL_CONFIG)

    # Set HETP compile flag to standalone
    target_compile_definitions(HETPBuildProperties
        INTERFACE "HETP_test"
    )

    # Set CMAKE_BUILD_TYPE to Release by default
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release"
            CACHE STRING
            "Set the build type"
            FORCE
    	)
    endif()

    # Run directory
    set(RUNDIR "" CACHE PATH "Path(s) to run directory. Specifies install location.")
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Fake CMAKE_INSTALL_PREFIX" FORCE)
    set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)

    # Determine which executables should be built
    set(HETP_EXE_TARGETS "hetp_test"
        CACHE STRING "Executable targets that get built as a part of \"all\""
    )

endif()

#-----------------------------------------------------------------------------
# Add source code directory
#-----------------------------------------------------------------------------
add_subdirectory(src)

#-----------------------------------------------------------------------------
# Write HETPBuildProperties's configuration to a file
#-----------------------------------------------------------------------------
get_target_property(BT_DEFINITIONS  HETPBuildProperties
    INTERFACE_COMPILE_DEFINITIONS
)
get_target_property(BT_OPTIONS      HETPBuildProperties
    INTERFACE_COMPILE_OPTIONS
)
get_target_property(BT_LIBRARIES    HETPBuildProperties
    INTERFACE_LINK_LIBRARIES
)
get_target_property(BT_INCLUDES     HETPBuildProperties
    INTERFACE_INCLUDE_DIRECTORIES
)
file(WRITE ${CMAKE_BINARY_DIR}/HETPBuildProperties.txt
    "# This file shows the HETPBuildProperties's configuration.\n"
    "\n"
    "HETPBuildProperties::INTERFACE_COMPILE_DEFINITIONS:${BT_DEFINITIONS}\n"
    "HETPBuildProperties::INTERFACE_COMPILE_OPTIONS:${BT_OPTIONS}\n"
    "HETPBuildProperties::INTERFACE_LINK_LIBRARIES:${BT_LIBRARIES}\n"
    "HETPBuildProperties::INTERFACE_INCLUDE_DIRECTORIES:${BT_INCLUDES}\n"
)

#-----------------------------------------------------------------------------
# Copy build information files to each RUNDIR and INSTALLCOPY directory
#----------------------------------------------------------------------------
set(COMBINED_INSTALL_DIRS "")
list(APPEND COMBINED_INSTALL_DIRS ${RUNDIR})
list(APPEND COMBINED_INSTALL_DIRS ${INSTALLCOPY})

# Install to run directories
foreach(INSTALL_PATH ${COMBINED_INSTALL_DIRS})
  if(INSTALL_PATH IN_LIST RUNDIR)
    set(CHECK_IS_RUNDIR TRUE)
  else()
    set(CHECK_IS_RUNDIR FALSE)
  endif()

  # Convert INSTALL_PATH to absolute
  if(NOT IS_ABSOLUTE "${INSTALL_PATH}")
    get_filename_component(INSTALL_PATH "${INSTALL_PATH}" ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}")
  endif()
  install(FILES ${CMAKE_BINARY_DIR}/CMakeCache.txt DESTINATION ${INSTALL_PATH}/build_info)
  #install(PROGRAMS ${CMAKE_SOURCE_DIR}/CMakeScripts/summarize_build DESTINATION ${INSTALL_PATH}/build_info)
endforeach()
