# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

# Build everything as static libs.
set (BUILD_SHARED_LIBS OFF)

# Project and version numbers.
project (polywog)
set (POLYWOG_MAJOR_VERSION 0)
set (POLYWOG_MINOR_VERSION 6)
set (POLYWOG_PATCH_VERSION 0)
set (POLYWOG_VERSION "${POLYWOG_MAJOR_VERSION}.${POLYWOG_MINOR_VERSION}.${POLYWOG_PATCH_VERSION}")
message("-- Building polywog (v${POLYWOG_VERSION})")

# Use our own CMake stuff.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Look for Polymec and use its CMake settings.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${POLYMEC_PREFIX}/share/polymec/")
include(polymec)
if (NOT DEFINED POLYMEC_VERSION)
  message(FATAL_ERROR "Could not find polymec library. Please use the polymec=/path/to/polymec option.")
endif()
message("-- Found polymec library in ${POLYMEC_PREFIX} (v${POLYMEC_VERSION})")

# Pass POLYWOG_VERSION and friends as definitions to the compiler.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOLYWOG_VERSION=\\\"${POLYWOG_VERSION}\\\"")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOLYWOG_MAJOR_VERSION=${POLYWOG_MAJOR_VERSION}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOLYWOG_MINOR_VERSION=${POLYWOG_MINOR_VERSION}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOLYWOG_PATCH_VERSION=${POLYWOG_PATCH_VERSION}")

# Look for headers in the top directory.
include_directories(${PROJECT_SOURCE_DIR})

# Look for libraries in these directories.
link_directories("${POLYMEC_PREFIX}/lib")

enable_testing()

add_subdirectory(polywog)

# Generate a polywog.cmake file that contains installation information.
set(CMAKE_INSTALL_PREFIX ${POLYMEC_PREFIX})
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/polywog.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polywog.cmake"
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/polywog.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polywog_executable.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polywog_test.cmake DESTINATION share/polymec)

