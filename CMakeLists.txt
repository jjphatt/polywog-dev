# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

# Project and version numbers.
project (polywog)
set (POLYWOG_MAJOR_VERSION 0)
set (POLYWOG_MINOR_VERSION 7)
set (POLYWOG_PATCH_VERSION 0)
set (POLYWOG_VERSION "${POLYWOG_MAJOR_VERSION}.${POLYWOG_MINOR_VERSION}.${POLYWOG_PATCH_VERSION}")
message("-- Building polywog (v${POLYWOG_VERSION})")

# Use our own CMake stuff.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Look for headers first the top directory first.
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

# Look for Polymec and use its CMake settings.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${POLYMEC_PREFIX}/share/polymec/")
include(polymec)
if (NOT DEFINED POLYMEC_VERSION)
  message(FATAL_ERROR "Could not find polymec library. Please use the polymec=/path/to/polymec option.")
endif()
message("-- Found polymec library in ${POLYMEC_PREFIX} (v${POLYMEC_VERSION})")

# Build everything the way Polymec does.
set(CMAKE_C_COMPILER ${POLYMEC_C_COMPILER})
set(CMAKE_C_FLAGS ${POLYMEC_C_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS ${POLYMEC_EXE_LINKER_FLAGS})
set(BUILD_SHARED_LIBS ${POLYMEC_BUILD_SHARED_LIBS})
set(CMAKE_SHARED_LINKER_FLAGS ${POLYMEC_SHARED_LINKER_FLAGS})
set(CMAKE_BUILD_TYPE ${POLYMEC_BUILD_TYPE})
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${POLYMEC_PREFIX}/lib")

# Generate polywog_version.h with the right version numbers.
add_custom_target(update_version_h ALL
                  python "${POLYMEC_PREFIX}/share/polymec/update_version_h.py" polywog ${POLYWOG_VERSION} ${PROJECT_BINARY_DIR}/polywog/polywog_version.h)

# Generate polyamri.h at the top level.
add_custom_target(generate_top_header ALL
                  python ${PROJECT_SOURCE_DIR}/tools/generate_top_header.py ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)
                  
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
install(FILES ${PROJECT_BINARY_DIR}/include/polywog.h DESTINATION include/polymec)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/polywog.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polywog_executable.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polywog_test.cmake DESTINATION share/polymec)

