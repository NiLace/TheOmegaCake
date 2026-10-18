# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) The NL Sounds contributors
#
# The Omega Cake -- build the LV2 bundle (Linux, x86-64).
#
#   cmake -S . -B build && cmake --build build -j
#
# The bundle is build/TheOmegaCake.lv2/; copy that directory into ~/.lv2 (or another LV2_PATH directory) to install
# it. The baseline is generic x86-64: the bundle runs on any 64-bit Intel or AMD CPU, and the one hot loop carries a
# second copy for CPUs with AVX2 and FMA, chosen at load time.
cmake_minimum_required(VERSION 3.21)

project(TheOmegaCake
    VERSION     1.0.2
    DESCRIPTION "The Omega Cake - fuzz pedal modelling plugin (LV2), 45 models"
    LANGUAGES   CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)            # the UI decodes and scales its images on a few threads
pkg_check_modules(CAIRO      REQUIRED IMPORTED_TARGET cairo)
pkg_check_modules(CAIRO_XLIB REQUIRED IMPORTED_TARGET cairo-xlib)
pkg_check_modules(X11        REQUIRED IMPORTED_TARGET x11)
pkg_check_modules(FREETYPE   REQUIRED IMPORTED_TARGET freetype2)

set(TOCAKE_LV2_BUNDLE "${CMAKE_BINARY_DIR}/TheOmegaCake.lv2")

# The DSP: header-only engines behind one LV2 entry point. Only lv2_descriptor is exported.
add_library(theomegacake_lv2 MODULE src/lv2.cpp)
set_target_properties(theomegacake_lv2 PROPERTIES
    PREFIX ""
    OUTPUT_NAME "theomegacake"
    LIBRARY_OUTPUT_DIRECTORY "${TOCAKE_LV2_BUNDLE}"
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON)
target_include_directories(theomegacake_lv2 PRIVATE src)
target_compile_options(theomegacake_lv2 PRIVATE -Wall -Wextra -Wpedantic)

# The UI: an Xlib window painted with Cairo; the fonts are loaded from the bundle through FreeType.
add_library(theomegacake_ui MODULE src/ui_x11.cpp)
set_target_properties(theomegacake_ui PROPERTIES
    PREFIX ""
    OUTPUT_NAME "theomegacake_ui"
    LIBRARY_OUTPUT_DIRECTORY "${TOCAKE_LV2_BUNDLE}")
target_include_directories(theomegacake_ui PRIVATE src)
target_link_libraries(theomegacake_ui PRIVATE
    PkgConfig::CAIRO PkgConfig::CAIRO_XLIB PkgConfig::X11 PkgConfig::FREETYPE Threads::Threads)
target_compile_options(theomegacake_ui PRIVATE -Wall -Wextra -Wpedantic)
# Without it the headers' inline variables are STB_GNU_UNIQUE symbols, and the loader keeps the module mapped
# after a host's dlclose for the life of the process.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(theomegacake_ui PRIVATE -fno-gnu-unique)
endif()

# The bundle: the two modules, the TTL files, the fonts with their licence texts, and the panel images.
foreach(_ttl manifest.ttl theomegacake.ttl theomegacake_ui.ttl)
    configure_file(${CMAKE_SOURCE_DIR}/lv2/${_ttl} ${TOCAKE_LV2_BUNDLE}/${_ttl} COPYONLY)
endforeach()
foreach(_font Doto-Round-Bold.ttf SpaceMono-400.ttf SpaceMono-700.ttf OFL-Doto.txt OFL-SpaceMono.txt)
    configure_file(${CMAKE_SOURCE_DIR}/assets/fonts/${_font} ${TOCAKE_LV2_BUNDLE}/${_font} COPYONLY)
endforeach()
file(GLOB TOCAKE_PHOTO_PNGS ${CMAKE_SOURCE_DIR}/gui/photo/*.png)
foreach(_png ${TOCAKE_PHOTO_PNGS})
    get_filename_component(_name ${_png} NAME)
    configure_file(${_png} ${TOCAKE_LV2_BUNDLE}/photo/${_name} COPYONLY)
endforeach()
