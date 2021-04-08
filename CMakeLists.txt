cmake_minimum_required(VERSION 3.13)


# ==============================================================================
# Project info

project(ROADEF_Challenge_2020
        VERSION "0.0.1"
        DESCRIPTION "ROADEF/EURO Challenge 2020: Maintenance Planning Problem"
        HOMEPAGE_URL "https://github.com/andremaravilha/ROADEF-Challenge-2020")


# ==============================================================================
# C++ standard and other settings

set(CXX ON)
enable_language(CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

# If build type is not set, then use "Release"
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Set MS Visual Studio compiler with static runtime libraries, if MSVS is the compiler
if(MSVC AND MT)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()


# ==============================================================================
# Third-party dependencies

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(GUROBI REQUIRED)


# ==============================================================================
# Paths to search for headers and libraries

include_directories(
        ${PROJECT_SOURCE_DIR}/lib/
        ${GUROBI_INCLUDE_DIRS})


# ==============================================================================
# Source files

set(SOURCE_FILES
        src/main.cpp
        src/roadef2020.h src/roadef2020.cpp
        src/heur_greedy.h src/heur_greedy.cpp
        src/heur_benders.h src/heur_benders.cpp)


# ==============================================================================
# Targets

add_executable(challengeRTE ${SOURCE_FILES})
target_link_libraries(challengeRTE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
target_link_libraries(challengeRTE ${GUROBI_LIBRARY})
