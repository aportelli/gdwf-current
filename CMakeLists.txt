cmake_minimum_required(VERSION 3.11.0)

# set the app name here
set(APP_NAME grid-app)

project(
  ${APP_NAME}
  VERSION 1.0
  LANGUAGES C CXX)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(GNUInstallDirs)
include(FetchContent)

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

set(FETCHCONTENT_BASE_DIR ${CMAKE_BINARY_DIR}/deps)
FetchContent_Declare(
  FindGrid
  GIT_REPOSITORY https://github.com/aportelli/FindGrid.git
  GIT_TAG main)
FetchContent_MakeAvailable(FindGrid)
list(APPEND CMAKE_MODULE_PATH ${findgrid_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

add_subdirectory(src)
