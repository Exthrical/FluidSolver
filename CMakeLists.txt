cmake_minimum_required(VERSION 3.15)
project(FluidSolver2D CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options for configuration
option(BUILD_VERBOSE "Enable verbose build output" OFF)

if(BUILD_VERBOSE)
  set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

# ImGui directory (headers + sources). You can point this to a local Dear ImGui checkout
# containing imgui.h and the backends/ folder. If not provided, auto-detect third_party/imgui.
set(IMGUI_DIR "" CACHE PATH "Path to Dear ImGui root containing imgui.h and backends/")

# Auto-detect bundled ImGui if not provided
if(NOT IMGUI_DIR)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui/imgui.h")
    set(IMGUI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui")
    message(STATUS "Auto-detected IMGUI_DIR: ${IMGUI_DIR}")
  endif()
endif()

find_package(OpenGL REQUIRED)

# Try to find GLFW via config or module
find_package(glfw3 QUIET CONFIG)
if(NOT glfw3_FOUND)
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(GLFW QUIET glfw3)
  endif()
endif()

add_executable(FluidSolver
  src/main.cpp
  src/flip_solver.cpp
  src/renderer.cpp
)

target_include_directories(FluidSolver PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# If IMGUI_DIR is set, add includes and sources for core + backends used
if(IMGUI_DIR)
  message(STATUS "Using IMGUI_DIR: ${IMGUI_DIR}")
  target_include_directories(FluidSolver PRIVATE ${IMGUI_DIR} ${IMGUI_DIR}/backends)
  # ImGui core sources + GLFW/OpenGL2 backends
  target_sources(FluidSolver PRIVATE
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui_demo.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl2.cpp
  )
else()
  message(WARNING "IMGUI_DIR not set and bundled third_party/imgui not found. ImGui headers will not be available.")
endif()

# Link OpenGL
target_link_libraries(FluidSolver PRIVATE OpenGL::GL)

# Link GLFW
if(glfw3_FOUND)
  target_link_libraries(FluidSolver PRIVATE glfw)
elseif(GLFW_FOUND)
  target_include_directories(FluidSolver PRIVATE ${GLFW_INCLUDE_DIRS})
  target_link_libraries(FluidSolver PRIVATE ${GLFW_LIBRARIES})
else()
  message(WARNING "GLFW not found via CMake or pkg-config. Please install GLFW3 or provide a toolchain that exposes glfw3 target.")
endif()

# Windows specifics
if(WIN32)
  target_compile_definitions(FluidSolver PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

# macOS: link Cocoa frameworks when needed (GLFW may pull them in)
if(APPLE)
  find_library(COCOA_LIBRARY Cocoa)
  find_library(IOKIT_LIBRARY IOKit)
  find_library(COREVIDEO_LIBRARY CoreVideo)
  if(COCOA_LIBRARY AND IOKIT_LIBRARY AND COREVIDEO_LIBRARY)
    target_link_libraries(FluidSolver PRIVATE ${COCOA_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY})
  endif()
endif()

# Enable more warnings
if(MSVC)
  target_compile_options(FluidSolver PRIVATE /W4 /permissive-)
else()
  target_compile_options(FluidSolver PRIVATE -Wall -Wextra -Wpedantic)
endif()
