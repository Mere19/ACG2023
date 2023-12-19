# This is the main CMake build file used to compile Nori
cmake_minimum_required (VERSION 3.15..3.27)
project(nori LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(NORI_HEADLESS OFF CACHE BOOL "Compile in headless mode")
set(NORI_COMPILE_LIB OFF CACHE BOOL "Compile lib along the executable")


if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR} )
    message(FATAL_ERROR "In source build (${CMAKE_CURRENT_BINARY_DIR}) are not allowed! Please create a build directory")
endif()

# Please do not change this as it is used for testing
if (DEFINED NORI_CMAKE_INCLUDE)
  include(${NORI_CMAKE_INCLUDE})
else()
  add_subdirectory(ext ext_build)
endif()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  SYSTEM ${EXTERNAL_INCLUDE_DIRS}
)


# The following lines build the main executable. If you add a source
# code file to Nori, be sure to include it in this list.

set(NORI_SOURCE_FILES
  # Header files
  include/nori/bbox.h
  include/nori/bitmap.h
  include/nori/block.h
  include/nori/bsdf.h
  include/nori/bvh.h
  include/nori/camera.h
  include/nori/color.h
  include/nori/common.h
  include/nori/dpdf.h
  include/nori/frame.h
  include/nori/gui.h
  include/nori/integrator.h
  include/nori/emitter.h
  include/nori/kdtree.h
  include/nori/mesh.h
  include/nori/object.h
  include/nori/parser.h
  include/nori/proplist.h
  include/nori/photon.h
  include/nori/ray.h
  include/nori/render.h
  include/nori/rfilter.h
  include/nori/sampler.h
  include/nori/scene.h
  include/nori/shape.h
  include/nori/texture.h
  include/nori/timer.h
  include/nori/transform.h
  include/nori/vector.h
  include/nori/warp.h
  include/nori/lodepng.h
  include/nori/medium.h
  include/nori/phasefunction.h
  include/nori/volume.h

  # Source code files
  src/bitmap.cpp
  src/block.cpp
  src/bvh.cpp
  src/chi2test.cpp
  src/common.cpp
  src/consttexture.cpp
  src/imagetexture.cpp
  src/checkerboard.cpp
  src/diffuse.cpp
  src/gui.cpp
  src/independent.cpp
  src/main.cpp
  src/mesh.cpp
  src/obj.cpp
  src/object.cpp
  src/parser.cpp
  src/perspective.cpp
  src/proplist.cpp
  src/render.cpp
  src/rfilter.cpp
  src/scene.cpp
  src/shape.cpp
  src/ttest.cpp
  src/warp.cpp
  src/microfacet.cpp
  src/photon.cpp
  src/mirror.cpp
  src/dielectric.cpp
  src/photonmapper.cpp
  src/ppm.cpp
  src/sphere.cpp
  src/arealight.cpp
  src/av.cpp
  src/direct.cpp
  src/direct_ems.cpp
  src/direct_mats.cpp
  src/direct_mis.cpp
  src/disney_diffuse.cpp
  src/disney_metallic.cpp
  src/disney_sheen.cpp
  src/normals.cpp
  src/path_mats.cpp
  src/path_mis.cpp
  src/pointlight.cpp
  src/lodepng.cpp
  src/volpath.cpp
  src/homogeneous.cpp
  src/henyeyGreenstein.cpp
  src/heterogeneous.cpp
  src/isotropic.cpp
  src/nullBSDF.cpp
  src/constantFloatVolume.cpp
  src/perlinVolume.cpp
  src/emissiveMedium.cpp
  src/grid.cpp
  src/perlinNoise3D.cpp
  src/volpath_emission.cpp
  src/heightVolume.cpp
  src/disney.cpp
)

add_executable(nori ${NORI_SOURCE_FILES})

# The following lines build the warping test application
add_executable(warptest
  include/nori/warp.h
  src/warp.cpp
  src/warptest.cpp
  src/microfacet.cpp
  src/object.cpp
  src/proplist.cpp
  src/common.cpp
)

target_link_libraries(nori ${EXTERNAL_LIBS})
target_link_libraries(warptest ${EXTERNAL_LIBS})

if (NORI_COMPILE_LIB)
  add_library(libnori ${NORI_SOURCE_FILES})
  target_link_libraries(libnori ${EXTERNAL_LIBS})
endif()

# Please do not change this as it is used for testing
if (NORI_HEADLESS)
  target_compile_definitions(nori PUBLIC NORI_HEADLESS)
endif()

# Force colored output for the ninja generator
if (CMAKE_GENERATOR STREQUAL "Ninja")
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fcolor-diagnostics")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
  endif()
endif()
