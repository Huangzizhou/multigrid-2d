cmake_minimum_required(VERSION 3.1)
project(multigrid)

# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")

### Compile Options ###
option(WITH_TBB             "Use tbb for parallel running"        OFF)

#######################

### Configuration
set(PROJECT_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(ProjectDependencies)

# Add project files
file(GLOB SRCFILES src/*.cpp)

add_executable(${PROJECT_NAME}_2D ${SRCFILES})

# TBB
if(WITH_TBB)
    download_tbb()
    set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
    set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)

    add_subdirectory(${PROJECT_EXTERNAL}/tbb tbb)
    set_property(TARGET tbb_static tbb_def_files PROPERTY FOLDER "dependencies")

    target_compile_definitions(tbb_static PUBLIC -DUSE_TBB)
    target_include_directories(${PROJECT_NAME}_2D SYSTEM PUBLIC ${PROJECT_EXTERNAL}/tbb/include)
    add_library(TBB::tbb ALIAS tbb_static)

    target_link_libraries(${PROJECT_NAME}_2D PUBLIC TBB::tbb)
endif()

# Eigen
if(NOT TARGET Eigen3::Eigen)
    download_eigen()
    add_library(eigen INTERFACE)
    target_include_directories(eigen SYSTEM INTERFACE
        $<BUILD_INTERFACE:${PROJECT_EXTERNAL}/eigen>
        $<INSTALL_INTERFACE:include>
    )
    set_property(TARGET eigen PROPERTY EXPORT_NAME Eigen3::Eigen)
    add_library(Eigen3::Eigen ALIAS eigen)
endif()
target_link_libraries(${PROJECT_NAME}_2D PUBLIC Eigen3::Eigen)