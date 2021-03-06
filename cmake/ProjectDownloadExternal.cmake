include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
  set(PROJECT_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
  set(PROJECT_EXTRA_OPTIONS "")
endif()

function(custom_download_project name)
  download_project(
    PROJ         ${name}
    SOURCE_DIR   ${PROJECT_EXTERNAL}/${name}
    DOWNLOAD_DIR ${PROJECT_EXTERNAL}/.cache/${name}
    QUIET
    ${PROJECT_EXTRA_OPTIONS}
    ${ARGN}
  )
endfunction()

################################################################################

# libigl
function(download_libigl)
  custom_download_project(libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG        b0d7740e0b7e887a7e93601c4c557ecf762b389b
  )
endfunction()

# TBB
function(download_tbb)
  custom_download_project(tbb
    # GIT_REPOSITORY https://github.com/intel/tbb.git
    # GIT_TAG        2018_U5
    GIT_REPOSITORY https://github.com/wjakob/tbb.git
    GIT_TAG        344fa84f34089681732a54f5def93a30a3056ab9
  )
endfunction()

# AMGCL
function(download_amgcl)
    custom_download_project(amgcl
       GIT_REPOSITORY https://github.com/ddemidov/amgcl.git
       GIT_TAG        461a66ce6d197a3816218bf94ffd114a367c1ef1
    )
endfunction()

# EIGEN
set(LIBIGL_EIGEN_VERSION 3.3.7 CACHE STRING "Default version of Eigen.")
function(download_eigen)
	  custom_download_project(eigen
		    GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
		    GIT_TAG        ${LIBIGL_EIGEN_VERSION}
	  )
endfunction()
