cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(km3sim)

find_package(Geant4 REQUIRED)

include(${Geant4_USE_FILE})

include(ExternalProject)
ExternalProject_Add(docopt.cpp
  PREFIX ${PROJECT_SOURCE_DIR}/external/
  GIT_REPOSITORY https://github.com/docopt/docopt.cpp.git
  LOG_DOWNLOAD ON
  LOG_CONFIGURE ON
  LOG_BUILD ON
  LOG_INSTALL ON
)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/external/seaweed)
include_directories(${PROJECT_SOURCE_DIR}/external/docopt.cpp)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/src/*.h)

add_executable(km3sim km3sim.cc ${sources} ${headers})
target_link_libraries(km3sim ${Geant4_LIBRARIES})

set(KM3SIM_SCRIPTS
	run_km3sim.sh
)

#foreach(_script ${EXAMPLEB1_SCRIPTS})
#	configure_file(
#		${PROJECT_SOURCE_DIR}/${_script}
#		${PROJECT_BINARY_DIR}/${_script}
#		COPYONLY
#	)
#endforeach()

#install(TARGETS DESTINATION bin)
