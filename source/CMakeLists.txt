########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
#CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################

### LIBRARY #################################################################

# definitions to pass to the compiler

# ADD_DEFINITIONS( "-pedantic -Wunused-value -O2" )
# ADD_DEFINITIONS("-Wno-long-long -Wreturn-type -Werror")

ADD_DEFINITIONS( "-Wall -ansi -pedantic -std=c++11 -O3" )
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
include_directories( ./include )


FILE(GLOB_RECURSE headers "include/*.h")
FILE(GLOB_RECURSE sources "src/*.cc")

#add_custom_target(headers SOURCES ${toto})

# add library
#AUX_SOURCE_DIRECTORY( ./src/CaloObject library_sources )
#AUX_SOURCE_DIRECTORY( ./src/Algorithm library_sources )
#ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
#INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )



#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

add_library( ${PROJECT_NAME} SHARED ${sources} ${headers} )


install( TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib )
install( DIRECTORY ./include DESTINATION . FILES_MATCHING PATTERN "*.h")
