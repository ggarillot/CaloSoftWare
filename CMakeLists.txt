###############################################
# cmake file for building CaloSoftWare package
# @author Eté Rémi
# Copyright (c) CNRS / IPNL
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
#################################################

# run in build directory with : 
# cmake -DCLHEP_DIR=$CLHEP ../

# project name
PROJECT( CaloSoftWare )

SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

###############################

### DEPENDENCIES ############################################################
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )
INCLUDE( default_settings )
#MESSAGE( "PROJECT_SOURCE_DIR/cmake : ${PROJECT_SOURCE_DIR}/cmake")
FIND_PACKAGE( CLHEP REQUIRED )
FIND_PACKAGE( Eigen3 REQUIRED )

FOREACH( pkg CLHEP )
    #IF( ${pkg}_FOUND )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
ENDFOREACH()

INCLUDE_DIRECTORIES( ${EIGEN3_INCLUDE_DIR} )

# Define options

# ----- Install doc -----
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )
IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()

###############################
#Generate CaloSoftWareConfig.h 
# CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.h.cmake.in"
#   "${PROJECT_SOURCE_DIR}/source/include/${PROJECT_NAME}Config.h" @ONLY )
    
###############################
# Add source directory
ADD_SUBDIRECTORY( ./source )

###############################
# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

###############################
# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( ${PROJECT_NAME}Config.cmake ${PROJECT_NAME}ConfigVersion.cmake ${PROJECT_NAME}LibDeps.cmake )