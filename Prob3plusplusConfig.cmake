if(NOT TARGET Prob3plusplus::All)

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
# This will define the following variables
#
#    Prob3plusplus_FOUND
#
# and the following imported targets
#
#    Prob3plusplus::All
#

# Colorful messaging facilities
if(NOT COMMAND cmessage)
  if(NOT WIN32)
    string(ASCII 27 Esc)
    set(CM_ColourReset "${Esc}[m")
    set(CM_ColourBold "${Esc}[1m")
    set(CM_Red "${Esc}[31m")
    set(CM_Green "${Esc}[32m")
    set(CM_Yellow "${Esc}[33m")
    set(CM_Blue "${Esc}[34m")
    set(CM_Magenta "${Esc}[35m")
    set(CM_Cyan "${Esc}[36m")
    set(CM_White "${Esc}[37m")
    set(CM_BoldRed "${Esc}[1;31m")
    set(CM_BoldGreen "${Esc}[1;32m")
    set(CM_BoldYellow "${Esc}[1;33m")
    set(CM_BoldBlue "${Esc}[1;34m")
    set(CM_BoldMagenta "${Esc}[1;35m")
    set(CM_BoldCyan "${Esc}[1;36m")
    set(CM_BoldWhite "${Esc}[1;37m")
  endif()

  message(STATUS "Setting up colored messages...")

  function(cmessage)
    list(GET ARGV 0 MessageType)
    if(MessageType STREQUAL FATAL_ERROR OR MessageType STREQUAL SEND_ERROR)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_BoldRed}${ARGV}${CM_ColourReset}")
    elseif(MessageType STREQUAL WARNING)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_BoldYellow}${ARGV}${CM_ColourReset}")
    elseif(MessageType STREQUAL AUTHOR_WARNING)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_BoldCyan}${ARGV}${CM_ColourReset}")
    elseif(MessageType STREQUAL STATUS)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_Green}[INFO]:${CM_ColourReset} ${ARGV}")
    elseif(MessageType STREQUAL CACHE)        
      list(REMOVE_AT ARGV 0)
      message(-- "${CM_Blue}[CACHE]:${CM_ColourReset} ${ARGV}")
    elseif(MessageType STREQUAL DEBUG)
      list(REMOVE_AT ARGV 0)
      if(BUILD_DEBUG_MSGS)
        message("${CM_Magenta}[DEBUG]:${CM_ColourReset} ${ARGV}")
      endif()
    else()
      message(${MessageType} "${CM_Green}[INFO]:${CM_ColourReset} ${ARGV}")
    endif()
  endfunction()
endif()

get_filename_component(Prob3plusplus_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

find_path(Prob3plusplus_INCLUDE_DIR
  NAMES mosc3.h
  PATHS ${Prob3plusplus_CMAKE_DIR}
)

#find the library name
file(GLOB Prob3plusplus_LIBRARY_FQPATHS ${Prob3plusplus_CMAKE_DIR}/*.a)
list(GET Prob3plusplus_LIBRARY_FQPATHS 0 Prob3plusplus_LIBRARY_FQPATH)
get_filename_component(Prob3plusplus_LIBRARY 
  "${Prob3plusplus_LIBRARY_FQPATH}" NAME)

get_filename_component(Prob3plusplus_LIBRARY_DIR 
  "${Prob3plusplus_LIBRARY_FQPATH}" DIRECTORY)

string(REGEX REPLACE libThreeProb_\(.*\)\\.a \\1
       Prob3plusplus_VERSION ${Prob3plusplus_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Prob3plusplus
    REQUIRED_VARS 
      Prob3plusplus_INCLUDE_DIR
      Prob3plusplus_LIBRARY_DIR
      Prob3plusplus_LIBRARY
    VERSION_VAR Prob3plusplus_VERSION
)

if(Prob3plusplus_FOUND)

  cmessage(STATUS "Prob3plusplus Found: ${Prob3plusplus_CMAKE_DIR} ")
  cmessage(STATUS "    Prob3plusplus_INCLUDE_DIR: ${Prob3plusplus_INCLUDE_DIR}")
  cmessage(STATUS "    Prob3plusplus_LIBRARY: ${Prob3plusplus_LIBRARY}")
  cmessage(STATUS "    Prob3plusplus_VERSION: ${Prob3plusplus_VERSION}")

  if(NOT TARGET Prob3plusplus::All)
      add_library(Prob3plusplus::All INTERFACE IMPORTED)
      set_target_properties(Prob3plusplus::All PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${Prob3plusplus_INCLUDE_DIRS}"
          INTERFACE_LINK_DIRECTORIES "${Prob3plusplus_LIBRARY_DIR}"
          INTERFACE_LINK_LIBRARIES "${Prob3plusplus_LIBRARY}"
      )
  endif()

endif()

endif() # Only run this whole file if it hasn't been run already.