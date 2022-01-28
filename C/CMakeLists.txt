cmake_minimum_required(VERSION 3.12.3)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR})
set(INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)

#Configure profiling for gprof.
set(PROFILING FALSE)
if (PROFILING)
   message(STATUS "Adding profiling capabilities.")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
   SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
   SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

#Configure Boost.
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost) 

if(NOT Boost_FOUND)
   message(ERROR "Could not find Boost.")
endif()

project(HasegawaMima)
message(STATUS "Configuring ${PROJECT_NAME}")

add_subdirectory(HM)