cmake_minimum_required(VERSION 3.12.3)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR})
set(INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)

project(HasegawaMima)
message(STATUS "Configuring ${PROJECT_NAME}")

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost) 

if(NOT Boost_FOUND)
   message(ERROR "Could not find Boost.")
endif()

add_subdirectory(HM)
