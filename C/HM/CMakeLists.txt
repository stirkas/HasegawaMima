set(SUBPROJECT_NAME "HM")
message(STATUS "Configuring ${SUBPROJECT_NAME}")

file(GLOB_RECURSE SRC_FILES *.cpp *.c)
file(GLOB_RECURSE HDR_FILES *.hpp *.h)

add_executable(${SUBPROJECT_NAME} ${HDR_FILES} ${SRC_FILES})

#target_include_directories(${SUBPROJECT_NAME} PRIVATE ${BLAS_INC_DIR} ${LAPACKE_INC_DIR})

target_link_libraries(${SUBPROJECT_NAME} PRIVATE fftw3)
#target_link_libraries(${SUBPROJECT_NAME} PRIVATE ${LAPACK_LIBS} fftw3)

install(TARGETS ${SUBPROJECT_NAME} DESTINATION ${INSTALL_DIR})