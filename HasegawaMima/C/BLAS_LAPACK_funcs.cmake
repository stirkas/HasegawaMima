function(Load_BLAS_LAPACK)

   #CMake has a good interface for BLAS and LAPACK - find_package(BLAS), etc. -
   #but not for LAPACKE (C interface to fortran LAPACK) so I'll just do them all manually.
   if(DEFINED ENV{BLAS_INC} AND EXISTS $ENV{BLAS_INC})
      set(BLAS_INCLUDE_DIR $ENV{BLAS_INC} PARENT_SCOPE)
      message(STATUS "BLAS_INCLUDE_DIR found: $ENV{BLAS_INC}")
   else()
      message(SEND_ERROR "BLAS_INC env var needs to be changed in CMake, env var not defined, or path does not exist.")
   endif()
   if(DEFINED ENV{LAPACKE_INC} AND EXISTS $ENV{LAPACKE_INC})
      set(LAPACKE_INC_DIR $ENV{LAPACKE_INC} PARENT_SCOPE)
      message(STATUS "LAPACKE_INC_DIR found: $ENV{LAPACKE_INC}")
   else()
      message(SEND_ERROR "LAPACKE_INC env var needs to be changed in CMake, env var not defined, or path does not exist.")
   endif()

   #Libs handled differently. Glob files from dir and set parent scope var to files not dir.
   #Linux installs BLAS, LAPACK, and LAPACKE to same place. Otherwise need to separate logic...
   #TODO: Test for existence of static libs since not globbing, and also test for existence of fortran compiler (.so). (.a requires dependencies...)
   if(DEFINED ENV{LAPACK_LIB} AND EXISTS $ENV{LAPACK_LIB})
      set(LAPACK_LIB_DIR $ENV{LAPACK_LIB})

      set(LAPACK_STATIC_LIBS "${LAPACK_LIB_DIR}/liblapacke.a" "${LAPACK_LIB_DIR}/liblapack.a" "${LAPACK_LIB_DIR}/libblas.a" "-lgfortran")
      set(LAPACK_STATIC_LIBS ${LAPACK_STATIC_LIBS} PARENT_SCOPE)
      file(GLOB LAPACK_DYN_LIBS "${LAPACK_LIB_DIR}/libblas.so" "${LAPACK_LIB_DIR}/liblapack.so" "${LAPACK_LIB_DIR}/liblapacke.so")
      set(LAPACK_LIBS ${LAPACK_DYN_LIBS} PARENT_SCOPE)

      message(STATUS "Following static LAPACK libs found: ${LAPACK_STATIC_LIBS}")
      message(STATUS "Following dynamic LAPACK libs found: ${LAPACK_DYN_LIBS}")
   else()
      message(SEND_ERROR "LAPACK_LIB env var needs to be changed in CMake, env var not defined, or path does not exist.")
   endif()

endfunction()