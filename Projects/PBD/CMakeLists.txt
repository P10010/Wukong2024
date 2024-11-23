if (APPLE)
    set(DEPS boost_filesystem tbb 
        tbbmalloc tbbmalloc_proxy 
        igl::core
        suitesparseconfig metis cholmod amd camd ccolamd colamd
        gmp mpfr
        )
else()
    set(DEPS boost_filesystem tbb 
        tbbmalloc tbbmalloc_proxy 
        mkl_intel_lp64 igl::core
        mkl_sequential mkl_core 
        suitesparseconfig metis cholmod amd camd ccolamd colamd
        gmp mpfr
        )
endif()


file(GLOB HEADERS "include/*.h")
file(GLOB SOURCES "src/*.cpp" )


add_executable(PBD ${HEADERS} ${SOURCES})
add_library(PBDLib STATIC ${HEADERS} ${SOURCES})

if (APPLE)
    find_package(CHOLMOD REQUIRED)
    include_directories(${CHOLMOD_INCLUDES})
    find_package(TBB REQUIRED)
    target_link_libraries(PBD ${CHOLMOD_LIBRARIES} TBB::tbb)
endif()



target_include_directories(PBD PUBLIC ../../Libs/spectra/include)
target_include_directories(PBD PUBLIC ../../Deps/libigl/include)
target_link_libraries(PBD ${DEPS} polyscope)


target_include_directories(PBDLib PUBLIC ../../Libs/spectra/include)
target_include_directories(PBDLib PUBLIC ../../Deps/libigl/include)
target_link_libraries(PBDLib ${DEPS} polyscope)
