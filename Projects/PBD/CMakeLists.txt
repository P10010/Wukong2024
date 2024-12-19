if (APPLE)
    set(DEPS boost_filesystem
        igl::core
        suitesparseconfig metis cholmod amd camd ccolamd colamd
        gmp mpfr
        )
else()
    set(DEPS boost_filesystem
            igl::core
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
#    find_package(TBB REQUIRED)
    target_link_libraries(PBD ${CHOLMOD_LIBRARIES} TBB::tbb)
else()
    find_package(CHOLMOD REQUIRED)
    include_directories(${CHOLMOD_INCLUDES})
    target_link_libraries(PBD ${CHOLMOD_LIBRARIES})
endif()



target_include_directories(PBD PUBLIC ../../Libs/spectra/include)
target_include_directories(PBD PUBLIC ../../Deps/libigl/include)
target_link_libraries(PBD ${DEPS} polyscope)


target_include_directories(PBDLib PUBLIC ../../Libs/spectra/include)
target_include_directories(PBDLib PUBLIC ../../Deps/libigl/include)
target_link_libraries(PBDLib ${DEPS} polyscope)
