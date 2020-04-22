cmake_minimum_required(VERSION 3.0)
project(GLIKE VERSION 0.8.0)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)

# load ROOT libraries, check that minuit is installed
set(ROOT_REQUIRED_COMPONENTS Minuit)
find_package(ROOT REQUIRED ${ROOT_REQUIRED_COMPONENTS})

# access root macros
include_directories(${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include)

# generate glike libraries
root_generate_dictionary(
    glikeDict
    ${PROJECT_SOURCE_DIR}/include/Lkl.h
    ${PROJECT_SOURCE_DIR}/include/ParabolaLkl.h
    ${PROJECT_SOURCE_DIR}/include/PoissonLkl.h
    ${PROJECT_SOURCE_DIR}/include/JointLkl.h
    ${PROJECT_SOURCE_DIR}/include/Iact1dUnbinnedLkl.h
    ${PROJECT_SOURCE_DIR}/include/Iact1dBinnedLkl.h
    ${PROJECT_SOURCE_DIR}/include/IactEventListIrf.h
    ${PROJECT_SOURCE_DIR}/include/MIACTEventListIRF.h
    ${PROJECT_SOURCE_DIR}/include/FermiTables2016Lkl.h
    ${PROJECT_SOURCE_DIR}/include/GloryDuckTables2019Lkl.h
    ${PROJECT_SOURCE_DIR}/include/TemplateLkl.h
    LINKDEF ${PROJECT_SOURCE_DIR}/include/gLikeLinkDef.h
)
add_library(
    glike SHARED 
    ${PROJECT_SOURCE_DIR}/src/Lkl.cc
    ${PROJECT_SOURCE_DIR}/src/ParabolaLkl.cc
    ${PROJECT_SOURCE_DIR}/src/PoissonLkl.cc
    ${PROJECT_SOURCE_DIR}/src/JointLkl.cc
    ${PROJECT_SOURCE_DIR}/src/Iact1dUnbinnedLkl.cc
    ${PROJECT_SOURCE_DIR}/src/Iact1dBinnedLkl.cc
    ${PROJECT_SOURCE_DIR}/src/IactEventListIrf.cc
    ${PROJECT_SOURCE_DIR}/src/MIACTEventListIRF.cc
    ${PROJECT_SOURCE_DIR}/src/FermiTables2016Lkl.cc
    ${PROJECT_SOURCE_DIR}/src/GloryDuckTables2019Lkl.cc
    ${PROJECT_SOURCE_DIR}/src/TemplateLkl.cc
    glikeDict
)
target_link_libraries(glike ${ROOT_LIBRARIES})

# generate executables from the scripts
set(
    EXELIST 
    testPoissonLkl
    testJointLkl
    testDMSearches
)

foreach(exe ${EXELIST})
  add_executable(${exe} ${PROJECT_SOURCE_DIR}/scripts/${exe}.C)
  target_link_libraries (${exe} glike ${ROOT_LIBRARIES})
endforeach()

# install destination
install(
    TARGETS glike ${EXELIST}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/static
)
