cmake_minimum_required(VERSION 3.0)
project(GLIKE VERSION 0.8.0)

# where to save the .so (.dylib) and .pcm files
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib) 
# where to save the executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# load ROOT libraries, define a cache variable to enable FITSIO support 
# if USE_FITSIO is selected, Fitsio is added to the ROOT requirements   
set(ROOT_REQUIRED_COMPONENTS Minuit)
option(USE_FITSIO "Enable support for FITS files I/O" OFF)                      
if(USE_FITSIO)                     
    list(APPEND ROOT_REQUIRED_COMPONENTS Fitsio)
endif()   
find_package(ROOT REQUIRED ${ROOT_REQUIRED_COMPONENTS})
message("-- Required ROOT components are: ${ROOT_REQUIRED_COMPONENTS}")

# generate glike ROOT.pcm dictionaries
root_generate_dictionary(
    gLikeDict
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

# create a shared library with the generated dictionary
add_library(
    gLike SHARED 
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
    gLikeDict
)
# link gLike libraries against ROOT's ones
target_link_libraries(gLike ${ROOT_LIBRARIES})
# where to find the header files for the gLike libraries
target_include_directories(gLike PUBLIC ${PROJECT_SOURCE_DIR}/include)
# install the gLike libraries in /usr/local/lib
install(TARGETS gLike DESTINATION lib)

