cmake_minimum_required(VERSION 3.0)
project(GLIKE VERSION 0.8.0)

# sane defaults for CMAKE_INSTALL_*DIR
include(GNUInstallDirs)
# RPATH handling. This makes setting LD_LIBRARY_PATH unecessary
# as it will compile the path to the libraries into the libraries / executables
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@executable_path/../${CMAKE_INSTALL_LIBDIR}")
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# make layout in the build directory the same as install directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# load ROOT libraries, define a cache variable to enable FITSIO support 
# if USE_FITSIO is selected, Fitsio is added to the ROOT requirements   
set(ROOT_REQUIRED_COMPONENTS Minuit)
option(USE_FITSIO "Enable support for FITS files I/O" OFF)                      
if(USE_FITSIO)                     
    list(APPEND ROOT_REQUIRED_COMPONENTS FITSIO)
endif()   
find_package(ROOT REQUIRED ${ROOT_REQUIRED_COMPONENTS})
message("-- Required ROOT components are: ${ROOT_REQUIRED_COMPONENTS}")
# for ROOT versions < 6
include(${ROOT_USE_FILE})

# generate glike ROOT .pcm dictionaries
root_generate_dictionary(
    gLikeDict
    ${PROJECT_SOURCE_DIR}/include/Lkl.h
    ${PROJECT_SOURCE_DIR}/include/ParabolaLkl.h
    ${PROJECT_SOURCE_DIR}/include/PoissonLkl.h
    ${PROJECT_SOURCE_DIR}/include/JointLkl.h
    ${PROJECT_SOURCE_DIR}/include/Iact1dUnbinnedLkl.h
    ${PROJECT_SOURCE_DIR}/include/Iact1dBinnedLkl.h
    ${PROJECT_SOURCE_DIR}/include/IactEventListIrf.h
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
    ${PROJECT_SOURCE_DIR}/src/FermiTables2016Lkl.cc
    ${PROJECT_SOURCE_DIR}/src/GloryDuckTables2019Lkl.cc
    ${PROJECT_SOURCE_DIR}/src/TemplateLkl.cc
    gLikeDict
)
# make gLike headers public, for usage by other codes
set(GLIKE_HEADERS
    ${PROJECT_SOURCE_DIR}/include/Lkl.h
    ${PROJECT_SOURCE_DIR}/include/ParabolaLkl.h
    ${PROJECT_SOURCE_DIR}/include/PoissonLkl.h
    ${PROJECT_SOURCE_DIR}/include/JointLkl.h
    ${PROJECT_SOURCE_DIR}/include/Iact1dUnbinnedLkl.h
    ${PROJECT_SOURCE_DIR}/include/Iact1dBinnedLkl.h
    ${PROJECT_SOURCE_DIR}/include/IactEventListIrf.h
    ${PROJECT_SOURCE_DIR}/include/FermiTables2016Lkl.h
    ${PROJECT_SOURCE_DIR}/include/GloryDuckTables2019Lkl.h
    ${PROJECT_SOURCE_DIR}/include/TemplateLkl.h
)
set_target_properties(gLike PROPERTIES PUBLIC_HEADER "${GLIKE_HEADERS}")
# link gLike libraries against ROOT's ones
target_link_libraries(gLike ${ROOT_LIBRARIES})
# where to find the header files for the gLike libraries
target_include_directories(gLike PUBLIC ${PROJECT_SOURCE_DIR}/include)
# install the gLike libraries in /usr/local/lib
install(
    TARGETS gLike
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# copy scripts, rcfiles, DM and data directories in the build directory
add_custom_target(copy_dirs ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/scripts ${CMAKE_BINARY_DIR}/scripts
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/rcfiles ${CMAKE_BINARY_DIR}/rcfiles
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/DM ${CMAKE_BINARY_DIR}/DM
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/.rootrc ${CMAKE_BINARY_DIR}/.rootrc
)

# create jointLklDM executable
add_executable(jointLklDM ${PROJECT_SOURCE_DIR}/exec/jointLklDM.cc)
target_link_libraries(jointLklDM gLike)
install(TARGETS jointLklDM DESTINATION bin)

# create computeCLBands executable
add_executable(computeCLBands ${PROJECT_SOURCE_DIR}/exec/computeCLBands.cc)
target_link_libraries(computeCLBands gLike)
install(TARGETS computeCLBands DESTINATION bin)

# if requested by the user, generate html documentation
option(MAKE_DOCS "Generate html documentation" OFF)    
if(MAKE_DOCS)                     
    add_custom_target(html_doc ALL 
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/htmldoc
        COMMAND root -b -q ${CMAKE_BINARY_DIR}/scripts/dohtml.C 2>&1 > ${CMAKE_BINARY_DIR}/htmldoc/dohtml.log
        DEPENDS gLike copy_dirs
    )
endif() 
