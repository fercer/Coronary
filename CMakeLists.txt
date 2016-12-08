MESSAGE(CMAKE_VERSION: ${CMAKE_VERSION})
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

if(POLICY CMP0020)
  cmake_policy(SET CMP0020 NEW)
endif()



SET(SRC_DIR "./src")
SET(HEADERS ${SRC_DIR}/reconstructor_3D.h ${SRC_DIR}/filtros.h ${SRC_DIR}/IMGVTK.h ${SRC_DIR}/args_fercer.h)
SET(SOURCES ${SRC_DIR}/reconstructor_3D.cpp ${SRC_DIR}/filtros.cpp ${SRC_DIR}/IMGVTK.cpp ${SRC_DIR}/args_fercer.cpp)

# The GUI version is going to be build?
OPTION(BUILD_GUI "Build gui version." OFF)
OPTION(BUILD_VTK "Build vtk version." OFF)
IF(BUILD_GUI)
    SET(NOMBRE_PROYECTO CoronaryGUI)
    
ELSE(BUILD_GUI)
    SET(NOMBRE_PROYECTO Coronary)
    
    SET(SOURCES ${SOURCES} ${SRC_DIR}/main.cpp)
ENDIF(BUILD_GUI)


SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3 -DNDEBUG -lm -fopenmp -std=c99")
SET(CMAKE_C_FLAGS_DEBUG "-g -lm -fopenmp -std=c99")
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -DNDEBUG -fopenmp")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp")

# Project name for building
PROJECT( ${NOMBRE_PROYECTO} )


IF(BUILD_GUI)
    ADD_DEFINITIONS(-DBUILD_GUI_VERSION)
    
    FILE(GLOB UI_FILES *.ui)
    FILE(GLOB QT_WRAP *.h)
    FILE(GLOB CXX_FILES *.cpp)
ENDIF(BUILD_GUI)


FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
IF(VTK_LIBRARIES AND BUILD_VTK)
    MESSAGE("VTK libraries found in: ${VTK_DIR}")
    ADD_DEFINITIONS(-DBUILD_VTK_VERSION)
    IF(${VTK_VERSION} VERSION_GREATER "6" AND VTK_QT_VERSION VERSION_GREATER "4")
        IF(BUILD_GUI)
            # Instruct CMake to run moc automatically when needed.
            SET(CMAKE_AUTOMOC ON)
            FIND_PACKAGE(Qt5Widgets REQUIRED QUIET)
            QT5_WRAP_UI(UISrcs ${UI_FILES} )
            # CMAKE_AUTOMOC in ON so the MOC headers will be automatically wrapped.
            QT5_USE_MODULES(${NOMBRE_PROYECTO} Core Gui)
        ENDIF(BUILD_GUI)
        
        # Add the sources to the project:
        ADD_EXECUTABLE(${NOMBRE_PROYECTO} MACOSX_BUNDLE ${CXX_FILES} ${UISrcs} ${QT_WRAP} ${SOURCES} ${HEADERS})
        
        TARGET_LINK_LIBRARIES(${NOMBRE_PROYECTO} ${VTK_LIBRARIES})
    ELSE()
        IF(BUILD_GUI)
            FIND_PACKAGE(Qt4 REQUIRED)
            INCLUDE(${QT_USE_FILE})
            QT4_WRAP_UI(UISrcs ${UI_FILES})
            QT4_WRAP_CPP(MOCSrcs ${QT_WRAP})
        ENDIF(BUILD_GUI)
        
        # Add the sources to the project:
        ADD_EXECUTABLE(${NOMBRE_PROYECTO} MACOSX_BUNDLE ${CXX_FILES} ${UISrcs} ${MOCSrcs} ${SOURCES} ${HEADERS})
        
        IF(${VTK_VERSION} VERSION_LESS "6")
            TARGET_LINK_LIBRARIES(${NOMBRE_PROYECTO} ${VTK_LIBRARIES} QVTK)
        ELSE()
            TARGET_LINK_LIBRARIES(${NOMBRE_PROYECTO} ${VTK_LIBRARIES})
        ENDIF()
    ENDIF()
ELSE(VTK_LIBRARIES AND BUILD_VTK)
    MESSAGE("VTK libraries could not be found, a No-VTK version will be built")
    # Add the sources to the project:
    ADD_EXECUTABLE(${NOMBRE_PROYECTO} MACOSX_BUNDLE ${CXX_FILES} ${UISrcs} ${QT_WRAP} ${SOURCES} ${HEADERS})
ENDIF(VTK_LIBRARIES AND BUILD_VTK)

# Add FFTW3 libraries:
FIND_PATH (FFTW_INCLUDES fftw3.h)
SET(FFTW_LIBRARIES "${FFTW_INCLUDES}/libfftw3-3.lib;${FFTW_INCLUDES}/libfftw3l-3.lib;${FFTW_INCLUDES}/libfftw3f-3.lib")

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)
ADD_DEFINITIONS(-DMY_FFTW_PATH="${FFTW_INCLUDES}/fftw3.h")
TARGET_LINK_LIBRARIES(${NOMBRE_PROYECTO} ${FFTW_LIBRARIES} ) 


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})



OPTION(BUILD_GDCM "Build GDCM version." OFF)
IF(BUILD_GDCM)
    FIND_PACKAGE(GDCM REQUIRED)
    IF(GDCM_FOUND)
        ADD_DEFINITIONS(-DBUILD_GDCM_VERSION)
        INCLUDE(${GDCM_USE_FILE})
        MESSAGE("GDCM libraries found in: ${GDCM_DIR}")
        SET(GDCM_LIBRARIES 
        #     gdcmcharls
            gdcmCommon
            gdcmDICT		
            gdcmDSED
        #     gdcmexpat
        #     gdcmgetopt
            gdcmIOD
            gdcmjpeg12
            gdcmjpeg16
            gdcmjpeg8
            gdcmMEXD
            gdcmMSFF
        #     gdcmopenjpeg
        #     gdcmzlib
            socketxx)
        #     vtkgdcm)
        #     Rpcrt4)
        TARGET_LINK_LIBRARIES(${NOMBRE_PROYECTO} ${GDCM_LIBRARIES})
    ELSE(GDCM_FOUND)
        MESSAGE(FATAL_ERROR "Cannot find GDCM, did you set GDCM_DIR?")
    ENDIF(GDCM_FOUND)    
ENDIF(BUILD_GDCM)
