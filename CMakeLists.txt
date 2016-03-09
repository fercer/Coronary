MESSAGE(CMAKE_VERSION: ${CMAKE_VERSION})
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(NOMBRE_PROYECTO Coronary)
SET(HEADERS src/reconstructor_3D.h src/filtros.h src/IMGVTK.h src/pars_fercer.h)
SET(SOURCES src/main.cpp src/reconstructor_3D.cpp src/filtros.cpp src/IMGVTK.cpp src/pars_fercer.c)

SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3 -DNDEBUG -lm -fopenmp")
SET(CMAKE_C_FLAGS_DEBUG "-g -lm -fopenmp")
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -DNDEBUG -fopenmp")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp")

PROJECT( ${NOMBRE_PROYECTO} )
ADD_DEFINITIONS(-std=c99)
ADD_EXECUTABLE( ${NOMBRE_PROYECTO} ${SOURCES} ${HEADERS})

find_path (FFTW_INCLUDES fftw3.h)
find_library (FFTW_LIBRARIES NAMES fftw3)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

TARGET_LINK_LIBRARIES(${NOMBRE_PROYECTO} ${FFTW_LIBRARIES} ) 

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
IF(VTK_LIBRARIES)
  MESSAGE("Librerias VTK encontradas en: ${VTK_DIR}")
  TARGET_LINK_LIBRARIES(${NOMBRE_PROYECTO} ${VTK_LIBRARIES})
ELSE(VTK_LIBRARIES)
  MESSAGE(ERROR:"Librerias VTK no encontradas")
ENDIF(VTK_LIBRARIES)

FIND_PACKAGE(GDCM REQUIRED)
IF(GDCM_FOUND)
  INCLUDE(${GDCM_USE_FILE})
  MESSAGE("Librerias GDCM encontradas en: ${GDCM_DIR}")
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
     socketxx
     vtkgdcm)
#     Rpcrt4)
  TARGET_LINK_LIBRARIES(${NOMBRE_PROYECTO} ${GDCM_LIBRARIES})
ELSE(GDCM_FOUND)
  MESSAGE(FATAL_ERROR "Cannot find GDCM, did you set GDCM_DIR?")
ENDIF(GDCM_FOUND)