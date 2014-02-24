
cmake_minimum_required(VERSION 2.8.11)

project(SSPipeline)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

include(${CMAKE_CURRENT_SOURCE_DIR}/Qt4GenerateMocs.cmake)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")

find_package(Slicer COMPONENTS NoExtensionSupport REQUIRED)
include(${VTK_USE_FILE})
include(${CTK_USE_FILE})
include(${Teem_USE_FILE})

find_package(Qt REQUIRED)
include(${QT_USE_FILE})

include_directories(
  ${vtkITK_INCLUDE_DIRS}
  ${MRMLCore_INCLUDE_DIRS}
  )

QT4_GENERATE_MOCS(main.cpp)

add_executable(ssp
  main.cpp
  )
target_link_libraries(ssp
  vtkITK
  MRMLCore
  ${PYTHON_LIBRARIES}
  ${QT_LIBRARIES}
  ${CTK_LIBRARIES}
  ${VTK_LIBRARIES}
  )
