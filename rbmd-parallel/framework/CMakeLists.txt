﻿# 依赖vtkm
cmake_minimum_required(VERSION 3.8)
project(framework)

#安装cxxopts
include(./cmake/cxxopts.cmake)
setup_cxxopts()

#安装json
include(./cmake/jsoncpp.cmake)
setup_jsoncpp()

if(NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
set(VTKm_DIR "${CMAKE_CURRENT_LIST_DIR}/../../vtkm_out_float/lib/cmake/vtkm-1.9" CACHE PATH "VTK directory override" FORCE)
endif()

find_package(VTKm REQUIRED QUIET)

find_package(MPI REQUIRED QUIET)
set(CMAKE_CXX_COMPILER mpicxx)
set(MPI_CXX_SKIP_MPICXX TRUE)
set(target ${project_name})

set(host_src
parser/Configuration.cpp

application/Application.cpp
application/MDApplication.cpp

#action/Action.cpp
#action/CreateSystemAction.cpp
#action/CreateExecutionerAction.cpp
#action/AddOutputAction.cpp
#action/SetupDeviceAction.cpp
#action/CreateInitConditionAction.cpp

executioner/Executioner.cpp

output/Output.cpp
output/FileOutput.cpp
output/ConsoleOutput.cpp
output/Progress.cpp

model/FormatTable.cpp
#model/Register.cpp

#system/System.cpp
InitCondition/MeshFreeCondition.cpp
InitCondition/ModelFileInitCondition.cpp
InitCondition/MeshFreeFileInitCondition.cpp
hyperparameters/Coulomb.cpp
hyperparameters/Extend.cpp
hyperparameters/ForceField.cpp
hyperparameters/Neighbor.cpp
)

set(device_src
#system/MadelungSystem.cpp
#system/LennardJonesSystem.cpp
#system/MeshFreeSystem.cpp
#system/SalineSolutionSystem.cpp # Remove the "<vtkm/cont/Algorithm.h>" to move out of Devie_src
#system/H2OSystem.cpp
#system/NaClSystem.cpp
#system/EAMSystem.cpp
#system/PEOSystem.cpp
output/ThermoOutput.cpp
output/TempOutput.cpp
output/RDFOutput.cpp
output/MSDOutput.cpp
output/VACFOutput.cpp
output/worklet/OutPutWorklet.cpp
output/TrajectoryOutput.cpp
#system/MDSystem.cpp
system/worklet/SystemWorklet.cpp

InitCondition/LJInitCondition.cpp
InitCondition/MadelungInitCondition.cpp
locator/ContPointLocator.cxx
forceFunction/ContForceFunction.cpp
topology/ContTopology.cpp
system/InitGlobal.cpp
staticTable/ContStaticTable.cpp
run/ExecutionNPT.cpp
run/ExecutionNVT.cpp
run/ExecutionNVE.cpp
run/ExecutionMD.cpp
run/worklet/RunWorklet.cpp
parallel/src/CollectiveCommunicatior/CommunicatorNonBlockingMPIComponents.cpp
parallel/src/CollectiveCommunicatior/CommunicatorNonBlockingScheduler.cpp
parallel/src/CollectiveCommunicatior/CommunicatorSequential.cpp
parallel/src/DomainDecomposition/BaseNeighbourCommunicationScheme.cpp
parallel/src/DomainDecomposition/CommunicationBuffer.cpp
parallel/src/DomainDecomposition/CommunicationPartner.cpp
parallel/src/DomainDecomposition/Domain.cpp
parallel/src/DomainDecomposition/DomainDecompBase.cpp
parallel/src/DomainDecomposition/DomainDecompositionMPI.cpp
parallel/src/DomainDecomposition/DomainDecompSequential.cpp
parallel/src/DomainDecomposition/DomianServiceLocator.cpp
parallel/src/DomainDecomposition/IndirectNeighbourCommunicationScheme.cpp
parallel/src/DomainDecomposition/NeighborAcquirer.cpp
parallel/src/MDataStructure/Atom.cpp
parallel/src/MDataStructure/Cell.cpp
parallel/src/MDataStructure/CellManager.cpp
parallel/src/MDataStructure/LinkedCell.cpp
parallel/src/ZonalMethods/ZonalMethod.cpp
parallel/src/RBMDParallelUtil.cpp
)

# 头文件引入
add_library(_compiler_flags INTERFACE)
target_include_directories(_compiler_flags INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/action
  ${CMAKE_CURRENT_SOURCE_DIR}/application
  ${CMAKE_CURRENT_SOURCE_DIR}/executioner
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh
  ${CMAKE_CURRENT_SOURCE_DIR}/meshmodifier
  ${CMAKE_CURRENT_SOURCE_DIR}/output
  ${CMAKE_CURRENT_SOURCE_DIR}/parser
  ${CMAKE_CURRENT_SOURCE_DIR}/system
  ${CMAKE_CURRENT_SOURCE_DIR}/userobject
  ${CMAKE_CURRENT_SOURCE_DIR}/model
  ${CMAKE_CURRENT_SOURCE_DIR}/InitCondition
  ${CMAKE_CURRENT_SOURCE_DIR}/locator
  ${CMAKE_CURRENT_SOURCE_DIR}/forceFunction
  ${CMAKE_CURRENT_SOURCE_DIR}/topology
  ${CMAKE_CURRENT_SOURCE_DIR}/hyperparameters
  ${CMAKE_CURRENT_SOURCE_DIR}/run
  ${JSON_INSTALL_PATH}/include
  ${CXXOPTS_INSTALL_PATH}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/parallel/include
  ${CMAKE_CURRENT_SOURCE_DIR}/parallel/include/CollectiveCommunicatior
  ${CMAKE_CURRENT_SOURCE_DIR}/parallel/include/DomainDecomposition
  ${CMAKE_CURRENT_SOURCE_DIR}/parallel/include/MDataStructure
  ${CMAKE_CURRENT_SOURCE_DIR}/parallel/include/ZonalMethods
  )
  # 看来要在res 文件下建立cmake文件
  # 定义lib
add_library(libsemd STATIC ${host_src} ${device_src})
target_link_libraries(libsemd PUBLIC _compiler_flags)

#链接json库
#set(CMAKE_PREFIX_PATH ${JSON_INSTALL_PATH}/lib/cmake)
#find_package(jsoncpp REQUIRED)
target_link_directories(libsemd PUBLIC ${JSON_INSTALL_PATH}/lib)
target_link_libraries(libsemd PRIVATE jsoncpp)

#加载 VTK-m 
target_link_libraries(libsemd PUBLIC vtkm_cont vtkm_worklet vtkm_io)
target_link_libraries(libsemd PRIVATE ${VTK_LIBRARIES})
target_link_libraries(libsemd PUBLIC fmt::fmt)

target_link_libraries(libsemd PRIVATE MPI::MPI_CXX)

vtkm_add_target_information(libsemd DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS DEVICE_SOURCES ${device_src})

  # 定义target
add_executable(${target} main.cpp)
#file(COPY ${CMAKE_SOURCE_DIR}/framework/res/Gnear_table_geq2_2.txt DESTINATION ${CMAKE_BINARY_DIR}/framework/res)
#file(COPY ${CMAKE_SOURCE_DIR}/framework/res/Gnear_table_leq2_2.txt DESTINATION ${CMAKE_BINARY_DIR}/framework/res)
#file(COPY "../res" DESTINATION ${CMAKE_BINARY_DIR}/framework)
#file(COPY "../res" DESTINATION ${CMAKE_BINARY_DIR}/framework)



target_link_libraries(${target} PUBLIC libsemd ${MPI_CXX_LIBRARIES})
vtkm_add_target_information(${target} DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS DEVICE_SOURCES main.cpp)

set_target_properties(${target} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set_target_properties(${target} PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

install(TARGETS ${target} 
  DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
  )

#install(DIRECTORY "../res"
#        DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
        #FILES_MATCHING
        #PATTERN "Gnear_table_geq2_2.txt"
        #PATTERN "Gnear_table_leq2_2.txt")

#install(DIRECTORY "../res/"
#        DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
#        FILES_MATCHING
#        PATTERN "LJFluid.json"
#        PATTERN "rbmd.data"
#        PATTERN "rbmd.json")