cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)
add_definitions(-Wall -O2)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

include_directories("../../../src")
include_directories("../../../src/ErrorCorrectionCode")
include_directories("../../../src/lib")
include_directories("../../../src/channel")
include_directories(".")
include_directories("../trace_learnt_path_selector/egolay_24/common")

set(common_src
        ../../../src/lib/libMath.cpp
        ../../../src/lib/libParser.cpp
        ../../../src/channel/channelBase.cpp
        ../../../src/channel/AWGN.cpp
        ../../../src/ErrorCorrectionCode/ErrorCorrectionCodeBase.cpp
        ../../../src/ErrorCorrectionCode/KernalManager.cpp
        ../../../src/ErrorCorrectionCode/BCH.cpp
        ../../../src/ErrorCorrectionCode/PolarMultiKernal_PS_PED.cpp
    )

add_executable(golay_mclass_bler_sim
        ../trace_learnt_path_selector/egolay_24/common/golay_mclass_bler_sim.cpp
        ${common_src}
    )

add_executable(golay_trace_dataset
        ../trace_learnt_path_selector/egolay_24/common/golay_trace_dataset.cpp
        ${common_src}
    )

add_executable(golay_adaptive_m_dataset
        ../adaptive_path_selector/egolay_24/golay_adaptive_m_dataset.cpp
        ${common_src}
    )
