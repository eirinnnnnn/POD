cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)
add_definitions(-Wall -O2)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

get_filename_component(project_name ${CMAKE_CURRENT_LIST_DIR} NAME)

include_directories("../../src")
include_directories("../../src/ErrorCorrectionCode")
include_directories("../../src/lib")
include_directories("../../src/channel")

# ================================================================
# Common sources shared by normal simulator and AED inspector
# IMPORTANT: no main.cpp here, no AED_SC_inspect.cpp here
# ================================================================
set(common_src
        ../../src/lib/libMath.cpp
        ../../src/lib/libParser.cpp
        ../../src/channel/channelBase.cpp
        ../../src/channel/AWGN.cpp
        ../../src/channel/BSC.cpp
        ../../src/ErrorCorrectionCode/ErrorCorrectionCodeBase.cpp
        ../../src/ErrorCorrectionCode/AED_attemp_20260102.cpp
    )

set(common_header
        ../../src/lib/libDebug.h
        ../../src/lib/libMath.h
        ../../src/lib/libParser.h
        ../../src/channel/channelBase.h
        ../../src/channel/AWGN.h
        ../../src/channel/BSC.h
        ../../src/ErrorCorrectionCode/ErrorCorrectionCodeBase.h
        ../../src/ErrorCorrectionCode/AED_attemp_20260102.h
        ../../src/ErrorCorrectionCode/AED_relation_check.h
    )

# ================================================================
# Original executable
# ================================================================
add_executable(${project_name}
        main.cpp
        ${common_src}
        ${common_header}
    )

# ================================================================
# Diagnostic executable
# ================================================================
add_executable(AED_SC_inspect
        ../../src/ErrorCorrectionCode/AED_SC_inspect.cpp
        ${common_src}
        ${common_header}
    )

# ================================================================
# BSC (Bernoulli bit-flip channel) executable
# ================================================================
add_executable(${project_name}_BSC
        main_BSC.cpp
        ${common_src}
        ${common_header}
    )