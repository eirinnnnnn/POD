cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)
add_definitions(-Wall -O2)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

include_directories("../../src")
include_directories("../../src/ErrorCorrectionCode")
include_directories("../../src/lib")
include_directories("../../src/channel")
include_directories(".")
include_directories("common")

set(common_src
        ../../src/lib/libMath.cpp
        ../../src/lib/libParser.cpp
        ../../src/channel/channelBase.cpp
        ../../src/channel/AWGN.cpp
        ../../src/ErrorCorrectionCode/ErrorCorrectionCodeBase.cpp
        ../../src/ErrorCorrectionCode/AED_attemp_20260102.cpp
        ../../src/ErrorCorrectionCode/AED_relation_check.cpp
    )

add_executable(mclass_dataset
        past_experiments/mclass_dataset.cpp
        ${common_src}
    )

add_executable(mclass_eval
        tools/mclass_eval.cpp
    )

add_executable(mclass_selected_llr
        past_experiments/mclass_selected_llr.cpp
        ${common_src}
    )

add_executable(mclass_llr_proxy
        past_experiments/mclass_llr_proxy.cpp
        ${common_src}
    )

add_executable(mclass_branch_orders
        past_experiments/mclass_branch_orders.cpp
        ${common_src}
    )

add_executable(mclass_perturb_audit
        past_experiments/mclass_perturb_audit.cpp
        ${common_src}
    )

add_executable(mclass_trace_dataset
        common/mclass_trace_dataset.cpp
        ${common_src}
    )

add_executable(mclass_code_structure
        past_experiments/mclass_code_structure.cpp
        ${common_src}
    )

add_executable(mclass_gmatrix_dump
        past_experiments/mclass_gmatrix_dump.cpp
        ${common_src}
    )

add_executable(mclass_bler_sim
        common/mclass_bler_sim.cpp
        ${common_src}
    )
