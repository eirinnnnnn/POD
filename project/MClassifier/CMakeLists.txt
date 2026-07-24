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
        mclass_dataset.cpp
        ${common_src}
    )

add_executable(mclass_eval
        mclass_eval.cpp
    )

add_executable(mclass_selected_llr
        mclass_selected_llr.cpp
        ${common_src}
    )

add_executable(mclass_llr_proxy
        mclass_llr_proxy.cpp
        ${common_src}
    )

add_executable(mclass_branch_orders
        mclass_branch_orders.cpp
        ${common_src}
    )

add_executable(mclass_perturb_audit
        mclass_perturb_audit.cpp
        ${common_src}
    )

add_executable(mclass_trace_dataset
        mclass_trace_dataset.cpp
        ${common_src}
    )

add_executable(mclass_code_structure
        mclass_code_structure.cpp
        ${common_src}
    )

add_executable(mclass_gmatrix_dump
        mclass_gmatrix_dump.cpp
        ${common_src}
    )
