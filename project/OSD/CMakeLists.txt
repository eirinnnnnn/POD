cmake_minimum_required(VERSION 2.6)

#set(CMAKE_BUILD_TYPE Release)
add_definitions(-Wall -O2)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

include_directories("../../src")
include_directories("../../src/ErrorCorrectionCode")
include_directories("../../src/lib")
include_directories("../../src/channel")

set(src main.cpp
        ../../src/ErrorCorrectionCode/ErrorCorrectionCodeBase.cpp
        ../../src/ErrorCorrectionCode/OSD.cpp
        ../../src/lib/libMath.cpp
        ../../src/lib/libParser.cpp
        ../../src/channel/channelBase.cpp
        ../../src/channel/AWGN.cpp
    )

set(header 	../../src/lib/libDebug.h
            ../../src/ErrorCorrectionCode/ErrorCorrectionCodeBase.h
			../../src/ErrorCorrectionCode/OSD.h
            ../../src/lib/libMath.h
            ../../src/lib/libParser.h
            ../../src/channel/channelBase.h
            ../../src/channel/AWGN.h
    )

add_executable(OSD_performance ${src} ${header})