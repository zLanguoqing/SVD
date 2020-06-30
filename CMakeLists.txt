cmake_minimum_required(VERSION 2.8)
project(SVD)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})

# Check C++11 or C++0x support
add_definitions(-DCOMPILEDWITHC11)
if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -msse3 -std=c++11 -pthread -g -march=native -Wno-deprecated-declarations")
else( CMAKE_BUILD_TYPE MATCHES "Release" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -msse3 -std=c++11 -pthread -O3 -march=native -Wno-deprecated-declarations")
endif( CMAKE_BUILD_TYPE MATCHES "Debug" )

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

find_package(OpenCV 2.4.9 REQUIRED)

include_directories(
${PROJECT_SOURCE_DIR}
${PROJECT_SOURCE_DIR}/include
)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

add_library(${PROJECT_NAME} SHARED
src/SVD.cpp
)

target_link_libraries(${PROJECT_NAME}
${OpenCV_LIBS}
)

add_executable(main src/main.cpp)
target_link_libraries(main ${PROJECT_NAME})
