cmake_minimum_required(VERSION 3.3)
project(rpkm_BasicContig)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    fastareader.c++
    fastareader.h
    helper.c++
    helper.h
    matchoutputparser.c++
    matchoutputparser.h
    rpkm.c++
    rpkm.h
    types.h
    utilities.c++
    utilities.h)

add_executable(rpkm_BasicContig ${SOURCE_FILES})