cmake_minimum_required(VERSION 3.10)
project(newtest)

set(CMAKE_CXX_STANDARD 11)

set(EASY_BMP EasyBMP.h EasyBMP_BMP.h EasyBMP.cpp EasyBMP_DataStructures.h EasyBMP_VariousBMPutilities.h)
set(SOURCE_FILES main.cpp)

set(CMAKE_CXX_FLAGS "-pedantic -Wall")

add_executable(newtest ${SOURCE_FILES} ${EASY_BMP})