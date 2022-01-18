cmake_minimum_required(VERSION 3.21)
project(homomorphic_encryption_demo)

set(CMAKE_CXX_STANDARD 14)

add_executable(homomorphic_encryption_demo
        RLWE.cpp)
