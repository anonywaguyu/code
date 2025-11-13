cmake_minimum_required(VERSION 3.21)
project(TSMFVC)

set(CMAKE_CXX_STANDARD 14)

add_executable(TSMFVC main.cpp function.h TSMFVC.h GenData.h tool.h PKE.h)

target_link_libraries(TSMFVC libgmp.so libflint.so libpbc.so m.so librelic_s.a /usr/local/lib/libntl.so)