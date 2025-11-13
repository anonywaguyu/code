cmake_minimum_required(VERSION 3.21)
project(TSNIVC)

set(CMAKE_CXX_STANDARD 14)

add_executable(TSNIVC main.cpp function.h TSNIVC.h GenData.h tool.h PKE.h)

target_link_libraries(TSNIVC libgmp.so libflint.so libpbc.so m.so librelic_s.a /usr/local/lib/libntl.so)