cmake_minimum_required(VERSION 3.21)
project(DecisionTree_TSMFVC)

set(CMAKE_CXX_STANDARD 14)

add_executable(DecisionTree_TSMFVC main.cpp TSMFVC.h PKE.h Eval_DT.h tool.h)

target_link_libraries(DecisionTree_TSMFVC libgmp.so libflint.so libpbc.so m.so librelic_s.a /usr/local/lib/libntl.so)