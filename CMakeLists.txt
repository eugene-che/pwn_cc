cmake_minimum_required(VERSION 2.8)
project(pwn_cc C CXX)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++14 -Wall")

set(PWN_CC_SRC
  pwn_cc.hpp
  pwn_cc.cpp
)

add_library(pwn_cc ${PWN_CC_SRC})

find_package(Boost COMPONENTS unit_test_framework)
if(${Boost_UNIT_TEST_FRAMEWORK_FOUND})
  enable_testing()
  include_directories(${Boost_INCLUDE_DIRS})
  add_definitions(-DBOOST_TEST_DYN_LINK)

  add_executable(pwn_cc_test pwn_cc_test.cpp)
  target_link_libraries(pwn_cc_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} pwn_cc)
  add_test(pwn_cc_test pwn_cc_test)
endif()
