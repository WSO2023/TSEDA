cmake_minimum_required(VERSION 3.14)
project(Frame)

include_directories(${PROJECT_SOURCE_DIR}/include)
set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS "-fopenmp")
add_subdirectory(${PROJECT_SOURCE_DIR}/lib)

add_executable(Frame main.cpp lib/pugixml.cpp include/pugixml.hpp  include/pugiconfig.hpp lib/classAndVarDefine.cpp include/classAndVarDefine.h include/common.h lib/config.cpp include/config.h lib/tools.cpp include/tools.h lib/GenOperator.cpp include/GenOperator.h lib/GenerateAchrom.cpp include/GenerateAchrom.h   NGA.cpp NGA.h HGA.cpp HGA.h   TSEDA.cpp TSEDA.h LWSGA.cpp LWSGA.h HMEC.cpp HMEC.h CGA.cpp CGA.h HEFT.cpp HEFT.h HPSO.cpp HPSO.h ADBRKGA.cpp ADBRKGA.h)
