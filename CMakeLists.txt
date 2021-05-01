cmake_minimum_required (VERSION 3.8)
project(INS)

set(CMAKE_CXX_STANDARD 14)
add_library(INS INS.cpp)
add_library(get_state_transition get_state_transition.cpp)
add_executable(main main.cpp)
target_link_libraries(main
        INS
        get_state_transition
)