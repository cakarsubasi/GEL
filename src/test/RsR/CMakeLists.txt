cmake_minimum_required(VERSION 3.30)

add_executable(RsR_test RsR_test.cpp)
target_include_directories(RsR_test PRIVATE ../..)
set_target_properties(RsR_test PROPERTIES
        CXX_STANDARD 20)

target_link_libraries(RsR_test GEL)
add_test(NAME GEL.Test.RsR COMMAND RsR_test)

#add_test(NAME GEL.Test.RsRPython COMMAND ${PYTEST} "${CMAKE_CURRENT_SOURCE_DIR}/rsr.py")
#add_test(NAME GEL.Test.RsRPython COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/rsr.py")
