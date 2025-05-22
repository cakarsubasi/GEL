cmake_minimum_required(VERSION 3.30)

get_filename_component(GEL_INCLUDE_DIR ../.. ABSOLUTE)

add_executable(RsR_test RsR_test.cpp)
add_executable(OBJParser_test OBJParser_test.cpp)

set_target_properties(RsR_test OBJParser_test PROPERTIES
    CXX_STANDARD 20
    INCLUDE_DIRECTORIES "${GEL_INCLUDE_DIR};${DOCTEST_INCLUDE_DIR}"
    LINK_LIBRARIES GEL
)

add_test(NAME GEL.Test.RsR COMMAND RsR_test)
add_test(NAME GEL.Test.OBJParser COMMAND OBJParser_test)

#add_test(NAME GEL.Test.RsRPython COMMAND ${PYTEST} "${CMAKE_CURRENT_SOURCE_DIR}/rsr.py")
#add_test(NAME GEL.Test.RsRPython COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/rsr.py")
