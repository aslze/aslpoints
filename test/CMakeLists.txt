
add_executable(tests1 tests.cpp)
target_link_libraries(tests1 asls aslpoints)
add_test(NAME fitting COMMAND tests1)

add_executable(tests2 tests2.cpp)
target_link_libraries(tests2 asls aslpoints)
add_test(NAME fitting_robust COMMAND tests2)
