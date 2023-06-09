include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(vecteur_test vecteur_test.cxx)
add_test(NAME vecteur_test COMMAND vecteur_test)
target_link_libraries(vecteur_test gtest_main)
add_executable(grid_test_one_particle grid_test_one_particle.cxx)
add_test(NAME grid_test_one_particle COMMAND grid_test_one_particle)
target_link_libraries(grid_test_one_particle gtest_main)
add_executable(particule_test particule_test.cxx)
add_test(NAME particule_test COMMAND particule_test)
target_link_libraries(particule_test gtest_main)
add_executable(grid_test_two_particle grid_test_two_particle.cxx)
add_test(NAME grid_test_two_particle COMMAND grid_test_two_particle)
target_link_libraries(grid_test_two_particle gtest_main)
add_executable(test_universe_boundaries test_universe_boundaries.cxx)
add_test(NAME test_universe_boundaries COMMAND test_universe_boundaries)
target_link_libraries(test_universe_boundaries gtest_main)
add_executable(gravitation_test gravitation_test.cxx)
add_test(NAME gravitation_test COMMAND gravitation_test)
target_link_libraries(gravitation_test gtest_main)



include(GoogleTest)
gtest_discover_tests(vecteur_test)
gtest_discover_tests(grid_test_one_particle)
gtest_discover_tests(grid_test_two_particle)
gtest_discover_tests(particule_test)
gtest_discover_tests(test_universe_boundaries)
gtest_discover_tests(gravitation_test)






target_link_libraries(vecteur_test particules)
target_link_libraries(vecteur_test vecteur)
target_link_libraries(vecteur_test universe)
target_link_libraries(vecteur_test cell)



target_link_libraries(grid_test_one_particle particules)
target_link_libraries(grid_test_one_particle vecteur)
target_link_libraries(grid_test_one_particle universe)
target_link_libraries(grid_test_one_particle cell)


target_link_libraries(grid_test_two_particle particules)
target_link_libraries(grid_test_two_particle vecteur)
target_link_libraries(grid_test_two_particle universe)
target_link_libraries(grid_test_two_particle cell)

target_link_libraries(particule_test particules)
target_link_libraries(particule_test vecteur)
target_link_libraries(particule_test universe)
target_link_libraries(particule_test cell)


target_link_libraries(test_universe_boundaries particules)
target_link_libraries(test_universe_boundaries vecteur)
target_link_libraries(test_universe_boundaries universe)
target_link_libraries(test_universe_boundaries cell)


target_link_libraries(gravitation_test particules)
target_link_libraries(gravitation_test vecteur)
target_link_libraries(gravitation_test universe)
target_link_libraries(gravitation_test cell)