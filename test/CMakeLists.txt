# Créer un executable à partir de fichiers sources
add_executable(test_1  test_1.cxx)
add_executable(listVSset  listVSset.cxx)
add_executable(gravSystem  gravSystem.cxx)
add_executable(test_universe test_universe.cxx)
add_executable(test_potentiel test_potentiel.cxx)
add_executable(test_universe_algo test_universe_algo.cxx)
add_executable(test_universe_tp3 test_universe_tp3.cxx)
add_executable(test_application test_application.cxx)
add_executable(test_boundaries test_boundaries.cxx)
add_executable(test_interaction test_interaction.cxx)

# Fait le lien avec la bibliothèque créée précedemment
# Les entêtes sont automatiquement détectés.
target_link_libraries(test_1 particules)
target_link_libraries(test_1 vecteur)
target_link_libraries(test_1 universe)
target_link_libraries(test_1 cell)

target_link_libraries(listVSset particules)
target_link_libraries(listVSset vecteur)
target_link_libraries(listVSset universe)
target_link_libraries(listVSset cell)

target_link_libraries(gravSystem particules)
target_link_libraries(gravSystem vecteur)
target_link_libraries(gravSystem universe)
target_link_libraries(gravSystem cell)

target_link_libraries(test_universe particules)
target_link_libraries(test_universe vecteur)
target_link_libraries(test_universe universe)
target_link_libraries(test_universe cell)

target_link_libraries(test_potentiel particules)
target_link_libraries(test_potentiel vecteur)
target_link_libraries(test_potentiel universe)
target_link_libraries(test_potentiel cell)


target_link_libraries(test_universe_algo particules)
target_link_libraries(test_universe_algo vecteur)
target_link_libraries(test_universe_algo universe)
target_link_libraries(test_universe_algo cell)

target_link_libraries(test_universe_tp3 particules)
target_link_libraries(test_universe_tp3 vecteur)
target_link_libraries(test_universe_tp3 universe)
target_link_libraries(test_universe_tp3 cell)

target_link_libraries(test_application particules)
target_link_libraries(test_application vecteur)
target_link_libraries(test_application universe)
target_link_libraries(test_application cell)

target_link_libraries(test_boundaries particules)
target_link_libraries(test_boundaries vecteur)
target_link_libraries(test_boundaries universe)
target_link_libraries(test_boundaries cell)

target_link_libraries(test_interaction particules)
target_link_libraries(test_interaction vecteur)
target_link_libraries(test_interaction universe)
target_link_libraries(test_interaction cell)
