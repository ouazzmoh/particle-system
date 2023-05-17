# Créer un executable à partir de fichiers sources
add_executable(test_1  test_1.cxx)
add_executable(listVSset  listVSset.cxx)
add_executable(gravSystem  gravSystem.cxx)
add_executable(test_universe test_universe.cxx)
add_executable(test_potentiel test_potentiel.cxx)
add_executable(test_universe_algo test_universe_algo.cxx)
add_executable(test_universe_tp3 test_universe_tp3.cxx)
add_executable(test_application test_application.cxx)
add_executable(test_small test_small.cxx)
add_executable(test_one test_one.cxx)

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

target_link_libraries(test_small particules)
target_link_libraries(test_small vecteur)
target_link_libraries(test_small universe)
target_link_libraries(test_small cell)

target_link_libraries(test_one particules)
target_link_libraries(test_one vecteur)
target_link_libraries(test_one universe)
target_link_libraries(test_one cell)