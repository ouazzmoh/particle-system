### TP C++ : Rendu Intermediare
#### Ouazzani Chahdi Mohammed
#### Cherkaoui Ayoub

Pour compiler et executer : 
    cd build 
    cmake ..
    make 
    ./test/<nom_du_test>

Description des tests:
    gravSystem : Construire un univers classique, et simuler le mouvement 
    listVSset : Comparaison de performance des list et des set pour l'insertion (Remarque: Nous changerons vers une structure Vector au lieu de List pour stocker les particules)
    test_1 : Test des operations sur les vecteurs
    test_application :  Creer les 2 bloques de particules et lancer l'interaction avec la force de potentiel en considerant le rCut
    test_potential : Tracer le potentiel en fonction du rayon de coupure
    test_universe_algo :  Test de stromer verlet sur les forces potentiel
    test_universe_tp3 : Optimisation du calcul de force d'interaction (Questions TP3)
    test_universe : Teest du decoupage de l'espace



La version avec les templates pour géneraliser à des dimensions multiples n'est toujours pas stable dû à des problèmes de liens, alors nous avons soumis une version sans les templates.
