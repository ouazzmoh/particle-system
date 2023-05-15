# Créer une bibliothèque à partir d'un ensemble de fichiers de définitions (.cxx) sans main.
# Ces fichiers ne contiennent pas de fonction main
# Remplacer ... par les fichiers nécessaires 
add_library(particules particule.cxx)
add_library(vecteur vecteur.cxx)
add_library(universe universe.cxx)
add_library(cell cell.cxx)
