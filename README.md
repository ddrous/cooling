# Projet "cooling"
Visualisation du dispositif de refroidissement d’un microprocesseur en C++, et visualisation à l'aide de Gnuplot et Paraview.

![quelques images](/rapport/overview.png)

## Travail à effectuer et rapport
- `projet.pdf`  
- `rapport/rapport_cooling.pdf`

## Compilation
 À partir du répertoire `src`: **`g++ cooling.cpp stationnaire.cpp instationnaire.cpp -o cooling`**

## Execution
Les fichier de configuration sont regroupés dans le répertoire `config_files`: **`./src/cooling config_files/simu_1.cfg`**

## Visualisation
- La visulation avec Gnuplot se fait automatiquement (en supposant que Gnuplot est installé)
- Pour visualiser l'évolution de la température avec Paraview, il faut utiliser les fichiers .vtk du répertoire `paraview_data`
