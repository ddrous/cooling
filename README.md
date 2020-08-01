# Projet "cooling"
Simulation numérique du dispositif de refroidissement d’un microprocesseur à l'aide de Gnuplot et Paraview.

## Travail à effectuer et rapport
- `projet.pdf`  
- `rapport/rapport_cooling.pdf`

## Compilation
`g++ cooling.cpp stationnaire.cpp instationnaire.cpp -o cooling` (à partir du répertoire `src`)

## Execution
`./src/cooling config_files/simu_1.cfg` (les fichier de configuration sont regroupés dans le répertoire `config_files`)

## Visualisation
- La visulation avec Gnuplot se fait automatiquement (en supposant que Gnuplot est installé)
- Pour visualiser l'évolution de la température avec Paraview, il faut utiliser les fichiers .vtk du répertoire `paraview_data`
