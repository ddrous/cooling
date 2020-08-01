# Projet "cooling"
Visualisation du dispositif de refroidissement d’un microprocesseur à l'aide de gnuplot et Paraview.

## Travail à effectuer et rapport
``` projet.pdf ```  
``` rapport/rapport_cooling.pdf ```

## Compilation
``` g++ cooling.cpp stationnaire.cpp instationnaire.cpp -o cooling ``` 

## Execution
``` ./cooling config_files/simu_1.cfg ```

## Visualisation
- La visulation avec `gnuplot` se fait automatiquement (en supposant que `gnuplot` est installé)
- Pour visualiser l'évolution de la température avec Paraview, il faut utiliser les fichiers .vtk du répertoire `paraview_data`
