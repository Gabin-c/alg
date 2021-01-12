# Projet ALG 2020
***
## Mappeur sans gap, détection de SNPs
Gabin Coudray - David Gallien  
Master 2 BIS
***

Ce programme a pour but d'élaborer un mappeur sans gap pour la détection de SNPs. Il se compose de deux fichier.py 
distincts :
- Le premier, **index.py**, permet l'indexation du génome de référence
- Le second, **map.py**, permet le mapping des reads sur le génome de référence

Le fichier **tools_karkkainen_sanders.py** est utilisé pour l'indexation.

Les dossiers *coli*, *covid* et *mallMappingTest* contiennent les jeux des données que nous avons testé.

Le dossier *rapports* contient 2 fichiers markdown :

- **rapport_developpeur.md** qui explique la structure du programme et le fonctionnement des fonctions majeures
- **rapport_utilisateur.md** qui présente les résultats obtenus selon les paramètres utilisés et pour les différents jeux de données.
