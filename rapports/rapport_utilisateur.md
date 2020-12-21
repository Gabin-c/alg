# RAPPORT UTILISATEUR
Gabin Coudray - David Gallien  
Master 2 BIS
***

## PARTIE I : *PRINCIPES DE BASE DE L'ASSEMBLEUR*
### Principe de base
Ce mappeur permet de mapper des reads sur un génome de référence sans prendre en compte les gaps. Il a 
pour but de détecter des SNPs. Le programme repose sur le principe de la transformée de Burrows-Wheeler qui permet de 
réorganiser la séquence du génome de référence afin de la rendre plus simple à requêter.

Le programme se décompose en deux fichiers pythons :  
- Le premier fichier *index.py* permet de créer un FM index à partir de la séquence du génome de référence et de le 
stocker dans un fichier. Ce FM index se présente sous la forme d'un tuple contenant 4 éléments :
  - La transformée de Burrows-Wheeler
  - Le Suffix Array (sa)
  - Le nombre n de chaque caractère
  - Le rang r de chaque caractère 
  
```shell
FMI = ('BWT', [SA], {n}, [r])    
```
- Le deuxième fichier *map.py* a pour but de mapper les reads sur la séquence du génome de référence en utilisant le 
FM index créé précédemment. Cela dans le but de détecter des variants (SNPs) et de les stocker dans un fichier VCF sous
la forme suivante :
```shell script
    #REF: fichier du génome de référence
    #READS: fichier des reads
    #K: longueur du kmer
    #MAX_SUBST: maximum de substitution
    #MIN_ABUNDANCE: minimum d'abondance des SNPs
    POS REF ALT ABUNDANCE
```


### Utilisation de l'assembleur
Le programme fonctionne en ligne de commande grâce à deux commandes spécifiques. La première pour la création du FM indexe :

```shell script
python index.py --ref [genome_file.fa] --out [dumped_index.dp]
```
Elle prend donc en entrée (*'--ref'*) le fichier fasta de la séquence de référence et en sortie (*'--out'*) le fichier 
binaire de sortie contenant le FM index au format *dp*.


La deuxième commande permet le mapping et la création du fichier VCF contenant les différents SNPs.
```shell script
python map.py --ref[genome_file.fa] --index[dumped_index.dp] --reads[reads.fa] -k[k_value] --max_hamming[h_value] --min_abundance[m_value] --out snps.vcf
```
Cette commande prendre plusieurs arguments en compte :
- *--ref* : le fichier fasta de la séquence de référence
- *--index* : le fichier *dp* contenant le FM index
- *--reads* : le fichier fasta contenant les reads à mapper
- *-k* : la longueur du kmer à considérer pour l'ancrage des reads sur le génome
- *--max_hamming* : le nombre maximum de substitutions pour le mapping
- *--min_abundance* : le minimum d'abondance des SNPs à considérer pour le fichier VCF
- *--out* : le fichier VCF de sortie

Pour obtenir l'aide des deux programmes et ainsi savoir la commande à taper ainsi que les paramètres par défaut on peut 
faire comme ci-dessous : 
````shell
python index.py -h 
python map.py -h
````
 

## PARTIE II : *RÉSULTATS OBTENUS*
### Premiers tests
Pour commencer, nous avons testé notre mappeur sur un petit jeu de données afin de vérifier s'il fonctionnait correctement.
C'est un jeu de données avec une séquence de référence de 1000 nucléotides et 11 reads de 100 nucléotides. Cela permet de 
créer le FM index et de mapper les reads sur la séquence de référence en moins de 1 seconde. 
Nous avions des modèles de résultats avec différents paramètres pour le mapping (k et max_hamming) et différentes séquences
de référence. Il s'avère qu'avec notre programme nous retrouvons exactement les mêmes résultats pour chacun des 
différents résultats modèles.

Suite à cela, nous pouvons donc passer à un jeu de données beaucoup plus important qui est celui de *Escherichia coli*. 


### Résultats pour *Escherichia coli*
Nous avons donc un jeu de données correspondant à *Escherichia coli*. Il contient un fichier fasta des 150 000 premiers 
nucléotides du génome de *Escherichia coli* et un fichier de 30 000 reads. Pour mapper ces 30 000 reads, notre programme
a besoin de 120 secondes (2 minutes) quels que soient les paramètres pris en compte.  
Nous avons fait des tests avec différents paramètres afin de trouver ceux donnant les meilleurs résultats. Pour cela, nous
avons fait étape par étape en partant des paramètres suivants :
```shell
-k 20
--max_hamming 10
--min_abundance 1
```
Avec ces paramètres nous avons obtenu des résultats de faible qualité avec une précision de 3.6%. En effet, on retrouvait
26 902 SNPs faux positifs. Pour améliorer les résultats nous avons dans un premier temps augmenté l'abondance minimum petit
à petit. A chaque augmentation, la précision augmentait et le nombre de faux positifs diminuait. Nous avons donc trouvé 
le minimum d'abondance idéal qui est de 7. Il s'agit de l'abondance qui permet d'avoir le taux de précision le plus élevé
et le taux de vrai positif le plus élevé. Voici le résultat ci-dessous avec une précision de 99.7% et un recall de 99.1%.
Si on passe une abondance supérieure à 7, le nombre de faux négatifs augmente (recall diminue) et la précision ne diminue que
très légèrement.
```
Pour k = 20, max_hamming = 10, min_abundance = 7 :
Nb Variants 1004 in truth
Nb TP 995 in truth & in pred
Nb FP 3 not in truth & in pred
Nb FN 9 in truth & not in pred
Precision 99.7 %
Recall 99.1 %
```
Avec ces paramètres, nous avons pu voir que notre programme utilise environ 930 000kb de mémoire.
Ensuite, nous avons donc gardé une abondance de 7 tout en faisant varier le maximum de substitution. En augmentant ce maximum,
nous avons une précision qui diminue et en le diminuant on retrouve une augmentation du nombre de faux positifs. Nous avons 
donc gardé un maximum de substitution de 10 car l'utilisation mémoire n'était pas impactée par ce paramètre.  
Pour finir, nous avons fait varier la longueur du kmer. L'augmentation de la longueur du kmer, le temps de calcul diminue 
et l'utilisation mémoire diminue. En effet si on augmente à 50, on a un temps de calcul de 88 secondes et une utilisation mémoire 
de 580 000kb. Cependant on perd en qualité car on retrouve 136 faux négatifs.
En augmentant un peu moins la longueur du kmer (à 30), le temps n'est pas impacté mais l'utilisation mémoire diminue à 800 000kb. 
De plus, la précision reste à 99.7%. Cependant le nombre de faux négatif passe de 9 à 12.  
Il faut donc faire un choix entre diminution de l'utilisation mémoire ou diminution du nombre de faux négatifs. Nous avons fait 
le second choix en gardant une longueur de kmer de 20.
A l'inverse, en diminuant la longueur du kmer, on obtient des résultats identiques avec une précision identique. Seulement,
l'utilisation de la mémoire augmente fortement (plus de 1GB) et le temps de calcul augmente légèrement.

Pour résumer, nous avons donc gardé en paramètres par défaut :
 - Une longueur de kmer de 20
 - Un maximum de substitution de 10
 - Un minimum d'abondance de 7  

Il est possible d'augmenter la taille des kmer jusqu'à 30 sans trop impacter la qualité des résultats pour avoir une utilisation
mémoire plus faible. 

En voyant ces résultats nous pouvons dire que notre approche permet d'obtenir des résultats de qualité avec les paramètre 
par défaut. De plus, elle permet un mapping relativement rapide. En revanche, l'utilisation mémoire semble élevée et nous 
n'avons pas su la diminuer.



## PARTIE 3 : LES DONNÉES COVID
