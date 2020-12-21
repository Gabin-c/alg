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
Elle prendre donc en entrée (*'--ref'*) le fichier fasta de la séquence de référence et en sortie (*'--out'*) le fichier 
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
temps de calculs
mémoire utilisée
qualité des résultats
EN FONCTION DES VALEURS DES PARAMETRES

Cette partie inclura une discussion et des conclusions quant aux qualités et défauts de
l'approche proposée, et guidera l'utilisateur sur le choix des paramètres.

## PARTIE 3 : LES DONNÉES COVID
