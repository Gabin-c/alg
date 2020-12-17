<div align="center"><h1> RAPPORT DÉVELOPPEUR </h1> 
Gabin Coudray - David Gallien
<br>Master 2 BIS</div>
<hr>


<div align="justify">
Ce programme a pour but de mettre en place un mappeur sans gap pour la détection de 
SNPs.
Il se distingue en deux fichier.py distincts :
<li>Le premier permet l'indexation du génome de référence</li>
<li>Le second permet le mapping des reads sur le génome de référence</li>
<hr>


<h1>Premier programme : index.py</h1>

Le fichier d'indexation nommé <i>index.py</i> a pour but de créer un FM index 
composé de :
<li> La transformée de Burrows-Wheeler (BWT) du génome de référence </li>
<li> Le suffix array (SA)</li>
<li> Le nombre de chaque caractère (n) </li>
<li> Le rang de chaque caractère (r)</li>
Les 3 derniers éléments de cet FM index permettent de parcourir et requêter
la transformée de Burrows Wheeler.


<h2>Fonctionnement du programme </h2>
La fonction <i>get_fmi(ref_fasta, output_file)</i> prend en entrée le génome de 
référence au format fasta ainsi que le nom choisi pour le fichier contenant le FM 
index créé. En sortie, la fonction retourne un fichier au format binaire contenant 
le FM index.

Le programme se déroule en deux étapes : dans un premier temps la transformée de 
Burrows-Wheeler va être déterminée à l'aide du génome de référence puis le rang et le
nombre de chaque caractère vont être identifiés à l'aide de la BWT.
A la fin de la première étape nous aurons le Suffix Array et la BWT et à la fin de 
la deuxième étape nous aurons les deux derniers éléments de l'index r et n.

<h3> Première étape </h3>
<h4> get_bwt(fasta) </h4>
La fonction <i>get_bwt(fasta: str)</i> prend en entrée le génome de référence au 
format fasta et permet d'avoir en sortie la BWT ainsi que le suffix array.
<br>Elle fonctionne de la manière suivante : 

Dans un premier temps la fonction <i>get_seq(fasta: str)</i> traite le fichier du 
génome de réference. Elle évite la première ligne commençant par ">" puis stock
la séquence du génome de référence (s). Enfin, elle créé le suffix array à l'aide 
la fonction <i>kark_sort()</i> qui est importée de <i>tools_karkkainen_sander.py</i>.
Ensuite, avec la séquence et le suffix array, on peut obtenir la BWT.
Enfin, la fonction retourne le suffix array et la transformée de 
Burrows-Wheeler, qui seront stockés dans le FM index.

<h3> Deuxième et dernière étape</h3>
<h4> get_r_n(bwt) </h4>
  
La fonction <i>get_r_n(bwt)</i> va calculer à partir de la transformée 
Burrows-Wheeler le rang r et le nombre n pour chaque caractère du génome de référence.

Dans un premier temps la fonction initialise un dictionnaire n avec les 4 nucléotides
et une liste r pour stocker les rangs.
Suite à cela, la BWT va être parcourue et on va implémenter au fur et à mesure 
notre dictionnaire n et notre list r. On obtiendra donc en sortie les deux derniers
élément de notre FM index : n et r. 

Ce programme doit être exécuté avec la commande suivante :

    python index.py --ref [genome_file.fa] --out [dumped_index.dp]
<hr>

![alt text](index_diag.png "Diagramme représentant le programme index.py")

<hr>
<h1>Deuxième programme : map.py </h1> 
L'objectif de ce progamme est l'alignement des reads sur la séquence de référence 
tout en choisissant les paramètres de l'alignement : nombre de substitutions maximum, 
longueur du k-mer et minimum d'abondance du SNP.

</div>
<h2> mapping(ref, index, reads, k, max_hamming, min_abundance, out_file) </h2>
<div align="justify">

La fonction <i>mapping()</i> prend en paramètres le génome de référence (ref), un
FM index (index) (voir index.py), la longueur du k-mer (k), le nombre maximum de 
substitutions (max_hamming), l'abondance minimum pour un SNP (min_abundance) et 
le nom du fichier vcf de sortie (out_file).

Cette fonction va retourner un tableau de SNP au format vcf obtenu grâce à 
l'alignement des reads sur le génome.

La fonction se décompose en 3 étapes :

<li>Création d'un dictionnaire contenant les positions d'alignement pour chaque 
k-mer de chaque read</li>
<li>Création d'un dictionnaire contenant la meilleur position d'alignement pour 
chaque read</li>
<li>Remplissage de la table des SNPs</li>

<h3> Etape 1 : Dictionnaire <i>{read: {kmer: positions}}</i> </h3>
Dans un premier temps le fichier contenant le FM index va être ouvert afin de rendre 
accessible les différents éléments du FM index.


<h4> get_kmer_position(k, reads, index) </h4>
Suite à cela les positions d'alignements des k-mers pour chaque read vont etre 
stockés dans un dictionnaire par la fonction <i>get_kmer_position(k, reads)</i>
avec 'k' la longueur des k-mers et 'reads' le fichier contenant les reads. 
Après l'utilisation de cette fonction on obtient donc un dictionnaire au format 
{Read: {K-mer: position}} pour chaque sens du brin.

Cette fonction permet d'obtenir la ou les positions de chaque kmer sur le brin sens 
et son reverse complement.

<h3> Etape 2 : Dictionnaire <i>{read:position}</i> </h3>
Une fois qu'on a toutes les positions possibles pour un kmer, la fonction passe à 
l'étape de la localisation de la meilleure position d'alignement de chaque read sur
le génome de référence pour cela <i>mapping()</i> utilise la trame suivante :
<br> Chaque position associée au read du dictionnaire va être lue pour les brins sens
et antisens. Ensuite un score d'alignement est calculé sans gap puis stocké. Enfin,
la position avec le score le plus élevé est gardée en prenant en compte le nombre de
substitutions maximum. On obtient un dictionnaire final qui contient la meilleure 
position d'alignement pour chaque read.

<h3> Etape 3 : Table SNP </h3>
Une fois le dictionnaire des meilleures positions d'alignement obtenu, on passe à 
l'étape de remplissage de la table de SNPs.
Pour cela, on a deux fonctions. La première 
<i>fill_vcf(mat, dict_final, sequence_initiale)</i>
pour remplir ce tableau et la seconde <i>order_vcf(tab_vcf)</i>
pour ordonner les SNPs selon l'ordre croissant des positions.
<br> Ce tableau est stocké dans un fichier au format vcf de la forme suivante :

    #REF: fichier du génome de référence
    #READS: fichier des reads
    #K: longueur du kmer
    #MAX_SUBST: maximum de substitution
    #MIN_ABUNDANCE: minimum d'abondance des SNPs
    POS REF ALT ABUNDANCE

Ce programme doit être exécuté avec la commande suivante :

    python map.py --ref[genome_file.fa] --index[dumped_index.dp] --reads[reads.fa] -k[k_value] --max_hamming[h_value] --min_abundance[m_value] --out snps.vcf

</div>
<hr>

![alt text](map_diag.png "Diagramme représentant le programme map.py")

<hr>



  





