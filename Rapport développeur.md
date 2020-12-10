Le programme se distingue en deux fichier .py distinct, le premier traite de l'indexation du génome de référence et le second du mapping des reads sur le génome de référence


<h1>index.py</h1>

Le fichier d'indexation nommé index.py a pour but de crée un index composé de la transformé de Burrows wheeler du génome de référence, ainsi que les informations permettant l'utilisation de cette transformée de burrows wheeler.

<h2>get_fmi(ref_fast, output_file) </h2>

La fonction get_fmi(ref_fasta, output_file ) permet de crée cet index, elle prend en entrée le génome de référence au format fasta, ainsi que le nom choisis pour nommé l'index crée, en sortie la fonction retourne un index au format pickle contenant la transformée de burrows wheeler (bwt), le suffix array (sa) correspondant ainsi que le rang (r) et le nombre (n) de chaque caractère.

La fonction fonctionne en deux étapes. Dans un premier temps la transformée de burrows wheeler va être déterminer à l'aide du génome de référence puis le rang et le nombre de chaque caractère va être identifier à l'aide de la bwt.
A la fin de la première étape nous aurons sa et bwt et à la fin de la deuxième étape nous aurons les deux derniers éléments de l'index r et n.

<h3> Première étape </h3>
<h4> get_bwt(fasta) </h4>
La fonction get_bwt(fasta: str) prend en entré le génome de référence au format fasta et permet d'avoir en sortie la bwt ainsi que le sa.

La fonction fonctionne de la manière suivante : 

Dans un premier temps la fonction get_seq( fasta: str ) traite le fichier du génome de réference, elle saute la première ligne consititué de la ligne informative 
commençant par ">" puis stocke la séquence du génome de référence (s) et crée le suffix array (sa) à l'aide la fonction kark_sort().

Puis avec sa et s, la transformée de burrows wheeler est formé.

La fonction retourne ensuite le suffix array sa et la transformée de burrows wheeler bwt, qui seront stocké dans le FM index.

<h3> Deuxième et dernière étape</h3>
<h4> get_r_n(bwt) </h4>
  
La fonction get_r_n(bwt) va calculer à partir de la transformée burrows wheeler le rang r et le nombre n pour chaque caractère du génome de référence.

Dans un premier temps la fonction initialise un dictionnaire n avec les 4 nucléotides et une liste r pour stocké les rangs.
Suite à ça, la bwt va être parcouru et va implémenter au fur et à mesure notre dictionnaire n et notre list r. On obtiendra en sortie un dictionnaire n et une liste r, nos deux derniers éléments de notre FM index.

<h1> map.py </h1> 

L'objectif de se progamme est l'alignement de reads sur la séquence de référence, tous en choisissant les paramètres de l'alignement : nombre de substitutions max, longueur du k-mer.

<h2> mapping(ref, index, reads, k, max_hamming, min_abundance, out_file) </h2>

La fonction mapping() prend en paramètres le génome de référence (ref), un FM index (index) (voir index.py), la longueur du k-mer (k), le nombre maximum de substitution (max_hamming), l'abondance minimum pour un SNP (min_abundance) et le nom du fichier vcf de sortie (out_file).

Cette fonction va retourner un tableau de SNP au format vcf obtenu grâce à l'alignement des reads sur le génome.

La fonction fonctionne de la façon suivante:

Dans un premier temps le fichier index va être ouvert afin de rendre accessible les différents éléments de l'index.


Suite à cela les positions d'alignements des k-mers pour chaque read vont etre stockés dans un dictionnaire par la fonction get_kmer_position(k , reads)
avec k la longueur des k-mers et reads le fichier de reads. Après l'utilisation de cette fonction on optient un dictionnaire au format { Read : { K-mer : position } }
    <h3> get_kmer_position( k , reads ) </h3>





  





