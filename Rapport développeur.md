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
<h4> get_r_n(bwt) <h4>
  
La fonction get_r_n(bwt) va calculer à partir de la transformée burrows wheeler le rang r et le nombre n pour chaque caractère du génome de référence.

Dans un premier temps la fonction initialise un dictionnaire n avec les 4 nucléotides et une liste r pour stocké les rangs.
Suite à ça, la bwt va être parcouru et va implémenter au fur et à mesure notre dictionnaire n et notre list r. On obtiendra en sortie un dictionnaire n et une liste r, nos deux derniers éléments de notre FM index.

<h3> Sortie get_fmi() </h3>

A la fin des deux étapes les 4 éléments bwt, sa , n et r sont stocké dans le FM index et on obtient à la sortie du programme un FMindex.db contenant les 4 élements.
  





