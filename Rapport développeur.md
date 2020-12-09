Le programme se distingue en deux fichier .py distinct, le premier traite de l'indexation du génome de référence et le second du mapping des reads sur le génome de référence

Le fichier d'indexation nommé index.py a pour but de crée un index composé de la transformé de Burrows wheeler du génome de référence, ainsi que les informations permettant l'utilisation de la bwt

<h1>get_fmi(ref_fast, output_file) </h1>

<p>La fonction get_fmi(ref_fasta, output_file ) permet de crée cet index, elle prend en entrée le génome de référence au format fasta, ainsi que le nom que l'utilisateur choisis 
de prendre pour nommé l'index crée, en sortie la fonction retourne un "nom".db au format pickle contenant la bwt, le rang et le nombre de chaque caractère du génome
ainsi que le suffixe array.</p>

Cette fonction fonctionne de la façon suivante. 

Dans un premier temps la transformé de burrows wheeler va etre former à l'aide du génome de réference. La fonction get_bwt(fasta: str) prend en entré donc le génome de 
référence et permet d'avoir en sortie la bwt. 
La fonction fonctionne de la manière suivante : 

la fonction get_seq( fasta: str ) traite le fichier du génome de réference, elle saute la première ligne consititué de la ligne informative 
commençant par ">" puis stocke la séquence du génome de référence (s) et crée le suffix array (sa) à l'aide la fonction kark_sort().

Ensuite à l'aide de sa et de s la transformée de burrows wheeler est calculé.

Arrivé dans cette partie de la fonction get_fmi on a déjà le sa et la bwt, il ne manque plus que le rang et nombre pour chaque caractère.

La fonction get_r_n(bwt) va calculer à partir de la transformée burrows wheeler le rang ainsi r et le nombre n pour chaque caractère du génome de référence.
En sortie on aura 




