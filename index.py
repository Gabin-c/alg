#!/usr/bin/env python
# -*- coding: utf-8 -*-

import tools_karkkainen_sanders as tks
import getopt
import sys
import pickle


def get_seq(fasta: str):
    """
    Permet d'avoir la sequence sans prendre en compte la premiere ligne du fichier fasta commençant par ">".
    Elle va ouvrir le fichier entré en paramètre puis lire la deuxième ligne et en faire le suffixe array

    :param fasta: sequence fasta de reference
    :return: La sequence du genome de reference et le suffixe array de cette sequence
    """
    with open(fasta) as fasta_file:  # ouverture du fichier fasta
        for line in fasta_file:  # lecture du fichier fasta
            if line[0] != ">":  # Première ligne commençant pas ">" ignoré
                s = line.strip()+"$"  # stockage de la séquence
                sa = tks.simple_kark_sort(s)  # stockage du suffix array sa
                return s, sa


def get_bwt(fasta: str):
    """
    Permet d'obtenir la transformée de Burrows Wheeler à partir du fichier de référence en entrée.
    Grâce à la fonction get_seq() on obtient la suffixe array qui va permettre de construire la BWT

    :param fasta: séquence fasta de référence
    :return: transformée de BW, suffix array
    """
    bwt = "" 
    sequence = get_seq(fasta)  # appel de get_seq() pour obtenir s et sa
    s = sequence[0]
    sa = sequence[1]
    for i in range(len(sa)):  # création de la burrows wheeler
        if sa[i] == 0: 
            bwt += "$"
        else:
            bwt += s[sa[i] - 1] 
    return bwt, sa


def get_r_n(bwt):
    
    """
    Permet d'obtenir à partir de la transformée de burrows wheeler, le dictionnaire n du nombre d'occurrences pour
    chaque nucléotides et la liste r des rangs de chaque caractères dans la séquence de références.

    :param bwt: la transformée de burrows wheeler
    :return: liste des rangs r de chaque nucléotides dans la séquence de références, dictionnaire n du nombre
            d'occurrences de chaque nucléotide
    """
    n = {"$": 0, "A": 0, "C": 0, "G": 0, "T": 0}  # clé = nucléotides + $, valeur = nombre d'occurrences
    r = []  # pour chaque lettre son rang dans la séquence de référence
    for letter in bwt:  # implémentation de n et r
        if letter not in n:
            n[letter] = 0
        n[letter] += 1
        r.append(n[letter])
    return r, n


def get_fmi(ref_fasta, output_file):
    """
    Créer un FM-index avec :
        L : BWT
        sa : suffix array 
        n : nombre de chaque caractère
        r : rang de chaque caractère
    :param ref_fasta : sequence fasta de reference
    :param output_file : fichier de sortie contenant le FMI
    :return: 
    """
    res_bwt = get_bwt(ref_fasta)  # obtention de : transformée de burrows wheeler + suffix array
    bwt = res_bwt[0]  # transformée de burrows wheeler
    sa = res_bwt[1]  # suffix array
    res_n_r = get_r_n(res_bwt[0])  # obtention de n et r
    n = res_n_r[1]  # n
    r = res_n_r[0]  # r

    with open(output_file, "wb") as f1:  # stockage de bwt, sa , r et n dans un pickle du nom de "output_file"
        pickle.dump((bwt, sa, n, r), f1) 
    return bwt, sa, n, r


if __name__ == "__main__":
    ref_file = ''
    out_file = ''
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "r:o:h", ["ref=", "out="])
    except getopt.GetoptError as err:
        print('index.py --ref [genome_file.fa] --out [dumped_index.dp]')
        sys.exit(2)

    for option, arg in opts:
        if option in "-h":
            print('Usage : \n index.py --ref [genome_file.fa] --out [dumped_index.dp]')
            sys.exit()
        elif option in ("-r", "--ref"):
            ref_file = arg
        elif option in ("-o", "--out"):
            out_file = arg

    get_fmi(ref_file, out_file)
    print('Reference file is ', ref_file)
    print('Output file is ', out_file)


