#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Indexer le génome de référence fourni en fasta (FM-index)

import tools_karkkainen_sanders as tks
# import getopt
# import sys
import pickle

# import du fichier fasta de référence


def get_seq(fasta: str):
    """
    Fonction pour avoir la sequence sans prendre en compte la premiere ligne du fichier fasta commençant par ">".
    Elle va ouvrir le fichier entré en paramètre puis lire la deuxième ligne et en faire le suffixe array

    :return: La séquence du génome de référence et le suffixe array de cette séquence
    """
    with open(fasta) as fasta_file:
        for line in fasta_file:
            if line[0] != ">":
                s = line.strip()
                sa = tks.simple_kark_sort(s)
                return s, sa


# Test de la fonction
get_seq("smallMappingTest/reference.fasta")

# Pour avoir la séquence :
# get_seq("smallMappingTest/reference.fasta")[0]
# Pour avoir la suffix array :
# get_seq("smallMappingTest/reference.fasta")[1]


def get_bwt(fasta: str) -> str:
    """
    Fonction pour obtenir la transformée de Burrows Wheeler à partir du fichier de référence en entrée.
    Grace à la fonction get_seq() on obtient la suffixe array qui va permettre la BWT

    :return: transformée de BW
    """
    bwt = ""
    s = get_seq(fasta)[0]
    sa = get_seq(fasta)[1]
    for i in range(len(sa)):
        if sa[i] == 0:
            bwt += "$"
        else:
            bwt += s[sa[i] - 1]
    return bwt


# Test de la fonction :
bwt = get_bwt("smallMappingTest/reference.fasta")
for i in range(len(get_seq("smallMappingTest/reference.fasta")[0])):
    print(f"{i}\t{get_seq('smallMappingTest/reference.fasta')[1][i]}\t{bwt[i]}")

# Ecriture de la bwt dans un fichier
pickle.dump(bwt, open("dumped_index.dp", "wb"))
