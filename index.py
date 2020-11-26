#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Indexer le genome de reference fourni en fasta (FM-index)

import tools_karkkainen_sanders as tks
# import getopt
# import sys
import pickle


def get_seq(fasta: str):
    """
    Fonction pour avoir la sequence sans prendre en compte la premiere ligne du fichier fasta commençant par ">".
    Elle va ouvrir le fichier entre en parametre puis lire la deuxieme ligne et en faire le suffixe array

    :param fasta: sequence fasta de reference
    :return: La sequence du genome de reference et le suffixe array de cette sequence
    """
    with open(fasta) as fasta_file:
        for line in fasta_file:
            if line[0] != ">":
                s = line.strip()
                sa = tks.simple_kark_sort(s)
                return s, sa


# Test de la fonction
get_seq("smallMappingTest/reference.fasta")

# Pour avoir la sequence :
# get_seq("smallMappingTest/reference.fasta")[0]
# Pour avoir la suffix array :
# get_seq("smallMappingTest/reference.fasta")[1]


def get_bwt(fasta: str) -> str:
    """
    Fonction pour obtenir la transformee de Burrows Wheeler à partir du fichier de reference en entree.
    Grace à la fonction get_seq() on obtient la suffixe array qui va permettre la BWT

    :param fasta: sequence fasta de reference
    :return: transformee de BW
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


s = get_seq("smallMappingTest/reference.fasta")[0]
sa = get_seq("smallMappingTest/reference.fasta")[1]

# Test de la fonction :
bwt = get_bwt("smallMappingTest/reference.fasta")
print("i\tsa[i]\tbwt[i]\tF")
for i in range(len(get_seq("smallMappingTest/reference.fasta")[0])):
    print(f"{i}\t{get_seq('smallMappingTest/reference.fasta')[1][i]}\t{bwt[i]}\t{s[sa[i]]}")


def get_fmi(ref_fasta):
    """
    Fonction qui creer le FMindex avec :
        L : BWT
        F : Ordre alphabetique
    :param ref_fasta: sequence fasta de reference
    :return:
    """
    for i in range(len(get_seq(ref_fasta)[0])):
        L = get_bwt(ref_fasta)[i]
        F = get_seq(ref_fasta)[0][get_seq(ref_fasta)[1][i]]
        return L, F


my_fmi = get_fmi("smallMappingTest/reference.fasta")

# Ecriture du FMI dans un fichier :
with open('dumped_index.dp', "wb") as f1:
    pickle.dump(my_fmi, f1)

# Verification du fichier :
verif = None
print(verif)
with open('dumped_index.dp', "rb") as f1:
    verif = pickle.load(f1)
type(verif)
print(verif)  # il n'y a que 1 tuple donc probleme dans la fonction get_fmi()
