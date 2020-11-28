#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Indexer le genome de reference fourni en fasta (FM-index)

import tools_karkkainen_sanders as tks
import getopt
import sys
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
                s = line.strip()+"$"
                sa = tks.simple_kark_sort(s)
                return s, sa


# Test de la fonction
# get_seq("smallMappingTest/reference.fasta")

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


'''
test :
s = get_seq("smallMappingTest/reference.fasta")[0]
sa = get_seq("smallMappingTest/reference.fasta")[1]
'''

# Test de la fonction :
"""
bwt = get_bwt("smallMappingTest/reference.fasta")
print("i\tsa[i]\tbwt[i]\tF")
for i in range(len(s)):
    print(f"{i}\t{sa[i]}\t{bwt[i]}\t{s[sa[i]]}")
"""


def get_n(fasta: str) -> {}:
    n = {"$": 0, "A": 0, "C": 0, "G": 0, "T": 0}  # key = letter, value = number of occurrences
    for letter in get_bwt(fasta):
        n[letter] += 1
    return n


def get_r(fasta: str) -> {}:
    n = {}  # key = letter, value = number of occurrences
    r = []  # for each i: rank of the i^th value in bwt
    for letter in get_bwt(fasta):
        if letter not in n:
            n[letter] = 0
        n[letter] += 1
        r.append(n[letter])
    return r


'''
n = get_n("smallMappingTest/reference.fasta")
print(n)
r = get_r("smallMappingTest/reference.fasta")
print(r)
print(n)
'''


def get_fmi(ref_fasta, output_file):
    """
    Fonction qui creer le FMindex avec :
        L : BWT
        sa :
        n : nombre de chaque caractere
        r : rang de chaque caractere
    :param ref_fasta : sequence fasta de reference
    :param output_file : fichier de sortie contenant le FMI
    :return:
    """
    bwt = ""
    sa = get_seq(ref_fasta)[1]
    n = get_n(ref_fasta)
    r = get_r(ref_fasta)
    for i in range(len(get_seq(ref_fasta)[0])):
        bwt += get_bwt(ref_fasta)[i]
    with open(output_file, "wb") as f1:
        pickle.dump((bwt, sa, n, r), f1)
    return bwt, sa, n, r



# verification fonction get_fmi
my_fmi = get_fmi("smallMappingTest/reference.fasta",'dumped_index.dp')
my_fmi[0] # bwt
my_fmi[1] # sa
my_fmi[2] # n
my_fmi[3] # r 


'''
# Verification du fichier :
verif = None
print(verif)
with open('dumped_index.dp', "rb") as f1:
    verif = pickle.load(f1)
type(verif)
print(verif)
'''

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
            print('index.py --ref [genome_file.fa] --out [dumped_index.dp]')
            sys.exit()
        elif option in ("-r", "--ref"):
            ref_file = arg
        elif option in ("-o", "--out"):
            out_file = arg
    print('Reference file is ', ref_file)
    print('Output file is ', out_file)

    get_fmi(ref_file, out_file)

# python index.py --ref smallMappingTest/reference.fasta --out dumped_index.dp
