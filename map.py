#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pickle

# Mapping des reads sur le génome de référence

with open('dumped_index.dp', "rb") as f1:
    my_fmi = pickle.load(f1)

# my_fmi[0] # bwt
# my_fmi[1] # sa
# my_fmi[2] # n
# my_fmi[3] # rank


def left_first(alpha: chr, k: int, n: {}) -> int:
    assert k <= n[alpha], f"Cannot ask for the {k}^th {alpha}, it does not exist"
    if alpha == "$":
        return 0
    if alpha == "A":
        return n["$"] + k - 1
    if alpha == "C":
        return n["$"] + n["A"] + k - 1
    if alpha == "G":
        return n["$"] + n["A"] + n["C"] + k - 1
    if alpha == "T":
        return n["$"] + n["A"] + n["C"] + n["G"] + k - 1
    raise ValueError(f"Character {alpha} not in the bwt")


def get_down(bwt: str, alpha: chr, start: int, stop: int) -> int:
    """
    Detects the first occurrence of alpha in bwt for i in [start, stop].

    From start go down in the bwt as long as bwt[line] != alpha and line <= stop
      - if bwt[line] == alpha, returns the corresponding line
      - if line > stop: returns -1
    """
    line = start
    while line <= stop:
        if bwt[line] == alpha:
            return line
        line += 1
    return -1


def get_up(bwt: str, alpha: chr, start: int, stop: int) -> int:
    line = stop
    while line >= start:
        if bwt[line] == alpha:
            return line
        line -= 1
    return -1


def get_occurrences(pattern: str, bwt: str, n: {}, r: [], sa: [int]) -> bool:
    """
    Returns pattern occurrences in the text coded in the bwt
    """
    start = 0
    stop = len(bwt)-1
    # read the pattern from right to left
    for pos_pattern in range(len(pattern)-1, -1, -1):
        current_char = pattern[pos_pattern]
        new_start = get_down(bwt, current_char, start, stop)
        if new_start == -1:
            return []
        new_stop = get_up(bwt, current_char, start, stop)
        start = left_first(bwt[new_start], r[new_start], n)
        stop = left_first(bwt[new_stop], r[new_stop], n)
    res = []
    for i in range(start, stop+1):
        res.append(sa[i])
    return res


get_occurrences("T", my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])


'''
def map(ref, index, reads, k, max_hamming, min_abundance, out ):
    for line in reads:
        if line[0] != ">":
            read = line
            start = 0
            end = k - 1
        for i in line:
            while end <= len(read):
                occ[i] = get_occurrences(read[start:end], my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])
                start += start
                end += end
                i += i
'''


reads = open("smallMappingTest/reads.fasta") # fichier de reads
k = 30 # k-mer
occ = [] # list des position des k-mer du read sur le génome de reference
for line in reads:
    if line[0] != ">":
        read = line.strip()  # on garde toutes les lignes du fichier fasta qui ne commence pas par ">" (ce sont nos reads)
        print("Read : ", read)
        start = -1  # -1 car sinon ça commence à la deuxième lettre du read (je ne sais pas pourquoi)
        end = k - 1 # longueur du k-mer pour stoper
        for i in range(start, len(read), end):  # on recherche le k-mer sur le genome de reference en changeant de pas sur toute la longueur du read (un pas = longueur du kmer)
            while end <= len(read):
                start += 1
                end += 1
                if len(read[start:end]) == k: # car sinon on peut avoir le dernier k-mer qui ne fait pas la taille k demandée
                    occ = ([get_occurrences(read[start:end], my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])]) # recherche des occurrences du kmer
                    if occ != [[]]: # pour ne pas afficher les listes vides
                        print("Le k-mer : ", read[start:end], "match sur la position : ", occ)


