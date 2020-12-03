#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pickle
import pandas as pd
import numpy as np
import fill_matrix as matrix

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
    for occ in range(start, stop+1):
        res.append(sa[occ])
    return res


# get_occurrences("T", my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])


'''
def map(ref, index, reads, k, max_hamming, min_abundance, out ):
'''

def bwt_2_seq(bwt: str, n: {}, r:[]) -> str:
    """
    Fonction qui retourne la sequence initiale à partir de la BWT
    :param bwt: Transformée de BW : my_fmi[0]
    :param n: nombre de chaque caractere
    :param r: rang de chaque caractere
    :return: La séquence d'origine
    """
    sequence_reconstructed = ""
    line = 0
    while True:
        if bwt[line] == "$":
            break
        sequence_reconstructed = bwt[line] + sequence_reconstructed
        line = left_first(bwt[line], r[line], n)
    return sequence_reconstructed



def get_kmer_position(k, read_fasta):
    """
    Retourne un dictionnaire contenant les kmer et leurs positions possibles pour chaque read
    :param k: longueur du kmer
    :param read_fasta: fichier fasta de reads
    :return:
    """
    with open(read_fasta) as reads_file:
        kmer_position = {}  # dictionnaire vide
        for line in reads_file:
            if line[0] != ">":
                read_line = line.strip()  # retire les lignes qui ne sont pas des reads
                start = -1  # -1 car sinon ça commence à la deuxième lettre du read (je ne sais pas pourquoi)
                end = k - 1  # longueur du k-mer
                kmer_position[read_line] = {}  # dictionnaire avec les reads
                for kmer in range(start, len(read_line), k):  # recherche du k-mer sur le genome de reference
                    while end <= len(read_line)-1:
                        start += 1
                        end += 1
                        occ = []
                        kmer_position[read_line][read_line[start:end]] = occ  # dictionnaire kmer : position
                        if len(read_line[start:end]) == k:  # possible k-mer qui ne fait pas la taille k demandée
                            occ = get_occurrences(read_line[start:end], my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])
                            kmer_position[read_line][read_line[start:end]] += occ
    return kmer_position


dict_kmer_position = get_kmer_position(100, "smallMappingTest/reads.fasta")

sequence_initiale = bwt_2_seq(my_fmi[0], my_fmi[2], my_fmi[3])
key_dict = list(dict_kmer_position.keys())
i = 0
hamming = 50
dict_final = {}
position_finale = 0
for kmer in dict_kmer_position.values():
    pos_r = 0
    score = 0
    fillH = 0
    print("read : ", key_dict[i])
    for pos in kmer.values():
        if pos != []:
            for position in pos:
                read = key_dict[i]
                # print("read : ", read)
                print(position)
                print("sequence : ", sequence_initiale[(position - pos_r):(position - pos_r + len(read))])
                dm = matrix.DynamicMatrix(read, sequence_initiale[(position - pos_r):(position - pos_r + len(read))], 1, 0,
                                          0)
                fillH = dm.fillH(1)
                print(fillH)
                if score < fillH and fillH > (len(read) - hamming):
                    position_finale = (position - pos_r)
                    #print(position_finale)
                    score = fillH
            dict_final[key_dict[i]] = position_finale
        pos_r += 1
    i += 1
dict_final


key_read = list(dict_final.keys())
# aligner notre read sur le genome de ref en comptant les substition
# creation d'un tableau avec les 4 colonnes position / ref / alt / abondance (numpy)
        #chaque read on va
        # si mismatche existe déja dans le tableau + 1 dans l'abondance.
        # si le mismatch n'existe pas, on crée la ligne correspondant.
        # si il n'y a match suivant
    # A la fin on garde les lignes qui respecte l'abondance.


column = ["POS", "REF", "ALT", "ABUNDANCE"]

data_vcf = pd.DataFrame(columns = column)
data_vcf
import pandas as pd
data_vcf['POS']
data_vcf


sequence_initiale

for cle, valeur in dict_final.items():
    pos_read = 0
    for read, ref in cle , sequence_initiale:
        if read != ref:
            if data_vcf.loc[data_vcf['POS'] == (valeur + pos_read)].empty == False:
                if data_vcf.loc[data_vcf['ALT'][data_vcf['POS'] == (valeur + pos_read)].empty == True:
                    new_line = {'POS': (valeur + pos_read), 'REF': ref, 'ALT': read, 'ABUNDANCE': 1}
                    data_vcf.append(new_line, ignore_index=True)
                else:
                    data_vcf['ABUNDANCE'][data_vcf['POS'] == (valeur + pos_read) and data_vcf['ALL'] == read ] += 1
            else:
                new_line = {'POS':(valeur + pos_read),'REF':ref ,'ALT':read ,'ABUNDANCE':1}
                data_vcf.append(new_line,ignore_index=True)
        pos_read += 1


