#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pickle
import pandas as pd
import fill_matrix_opti as matrix
# import getopt
# import sys
import argparse
import time


# Ce programme a pour but de mapper des reads sur une séquence de référence afin d'avoir un fichier vcf de la forme :
# POSITION / REFERENCE / ALTERNATIF / ABONDANCE

def get_my_fmi(index):
    """
    Fonction qui lit le fichier contenant le FM index et le retourne pour pouvoir l'utiliser

    :param index: fichier contenant le FM index
    :return: le FM index
    """
    with open(index, "rb") as f1:
        my_fmi = pickle.load(f1)
        return my_fmi


def left_first(alpha: chr, k: int, n: {}) -> int:
    """
    Fonction pour connaitre la position des suffixe commençant par un caractère

    :param alpha: Caractère ($, A, C, G ou T) dont on veut connaitre la ligne
    :param k: Rang du caractère demandé
    :param n: Rang de chaque caractère
    :return:  renvoie la ligne l telle que sa[l] est la position du k ième suffixe (ordonnés dans l’ordre lexico)
              débutant par le caractère alpha
    """
    # message d'erreur si alpha n'existe plus :
    assert k <= n[alpha], f"Cannot ask for the {k}^th {alpha}, it does not exist"
    # return la ligne pour chaque caractère :
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
    raise ValueError(f"Character {alpha} not in the bwt")  # pour n'avoir que les caractères de la bwt


def get_down(bwt: str, alpha: chr, start: int, stop: int) -> int:
    """
    Détecte la première occurrence d'alpha dans la BWT pour i dans [start, stop].

    A partir du départ, descendre dans la bwt aussi longtemps que bwt[line] != alpha et line <= stop
      - si bwt[line] == alpha, renvoie la ligne correspondante
      - si ligne > stop : renvoie -1

    :param bwt: burrows wheeler présente dans le FMI
    :param alpha: caractère ($, A, C, G ou T)
    :param start:  ligne de départ
    :param stop: ligne de fin
    :return: la ligne correspondante au caractère alpha pour bwt[line] != alpha et line <= stop
    """
    line = start
    while line <= stop:
        if bwt[line] == alpha:
            return line
        line += 1
    return -1


def get_up(bwt: str, alpha: chr, start: int, stop: int) -> int:
    """
        Détecte la première occurrence d'alpha dans la BWT pour i dans [start, stop].

    A partir du départ, descendre dans la bwt aussi longtemps que bwt[line] != alpha et line >= start
      - si bwt[line] == alpha, renvoie la ligne correspondante
      - si ligne > stop : renvoie -1

    :param bwt: burrows wheeler présente dans le FMI
    :param alpha: caractère ($, A, C, G ou T)
    :param start: ligne de départ
    :param stop: ligne de fin
    :return: la ligne correspondant au caractère alpha pour bwt[line] != alpha et line >= start
    """
    line = stop
    while line >= start:
        if bwt[line] == alpha:
            return line
        line -= 1
    return -1


def get_occurrences(pattern: str, bwt: str, n: {}, r: [], sa: [int]) -> []:
    """
    Retourne les positions des occurences du pattern dans la séquence de référence à l'aide de bwt, sa , n , r.
    On obtient en sortie une liste des positions des occurences.
    
    :param pattern: séquencce à tester
    :param bwt: Transformée de BW : my_fmi[0]
    :param n: nombre de chaque caractere
    :param r: rang de chaque caractere
    :param sa: suffix array
    :return: liste d'occurence du pattern dans la bwt
    """
    start = 0
    stop = len(bwt)-1
    
    # lit le pattern de droite à gauche
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


def bwt_2_seq(bwt: str, n: {}, r: []) -> str:
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


def get_kmer_position(k: int, reads: str, index) -> {}:
    """
    Retourne un dictionnaire contenant les kmer et leurs positions possibles pour chaque read :
    { Read : { k-mer : [position] } }

    :param k: longueur du kmer
    :param reads: fichier fasta de reads
    :param index: fichier contenant la FM index
    :return: un dictionnaire { Read : { k-mer : [position] } }
    """
    my_fmi = get_my_fmi(index)  # stockage du FM index
    with open(reads, 'r') as reads_file:  # lecture du fichier de reads
        kmer_position = {}  # dictionnaire vide
        for line in reads_file:
            if line[0] != ">":
                read_line = line.strip()  # retire les lignes qui ne sont pas des reads
                start = -1  # -1 car sinon ça commence à la deuxième lettre du read
                end = k - 1  # longueur du k-mer
                kmer_position[read_line] = {}  # dictionnaire avec les reads
                for kmers in range(start, len(read_line), k):  # recherche du k-mer sur le genome de reference
                    while end <= len(read_line)-1:
                        start += 1
                        end += 1
                        occ = []
                        kmer_position[read_line][read_line[start:end]] = occ  # dictionnaire kmer : position
                        if len(read_line[start:end]) == k:  # ne garde que les kmer faisant la taille demandée
                            occ = get_occurrences(read_line[start:end], my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])
                            kmer_position[read_line][read_line[start:end]] += occ
                            # ajout des occurences dans le dictionnaire
    return kmer_position


def mapping(ref, index, reads: str, k: int, max_hamming: int, min_abundance: int, out_file):
    """
    Map les différents reads sur la séquence de référence et retourne un fichier vcf avec la position des SNPs,
    les allèles de référence et alternatif et l'abondance du snp.

    :param ref: fichier contenant la séquence de référence
    :param index: fichier contenant le FM index
    :param reads: fichier contenant les reads
    :param k: longueur du kmer pour l'ancrage
    :param max_hamming: maximum de substitutions
    :param min_abundance: minimum d'abondance de substitutions
    :param out_file: fichier vcf de sortie
    :return: Un fichier vcf contenant des informations sur la position et l'abondance des snps
    """
    # stockage du FM index
    with open(index, "rb") as f2:
        fmi = pickle.load(f2)

    # CREATION DU DICTIONNAIRE {Read:{k-mer:[position]}} contenant les positions des occurences k-mers pour chaque read
    dict_kmer_position = get_kmer_position(k, reads, index)
    
    # Recherche des MEILLEURES POSITIONS D'ALIGNEMENT
    dict_final = {}  # Créer un dictionnaire qui contiendra les meilleures positions d'alignement pour chaque read.
    key_dict = list(dict_kmer_position.keys())  # stockage des différentes reads à alignés
    sequence_initiale = bwt_2_seq(fmi[0], fmi[2], fmi[3])  # récupération de la séquence initiale
    
    i = 0  # permet de déterminer la read que l'on aligne sur le génome
    for kmer in dict_kmer_position.values():  # lecture de chaque dictionnaire associé aux reads
        #  lecture de chaque position associée aux k-mers
        pos_r = 0  # position du k-mer sur le read. Le premier kmer commence toujours sur le premier nucléotide du read
        score = 0  # score de l'alignement entre le read et son ancrage sur le génome de référence
        for pos in kmer.values(): 
            if pos:  # si il y a bien une position associé au k-mer
                for position in pos:  # pour toute les positions du k-mer
                    read = key_dict[i]  # read numéro i

                    # CREATION D'UNE MATRICE D'ALIGNEMENT initialisée à 0 entre le read et son ancrage sur le génome
                    # La position du read sur le génome est déterminée par la position de l'alignement du k-mer sur le
                    # génome et la position du k-mer sur le read (pos_r).
                    # Le read s'alignera de la position "position - pos_r" à cette même position + la  longueur du read.
                    dm = matrix.DynamicMatrix(read,
                                              sequence_initiale[(position - pos_r):(position - pos_r + len(read))])
                    fill_mat = dm.fillH()  # Calcul du score d'alignement

                    # VERIFICATION DU SCORE
                    # si le score est plus grand que le précédent et respecte le nombre max de substitutions autorisé
                    # alors le score prédédent et la meilleur position d'alignement remplacé
                    if score < fill_mat and fill_mat >= (len(read) - max_hamming):
                        position_finale = (position - pos_r)
                        score = fill_mat
                        dict_final[key_dict[i]] = position_finale 
            pos_r += 1
        i += 1
    # Dictionnaire dict_final {Read: Position} ne contenant que la meilleure position d'alignement pour chaque read.
        
    # CREATION DE LA TABLE VCF
    column = ["POS", "REF", "ALT", "ABUNDANCE"] 
    data_vcf = pd.DataFrame(columns=column)
    
    # REMPLISSAGE DE LA TABLE VCF
    for cle, valeur in dict_final.items():
        pos_read = 0
        for read2, refer in zip(cle, sequence_initiale[valeur:(valeur+100)]):
            if read2 != refer:
                if data_vcf.loc[data_vcf['POS'] == (valeur + pos_read)].empty is False:
                    # Si la position détectée existe déjà
                    if data_vcf.loc[data_vcf['POS'] == (valeur + pos_read)]['ALT'].empty is True:
                        # Si l'allèle alternatif à cette position n'existe pas
                        # Nouvelle ligne dans le tableau avec cet allèle alternatif à cette position
                        new_line = {'POS': (valeur + pos_read), 'REF': refer, 'ALT': read2, 'ABUNDANCE': 1}
                        data_vcf = data_vcf.append(new_line, ignore_index=True)
                    else:
                        # +1 à l'abondance si la position existe déjà et si l'allèle alternatif existe déjà
                        if data_vcf['ABUNDANCE'][data_vcf['POS'] == (valeur + pos_read)][data_vcf['ALT'] == read2].empty is False:
                            data_vcf['ABUNDANCE'][data_vcf['POS'] == (valeur + pos_read)] += 1
                else:
                    # Si la position détectée n'existe pas encore : nouvelle ligne avec les valeurs
                    new_line = {'POS': (valeur + pos_read), 'REF': refer, 'ALT': read2, 'ABUNDANCE': 1}
                    data_vcf = data_vcf.append(new_line, ignore_index=True)
            pos_read += 1

    with open(out_file, 'w') as vcf:
        # Ecriture des 1eres lignes du fichier :
        vcf.write("#REF: " + ref + "\n""#READS: " + reads +
                  "\n"'#K: ' + str(k) + '\n''#MAX_SUBST: ' +
                  str(max_hamming) + '\n''#MIN_ABUNDANCE: ' + str(min_abundance) + '\n')
    data_vcf = data_vcf.sort_values(by='POS')  # Trier en ordre croissant les positions
    # Ajout du tableau dans le fichier en prenant en compte le minimum d'abondance :
    data_vcf[data_vcf.ABUNDANCE >= min_abundance].to_csv(out_file, index=None, sep='\t', mode='a', header=False)


tps1 = time.time()
# mapping('smallMappingTest/reference.fasta', 'dumped_index.dp', 'smallMappingTest/reads.fasta', 19, 5, 1, 'snps10.vcf')
mapping('coli/ecoli_sample.fasta', 'dumped_index3.dp', 'coli/ecoli_mutated_reads_1000.fasta', 19, 6, 3, 'snps12.vcf')
tps2 = time.time()
print(tps2 - tps1)

parser = argparse.ArgumentParser()
parser.add_argument("--ref", help="Génome de référence", required=True)
parser.add_argument("--index", help="FM index", required=True)
parser.add_argument("--reads", help="Fichier de reads", required=True)
parser.add_argument("-k", help="Kmer", type=int, required=True)
parser.add_argument("--max_hamming", help="Maximum de substitution", type=int, required=True)
parser.add_argument("--min_abundance", help="Minimum abondance", type=int, required=True)
parser.add_argument("--out", help="Fichier vcf de sortie", required=True)

args = parser.parse_args()
# print(get_kmer_position(args.k, args.reads))
map(args.ref, args.index, args.reads, args.k, args.max_hamming, args.min_abundance, args.out)

# python map.py --ref smallMappingTest/reference.fasta --index dumped_index.dp --reads smallMappingTest/reads.fasta -k 19 --max_hamming 5 --min_abundance 3 --out snps.vcf
