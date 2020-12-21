#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pickle
import getopt
import sys
import time


# Ce programme a pour but de mapper des reads sur une séquence de référence en relevant les différentes substitutions 
# que peut contenir l'alignement et de sortir ses différences dans un fichier SNPs au format vcf :
# POSITION / REFERENCE / ALTERNATIF / ABONDANCE

def get_my_fmi(index):
    """
    Lit le fichier contenant le FM index et le retourne pour pouvoir l'utiliser

    :param index: fichier contenant le FM index
    :return: le FM index
    """
    with open(index, "rb") as f1:
        my_fmi = pickle.load(f1)
        return my_fmi


def left_first(alpha: chr, k: int, n: {}) -> int:
    """
    Permet de connaitre la position des suffixes commençant par un caractère donné

    :param alpha: Caractère ($, A, C, G ou T) dont on veut connaitre la ligne
    :param k: Rang du caractère demandé
    :param n: Rang de chaque caractère
    :return:  renvoie la ligne l telle que sa[l] est la position du k ième suffixe (ordonnés dans l’ordre lexico)
              débutant par le caractère alpha
    """
    # message d'erreur si alpha n'existe pas :
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

    A partir du départ, descends dans la bwt aussi longtemps que bwt[line] != alpha et line <= stop
      - si bwt[line] == alpha, renvoie la ligne correspondante
      - si ligne > stop : renvoie -1

    :param bwt: transformee de Burrows Wheeler
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

    A partir du départ, descends dans la bwt aussi longtemps que bwt[line] != alpha et line >= start
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
    
    :param pattern: séquence à tester
    :param bwt: Transformée de BW : my_fmi[0]
    :param n: nombre de chaque caractere
    :param r: rang de chaque caractere
    :param sa: suffix array
    :return: liste d'occurences du pattern dans la bwt
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
    Retourne la sequence initiale à partir de la BWT
    
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


def reverse_complement(seq):
    """
    Permet d'obtenir le reverse complément de la séquence passée en paramètre

    :param seq: Séquence de nucléotides
    :return: Reverse complement de seq
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in seq[::-1]])


def get_kmer_position(k: int, reads: str, index) -> {}:
    """
    Retourne un dictionnaire pour contenant les kmer et leurs positions possibles pour chaque read
    et son reverse complément :
    {Rang du read: {k-mer: [position]}} ainsi qu'une liste de reads et leurs reverse compléments

    :param k: longueur du kmer
    :param reads: fichier fasta de reads
    :param index: fichier contenant la FM index
    :return: un dictionnaires {Rang du read: {k-mer: [position]}} et une liste de reads
    """
    my_fmi = get_my_fmi(index)  # stockage du FM index
    with open(reads, 'r') as reads_file:  # lecture du fichier de reads

        kmer_position = {}  # dictionnaire vide
        kmer_reverse_position = {}
        nb_read = 0
        pos_sens = 1
        pos_rev = 2
        readfasta = []  # liste de reads vide

        for line in reads_file:
            if line[0] != ">":
                read_line = line.strip()  # retire les lignes qui ne sont pas des reads
                read_line_complement = reverse_complement(read_line)
                readfasta.append(read_line)  # ajout du read
                readfasta.append(read_line_complement)  # ajout du reverse complement
                start = -1  # -1 car sinon commence à la deuxième lettre du read
                end = k - 1  # longueur du k-mer
                kmer_position[pos_sens] = {}  # dictionnaire avec les reads
                kmer_reverse_position[pos_rev] = {}
                # read initial
                for kmers in range(start, len(read_line), k):  # recherche du k-mer sur le genome de reference
                    while end <= len(read_line)-1:
                        start += 1
                        end += 1
                        occ = []
                        occ_rev = []
                        kmer_position[pos_sens][read_line[start:end]] = occ  # dictionnaire kmer : position
                        kmer_reverse_position[pos_rev][read_line_complement[start:end]] = occ_rev
                        if len(read_line[start:end]) == k or len(read_line_complement[start:end]) == k:
                            # ne garde que les kmer faisant la taille demandée
                            occ = get_occurrences(read_line[start:end], my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])
                            kmer_position[pos_sens][read_line[start:end]] += occ
                            occ_rev = get_occurrences(read_line_complement[start:end], my_fmi[0],
                                                      my_fmi[2], my_fmi[3], my_fmi[1])
                            kmer_reverse_position[pos_rev][read_line_complement[start:end]] += occ_rev
                pos_sens += 2
                pos_rev += 2
                nb_read += 1
        kmer_position = {**kmer_position, **kmer_reverse_position}  # fusion des deux dictionnaires
        kmer_position = dict(sorted(kmer_position.items()))
        return kmer_position, readfasta


def fill_vcf(mat, dict_final, sequence_initiale, list_read):
    """
    RempliT une matrice au format vcf : POSITION / REFERENCE / ALTERNATIF / ABONDANCE à
    partir d'un dictionnaire Rang-read:position et sa liste de reads associé ainsi que d'une sequence de référence.
   
    
    :param mat: matrice [[],[],[],[]]
    :param dict_final: dictionnaire Rang-read:position
    :param sequence_initiale: Sequence de référence
    :param list_read : Liste de reads +/- ordonnée
    :return: la matrice mat remplit avec les SNPs
    
    """
    for cle, valeur in dict_final.items():  # Parcours du dictionnaire Rang-read : position
        pos_read = 0

        for read, refer in zip(list_read[(cle-1)], sequence_initiale[valeur:(valeur + 100)]):
            # Parcours de chaque nucléotide du read et de la sequence où le read s'aligne
            # Si il y a une substitution.
            if read != refer:
                # Si la substitution a déja été pris en compte
                if (valeur + pos_read) in mat[0]:
                    y = mat[0].index((valeur + pos_read))
                    mat[3][y] += 1
                # Si la substitution n'a pas encore été pris en compte
                else:
                    mat[0].append((valeur + pos_read))
                    mat[1].append(refer)
                    mat[2].append(read)
                    mat[3].append(1)
            pos_read += 1
    return mat


def order_vcf(tab_vcf):
    """
    Permet de mettre en ordre en fonction de la position une table vcf au format :
    POSITION / REFERENCE / ALTERNATIF / ABONDANCE
    
    :param tab_vcf: matrice [[],[],[],[]] au format vcf
    :return: une nouvelle matrice rangé par position
    
    """

    pos_sorted = sorted(tab_vcf[0])
    new_mat = [[], [], [], []]
    for i in pos_sorted:
        for y, u in enumerate(tab_vcf[0]):
            if u == i:
                new_mat[0].append(i)
                new_mat[1].append(tab_vcf[1][y])
                new_mat[2].append(tab_vcf[2][y])
                new_mat[3].append(tab_vcf[3][y])
    return new_mat


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

    # CREATION DU DICTIONNAIRE {Rang-read:{k-mer:[position]}}
    # contenant les positions des occurences k-mers pour chaque read
    # ainsi que la liste ordonnée des reads +/-
    kmer_position = get_kmer_position(k, reads, index)

    # RECHERCHE DES MEILLEURES POSITIONS D'ALIGNEMENT
    dict_final = {}  # Créer un dictionnaire qui contiendra les meilleures positions d'alignement pour chaque read +/-.
    sequence_initiale = bwt_2_seq(fmi[0], fmi[2], fmi[3])  # récupération de la séquence initiale

    i = 0  # permet de déterminer la read que l'on aligne sur le génome
    for kmer in kmer_position[0].values():  # lecture de chaque dictionnaire associé aux reads
        #  lecture de chaque position associée aux k-mers
        pos_r = 0  # position du k-mer sur le read. Le premier kmer commence toujours sur le premier nucléotide du read
        score = 0  # score de l'alignement entre le read et son ancrage sur le génome de référence
        for pos in kmer.values():
            if pos:  # si il y a bien une position associé au k-mer
                for position in pos:  # pour toute les positions du k-mer
                    read = kmer_position[1][i]  # read numéro i

                    # ALIGNEMENT entre le read et son ancrage sur le génome
                    # La position du read sur le génome est déterminée par la position de l'alignement du k-mer sur le
                    # génome et la position du k-mer sur le read (pos_r).
                    # Le read s'alignera de la position "position - pos_r" à cette même position + la  longueur du read.
                    # On compare chaque nucléotide du read avec celle de son angrage et on ajoute +1 au score si
                    # les deux nucléotides sont égales
                    score_ali = 0
                    for Nucread, Nucseq in zip(read,
                                               sequence_initiale[(position - pos_r):(position - pos_r + len(read))]):
                        if Nucread == Nucseq:
                            score_ali += 1

                    # VERIFICATION DU SCORE
                    # si le score est plus grand que le précédent et respecte le nombre max de substitutions autorisé
                    # alors le score prédédent et la meilleur position d'alignement remplacé
                    if score < score_ali and fill_mat >= (len(read) - max_hamming):
                        position_finale = (position - pos_r)
                        score = score_ali
                        dict_final[(i+1)] = position_finale
            pos_r += 1
        i += 1

    # CREATION DE LA TABLE VCF
    mat = [[], [], [], []]  # initialisation de la table vcf
    tab_vcf = fill_vcf(mat, dict_final, sequence_initiale, kmer_position[1])  # remplissage de la table vcf
    tab_vcf = order_vcf(tab_vcf)
    
    # ECRITURE DU FICHIER VCF
    with open(out_file, 'w') as vcf:
        # Ecriture des 1eres lignes du fichier :
        vcf.write("#REF: " + ref + "\n""#READS: " + reads +
                  "\n"'#K: ' + str(k) + '\n''#MAX_SUBST: ' +
                  str(max_hamming) + '\n''#MIN_ABUNDANCE: ' + str(min_abundance) + '\n')
        
        # Ecriture des données de la table vcf en fonction de l'abondance minimum retenu
        i = 0
        while i < len(tab_vcf[0]):
            if tab_vcf[3][i] >= min_abundance:
                vcf.write(
                    str(tab_vcf[0][i]) + '\t' + tab_vcf[1][i] + '\t' + tab_vcf[2][i] + '\t' + str(
                        tab_vcf[3][i]) + '\n')
            i += 1


if __name__ == "__main__":
    reference = ''
    index_file = ''
    read_file = ''
    k_mers = 1
    hamming = 1
    abundance = 0
    output = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:], "k:h",
                                   ["ref=", "index=", "reads=", "max_hamming=", "min_abundance=", "out="])

    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)
    for option, arg in opts:
        if option in "-h":
            print('python map.py --ref[genome_file.fa] --index[dumped_index.dp] '
                  '--reads[reads.fa] -k[k_value] --max_hamming[h_value] --min_abundance[m_value] --out snps.vcf')
            sys.exit(2)
        elif option in "--reads":
            read_file = arg
        elif option in "-k":
            k_mers = int(arg)
        elif option in "--ref":
            reference = arg
        elif option in "--index":
            index_file = arg
        elif option in "--max_hamming":
            hamming = int(arg)
        elif option in "--min_abundance":
            abundance = int(arg)
        elif option in "--ref":
            ref_file = arg
        elif option in "--out":
            output = arg

t1 = time.time()
mapping(reference, index_file, read_file, k_mers, hamming, abundance, output)
t2 = time.time()
print(t2-t1)

# python map.py --ref smallMappingTest/reference.fasta --index smallMappingTest/dumped_index_small.dp --reads
# smallMappingTest/reads.fasta -k 20 --max_hamming 5 --min_abundance 1 --out smallMappingTest/snps_ref_k19_d5.vcf
# python map.py --ref coli/ecoli_sample.fasta --index coli/dumped_index_coli.dp --reads
# coli/ecoli_mutated_reads_1000.fasta -k 20 --max_hamming 10 --min_abundance 7 --out coli/snps_coli_k20_d10_ab7.vcf
