#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pickle
import fill_matrix as matrix
import getopt
import sys
import time


# Ce programme a pour but de mapper des reads sur une séquence de référence en relevant les différentes substitutions 
# que peut contenir l'alignement et de sortir ses différences dans un fichier SNPs au format vcf :
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

get_my_fmi('smallMappingTest/dumped_index_small.dp')
def left_first(alpha: chr, k: int, n: {}) -> int:
    """
    Fonction qui permet de connaitre la position des suffixe commençant par un caractère

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

    A partir du départ, descends dans la bwt aussi longtemps que bwt[line] != alpha et line <= stop
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

def reverse_complement(seq):
    """
    Fonction qui permet d'obtenir le reverse complément de la séquence passer en paramètre
    
    :param seq: Sequence de nucléotides
    :return: Reverse complement de seq
    """
    alt_map = {'ins': '0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases


def get_kmer_position(k: int, reads: str, index) -> {}:
    """
    Retourne un dictionnaire pour chaque brin contenant les kmer et leurs positions possibles pour chaque read :
    {Read: {k-mer: [position]}}

    :param k: longueur du kmer
    :param reads: fichier fasta de reads
    :param index: fichier contenant la FM index
    :return: deux dictionnaires {Read: {k-mer: [position]}}
    """
    my_fmi = get_my_fmi(index)  # stockage du FM index
    with open(reads, 'r') as reads_file:  # lecture du fichier de reads

        kmer_position = {}  # dictionnaire vide
        kmer_reverse_position = {}
        for line in reads_file:
            if line[0] != ">":
                read_line = line.strip()  # retire les lignes qui ne sont pas des reads
                read_line_complement = reverse_complement(read_line)
                start = -1  # -1 car sinon commence à la deuxième lettre du read
                start_reverse = -1
                end = k - 1  # longueur du k-mer
                end_reverse = k - 1
                kmer_position[read_line] = {}  # dictionnaire avec les reads
                kmer_reverse_position[read_line_complement] = {}
                # read initial
                for kmers in range(start, len(read_line), k):  # recherche du k-mer sur le genome de reference
                    while end <= len(read_line)-1:
                        start += 1
                        end += 1
                        occ = []
                        #occ_rev = []
                        kmer_position[read_line][read_line[start:end]] = occ  # dictionnaire kmer : position
                        #kmer_reverse_position[read_line_complement][read_line_complement[start:end]] = occ_rev
                        if len(read_line[start:end]) == k:  # ne garde que les kmer faisant la taille demandée
                            occ = get_occurrences(read_line[start:end], my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])
                            kmer_position[read_line][read_line[start:end]] += occ
                            # ajout des occurences dans le dictionaire
                        #elif len(read_line_complement[start_reverse:end_reverse]) == k:  # ne garde que les kmer faisant la taille demandée
                            #occ_rev = get_occurrences(read_line_complement[start:end], my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])

                            #kmer_reverse_position[read_line_complement][read_line_complement[start:end]] += occ_rev
                            # ajout des occurences dans le dictionaire
                #  read reverse complement 
                for kmers in range(start_reverse, len(read_line_complement), k):  # recherche du k-mer sur le genome de reference
                    while end_reverse <= len(read_line_complement)-1:
                        start_reverse += 1
                        end_reverse += 1
                        occ = []
                        kmer_reverse_position[read_line_complement][read_line_complement[start_reverse:end_reverse]] = occ  # dictionnaire kmer : position
                        if len(read_line_complement[start_reverse:end_reverse]) == k:  # ne garde que les kmer faisant la taille demandée
                            occ = get_occurrences(read_line_complement[start_reverse:end_reverse], my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])

                            kmer_reverse_position[read_line_complement][read_line_complement[start_reverse:end_reverse]] += occ
                            # ajout des occurences dans le dictionnaire

        return kmer_position, kmer_reverse_position


def fill_vcf(mat, dict_final, sequence_initiale):
    """
    Fonction qui permet le remplissage d'une matrice au format vcf : POSITION / REFERENCE / ALTERNATIF / ABONDANCE à partir d'un dictionnaire read:position et de d'une sequence de référence.
   
    
    :param mat: matrice [[],[],[],[]]
    :param dict_final: dictionnaire read:position
    :param sequence_intiale: Sequence de référence
    
    :return: la matrice mat remplit avec les SNPs
    
    """
    for cle, valeur in dict_final.items(): # Parcours du dictionnaire read:position
        pos_read = 0
        for read, refer in zip(cle, sequence_initiale[valeur:(valeur + 100)]): # Parcours de chaque nucléotide du read et de la sequence où le read s'aligne
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
    Fonction permettant de mettre en ordre en fonction de la position une table vcf au format : POSITION / REFERENCE / ALTERNATIF / ABONDANCE
    
    :param tab_vcf: matrice [[],[],[],[]] au format vcf
    :return: une nouvelle matrice rangé par position
    
    """

    pos_sorted = sorted(tab_vcf[0])
    new_mat = [[],[],[],[]]
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

    # CREATION DU DICTIONNAIRE {Read:{k-mer:[position]}} contenant les positions des occurences k-mers pour chaque read
    dict_kmer_position = get_kmer_position(k, reads, index)

    # RECHERCHE DES MEILLEURES POSITIONS D'ALIGNEMENT
    dict_final_sens = {}  # Créer un dictionnaire qui contiendra les meilleures positions d'alignement pour chaque read.
    dict_final_reverse = {}
    key_dict_sens = list(dict_kmer_position[0].keys())  # stockage des différentes reads à alignés
    key_dict_reverse = list(dict_kmer_position[1].keys())
    sequence_initiale = bwt_2_seq(fmi[0], fmi[2], fmi[3])  # récupération de la séquence initiale

    i = 0  # permet de déterminer la read que l'on aligne sur le génome
    for kmer in dict_kmer_position[0].values():  # lecture de chaque dictionnaire associé aux reads
        #  lecture de chaque position associée aux k-mers
        pos_r = 0  # position du k-mer sur le read. Le premier kmer commence toujours sur le premier nucléotide du read
        score = 0  # score de l'alignement entre le read et son ancrage sur le génome de référence
        for pos in kmer.values():
            if pos:  # si il y a bien une position associé au k-mer
                for position in pos:  # pour toute les positions du k-mer
                    read = key_dict_sens[i]  # read numéro i

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
                        dict_final_sens[key_dict_sens[i]] = position_finale
            pos_r += 1
        i += 1
    # Dictionnaire dict_final_sens {Read: Position} ne contenant que la meilleure position d'alignement pour chaque read du brin +.

    i = 0  # permet de déterminer la read que l'on aligne sur le génome
    for kmer in dict_kmer_position[1].values():  # lecture de chaque dictionnaire associé aux reads
        #  lecture de chaque position associée aux k-mers
        pos_r = 0  # position du k-mer sur le read. Le premier kmer commence toujours sur le premier nucléotide du read
        score = 0  # score de l'alignement entre le read et son ancrage sur le génome de référence
        for pos in kmer.values():
            if pos:  # si il y a bien une position associée au k-mer
                for position in pos:  # pour toute les positions du k-mer
                    read = key_dict_reverse[i]  # read numéro i
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
                        dict_final_reverse[key_dict_reverse[i]] = position_finale
            pos_r += 1
        i += 1
        # Dictionnaire dict_final_reverse {Read: Position} ne contenant que la meilleure position d'alignement pour chaque read du brin -.

  
    # CREATION DE LA TABLE VCF
    mat = [[], [], [], []] # initialisation de la table vcf
    tab_vcf = fill_vcf(mat, dict_final_sens, sequence_initiale) # Remplissage de la table vcf avec les reads brin + alignés
    tab_vcf = fill_vcf(tab_vcf, dict_final_reverse, sequence_initiale) # Remplissage de la matrice avec les reads brin - alignés
    tab_vcf = order_vcf(tab_vcf) # Ragement de la matrice vcf
    
    # ECRITURE DU FICHIER VCF
    with open(out_file, 'w') as vcf:
        # Ecriture des 1eres lignes du fichier :
        vcf.write("#REF: " + ref + "\n""#READS: " + reads +
                  "\n"'#K: ' + str(k) + '\n''#MAX_SUBST: ' +
                  str(max_hamming) + '\n''#MIN_ABUNDANCE: ' + str(min_abundance) + '\n')
        
        # Ecriture des données de table vcf en fonction de l'abondance minimum retenu
        i = 0
        while i < len(tab_vcf[0]):
            if tab_vcf[3][i] >= min_abundance:
                vcf.write(
                    str(tab_vcf[0][i]) + '\t' + tab_vcf[1][i] + '\t' + tab_vcf[2][i] + '\t' + str(
                        tab_vcf[3][i]) + '\n')
            i += 1



# mapping('smallMappingTest/reference.fasta', 'dumped_index.dp', 'smallMappingTest/reads.fasta', 20, 5, 2, 'snp15.vcf')


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
            print('python map.py --ref[genome_file.fa] --index[dumped_index.dp] --reads[reads.fa] -k[k_value] --max_hamming[h_value] --min_abundance[m_value] --out snps.vcf')
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
# python map.py --ref smallMappingTest/reference.fasta --index smallMappingTest/dumped_index_small.dp --reads smallMappingTest/reads.fasta -k 19 --max_hamming 5 --min_abundance 1 --out smallMappingTest/snps_ref_k19_d5.vcf


# python map.py --ref coli/ecoli_sample.fasta --index coli/dumped_index_coli.dp --reads coli/ecoli_mutated_reads_1000.fasta -k 20 --max_hamming 10 --min_abundance 1 --out coli/snps_coli_k20_d10_ab1.vcf

