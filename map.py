#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pickle
import pandas as pd
import fill_matrix as matrix
import getopt
import sys

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


def get_occurrences(pattern: str, bwt: str, n: {}, r: [], sa: [int]) -> []:
    """
    Retourne les positions des occurences du pattern dans la séquence de référence à l'aide de bwt, sa , n , r.
    On obtient en sortie une liste des positions des occurences.
    
    :param pattern: séquencce
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


# get_occurrences("T", my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])


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


def get_kmer_position(k, reads):
    """
    Retourne un dictionnaire contenant les kmer et leurs positions possibles pour chaque read :
    { Read : { k-mer : [position] } }
    :param k: longueur du kmer
    :param reads: fichier fasta de reads
    :return: un dictionnaire { Read : { k-mer : [position] } }
    """
    with open(reads) as reads_file:
        kmer_position = {}  # dictionnaire vide
        for line in reads_file:
            if line[0] != ">":
                read_line = line.strip()  # retire les lignes qui ne sont pas des reads
                start = -1  # -1 car sinon ça commence à la deuxième lettre du read (je ne sais pas pourquoi)
                end = k - 1  # longueur du k-mer
                kmer_position[read_line] = {}  # dictionnaire avec les reads
                for kmers in range(start, len(read_line), k):  # recherche du k-mer sur le genome de reference
                    while end <= len(read_line)-1: 
                        start += 1
                        end += 1
                        occ = [] # INUTILE ???
                        kmer_position[read_line][read_line[start:end]] = occ  # dictionnaire kmer : position
                        if len(read_line[start:end]) == k:  # possible k-mer qui ne fait pas la taille k demandée
                            occ = get_occurrences(read_line[start:end], my_fmi[0], my_fmi[2], my_fmi[3], my_fmi[1])
                            kmer_position[read_line][read_line[start:end]] += occ
    return kmer_position


def mapping(ref, index, reads, k, max_hamming, min_abundance, out_file):
    with open(index, "rb") as f2:
        fmi = pickle.load(f2)
        
    # Création du dictionnaire { Read : { k-mer : [position] } } contenant
    # les positions des occurences k-mers pour chaque read
    dict_kmer_position = get_kmer_position(k, reads) 
    
    # Recherche des meilleurs positions d'alignements
    
    dict_final = {} #Création d'un dictionnaire qui ne contiendra que les meilleurs position d'alignement pour chaque read.
    key_dict = list(dict_kmer_position.keys()) # stockage des différentes reads à alignés
    sequence_initiale = bwt_2_seq(fmi[0], fmi[2], fmi[3]) # récupération de la séquence initiale
    
    i = 0 # permet de déterminer la read que l'on aligne sur le génome 
    for kmer in dict_kmer_position.values(): # lecture de chaque dictionnaire associé aux reads
        
        # lecture de chaque position associé aux k-mers
        pos_r = 0 # position du k-mer sur le read le premier kmer commence toujours sur le premier nucléotide du read
        score = 0 # score de l'alignement entre le read et son ancrage sur le génome de référence        
        for pos in kmer.values(): 
            
            if pos: # si il y a bien une position associé au k-mer 
                
                for position in pos: # pour toute les positions du k-mer
                    read = key_dict[i] # read numéro i
                    # Création d'une matrice d'alignement initialisé à 0 entre le read et son ancrage sur le génome
                    # la position du read sur le génome est déterminé par la position de l'alignement du k-mer sur le génome (position)
                    # et la position du k-mer sur le read (pos_r). le read s'alignera de la position "position - pos_r" à cette même position + la 
                    # longueur du read.
                    dm = matrix.DynamicMatrix(read,
                                              sequence_initiale[(position - pos_r):(position - pos_r + len(read))], 1,
                                              0, 0) #MODIFIER
                    
                    # Calcul du score d'alignement
                    fill_mat = dm.fillH(1) #MODIFIER
                    
                    # Vérification du score 
                    # si le score est plus grand que le précédent et respecte le nombre max de substitutions autorisé
                    # alors le score prédédent et la meilleur position d'alignement remplacé
                    if score < fill_mat and fill_mat >= (len(read) - max_hamming): 
                        position_finale = (position - pos_r)
                        score = fill_mat
                        dict_final[key_dict[i]] = position_finale 
            pos_r += 1
        i += 1
    # Optention d'un dictionnaire dict_final ( Read : Position ) ne contenant que la meilleur position d'alignement pour chaque read
        
    
    # Création et remplisage de la tab VCF
    
    ## création de la table VCF
    column = ["POS", "REF", "ALT", "ABUNDANCE"] 
    data_vcf = pd.DataFrame(columns=column)
    
    # Remplissage de la table VCF
    for cle, valeur in dict_final.items():
        pos_read = 0
        for read, refer in zip(cle, sequence_initiale[valeur:(valeur+100)]):
            if read != refer:
                if data_vcf.loc[data_vcf['POS'] == (valeur + pos_read)].empty is False:
                    if data_vcf.loc[data_vcf['POS'] == (valeur + pos_read)]['ALT'].empty is True:
                        new_line = {'POS': (valeur + pos_read), 'REF': refer, 'ALT': read, 'ABUNDANCE': 1}
                        data_vcf = data_vcf.append(new_line, ignore_index=True)
                    else:
                        if data_vcf['ABUNDANCE'][data_vcf['POS'] == (valeur + pos_read)][data_vcf['ALT'] == read].empty\
                                is False:
                            data_vcf['ABUNDANCE'][data_vcf['POS'] == (valeur + pos_read)] += 1
                else:
                    new_line = {'POS': (valeur + pos_read), 'REF': refer, 'ALT': read, 'ABUNDANCE': 1}
                    data_vcf = data_vcf.append(new_line, ignore_index=True)
            pos_read += 1

    with open(out_file, 'w') as vcf:
        vcf.write("#REF: " + ref + "\n""#READS: " + reads +
                  "\n"'#K: ' + str(k) + '\n''#MAX_SUBST: ' +
                  str(max_hamming) + '\n''#MIN_ABUNDANCE: ' + str(min_abundance) + '\n')
    data_vcf[data_vcf.ABUNDANCE >= min_abundance].to_csv(out_file, index=None, sep='\t', mode='a')

#  mapping('smallMappingTest/reference.fasta', 'dumped_index.dp', 'smallMappingTest/reads.fasta', 19, 5, 1, 'snps.vcf')


if __name__ == "__main__":
    reference = ''
    index_file = ''
    read_file = ''
    k_mers = 1
    hamming = 0
    abundance = 0
    out = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:], "k:h:", ["ref=", "index=", "reads=", "max_hamming=", "min_abundance=",
                                                          "out="])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)
    for option, arg in opts:
        if option in "-h":
            print('python map.py --ref[genome_file.fa] --index[dumped_index.dp] --reads[reads.fa] -k[k_value] '
                  '--max_hamming[h_value] --min_abundance[m_value] --out snps.vcf')
            sys.exit()
        elif option in "--ref":
            reference = arg
        elif option in "--index":
            index_file = arg
        elif option in "--reads":
            read_file = arg
        elif option in "-k":
            k_mers = arg
        elif option in "--max_hamming":
            hamming = arg
        elif option in "--min_abundance":
            abundance = arg
        elif option in "--ref":
            ref_file = arg
        elif option in "--out":
            out = arg

    map(reference, index_file, read_file, k_mers, hamming, abundance, out)


#  python map.py --ref smallMappingTest/reference.fasta --index dumped_index.dp --reads smallMappingTest/reads.fasta
#  -k 19 --max_hamming 5 --min_abundance 3 --out snps.vcf
