#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 13:24:49 2020

@author: Gabin
"""

## Quel type d’objet python est utilisé pour stocker la matrice de programmation dynamique ? Avec
## quelle commande accède-t-on à la valeur de la cellule en face de la i-ème lettre de S et la j-ième de T ?

# Objet python : matrice.
# Commande pour accéder : self.matrix[i][j]


class DynamicMatrix:
    '''stores a matrix |S|x|T| (|S|+1 lines and |T|+1columns), sequences S and T and the score system (match, mismatch, gap)
        defines some global alignment functions
        '''
    def __init__(self, S, T, match, mismatch, gap):
        ''' defines and stores initial values'''
        
        self.S=S
        self.T=T
        self.gap=gap
        self.match=match
        self.mismatch=mismatch
        
        self.matrix = [0 for i in range(len(S)+1)]
        for i in range(len(S)+1):
            self.matrix[i] = [0 for j in range(len(T)+1)]


# Écrire une méthode score qui prend en argument 2 caractères et qui renvoie le score d’un match si
# les deux caractères sont égaux et le score d’un mismatch sinon
    def score(self,a, b):

        if a == b:
            return self.match
        else:
            return self.mismatch

        
# Écrire une méthode initGlobal qui initialise la matrice pour l’alignement global (première ligne et
# première colonne).
    def initGlobal(self):
        m, n = len(self.S), len(self.T)
    
        for i in range(0, m + 1):
            self.matrix[i][0] = self.gap * i
        
        for j in range(0, n + 1):
            self.matrix[0][j] = self.gap * j
        return self.matrix
        
    
# Écrire une méthode fill qui remplit la matrice selon l’algorithme de Needleman-Wunsch (formule
# de récurrence avec les trois cases voisines en haut et à gauche), et renvoie le score du meilleur alignement
# global des deux séquences.
    def fill(self):
        m, n = len(self.S), len(self.T)
        self.matrix = self.initGlobal()
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                diagonal = self.matrix[i - 1][j - 1] + self.score(self.S[i - 1], self.T[j - 1])
                up = self.matrix[i - 1][j] + self.gap
                left = self.matrix[i][j - 1] + self.gap
                max_pointer = max(diagonal, up, left)
                self.matrix[i][j] = max_pointer
        return self.matrix
    
# Proposer une méthode printGlobalAln qui affiche un alignement de meilleur score de S contre
# T comme dans l’exemple qui suit (dans cet exemple, le système de score suivant est utilisé : match=2,
# mismatch=-1, gap=-1) et qui renvoie son pourcentage d’identité.
    def printGlobalAln(self):
        align1, align2 = '', ''
        i, j = len(self.S), len(self.T)  # Dernière cellule de la matrice
        t = 0
    
        while i > 0 or j > 0:  
            score_current = self.matrix[i][j]
            score_diagonal = self.matrix[i - 1][j - 1]
            score_up = self.matrix[i][j - 1]
            score_left = self.matrix[i - 1][j]
            
            match_score_val = self.score(
                self.S[i - 1], self.T[j - 1])
    
            if score_current == score_diagonal + match_score_val:
                align1 += self.S[i - 1]
                align2 += self.T[j - 1]
                i -= 1
                j -= 1
                if match_score_val == self.match:
                    t += 1
            elif score_current == score_left + self.gap:
                align1 += self.S[i - 1]
                align2 += '-'
                i -= 1
            elif score_current == score_up + self.gap:
                align1 += '-'
                align2 += self.T[j - 1]
                j -= 1                                  
    
        print(align1[::-1])
        print(align2[::-1])
        print("Pourcentage d'identité : ")
        print(t/max(len(self.S),len(self.T))*100)

## Implémenter une heuristique de l’alignement global qui consiste à contraindre l’alignement autour
## de la diagonale dans une bande d’épaisseur 2 × width + 1, width sera un paramètre. On se contentera de
## calculer le score de l’alignement global (sans l’afficher).  
    def fillH(self, width):

        for i in range(1, len(self.S)+1):
            for j in range(max(1, i-width), min(len(self.T)+1,i+width+1)):

                diagonal = self.matrix[i-1][j-1] + self.score(self.S[i-1], self.T[j-1])
                up = self.matrix[i-1][j] + self.gap
                left = self.matrix[i][j-1] + self.gap

                self.matrix[i][j] = max(diagonal, up, left)

        for u in range(1,len(self.S))[::-1]:
            if self.matrix[u][len(self.T)] != 0:
                print(self.matrix[u][len(self.T)])
                break
                

        
    


    
    def printMatrix(self):
        ''' prints the matrix'''
        width = 4
        vide = " "
        line = f"{vide:>{2*width}}"
        for j in range(0,len(self.T)):
            line += f"{self.T[j]:>{width}}"
        print(line)
        line = f"{vide:>{width}}"
        for j in range(0,len(self.T)+1):
            line += f"{self.matrix[0][j]:>{width}}"
        print(line)
        for i in range(1,len(self.S)+1):
            line = f"{self.S[i-1]:>{width}}"
            for j in range(0,len(self.T)+1):
                line += f"{self.matrix[i][j]:>{width}}"
            print(line)
            
            

dm = DynamicMatrix("AATGAATCAAT", "GATAG", +2, -1, -1)
dm.printMatrix()
dm.initGlobal()
dm.fill()
dm.printMatrix()
dm.printGlobalAln()


dm.fillH(4)


import random
 


def adn_random(length):
    length =  length
    adn = ""
    i = 0
 
    random.seed()
 
    while 1 :
        while i < length :
            adn += random.choice('ATGC')
            i += 1
        i = 0
 
        if length % 2 != 0 or adn[0:int(length / 2)] != adn[int(length / 2):length]:
            break
 
    return adn



adn_501 = adn_random(500)
adn_502 = adn_random(500)
adn_1001 = adn_random(1000)
adn_1002 = adn_random(1000)
adn_2001 = adn_random(2000)
adn_2002 = adn_random(2000)

import time

tps1 = time.time()
dm = DynamicMatrix(adn_501,adn_502, +2, -1, -2)
dm.fill()
tps2 = time.time()
print(tps2 - tps1)

tps1 = time.time()
dm = DynamicMatrix(adn_501,adn_502, +2, -1, -2)
dm.fillH(1)
tps2 = time.time()
print(tps2 - tps1)

tps1 = time.time()
dm = DynamicMatrix(adn_501,adn_502, +2, -1, -2)
dm.fillH(5)
tps2 = time.time()
print(tps2 - tps1)

tps1 = time.time()
dm = DynamicMatrix(adn_501,adn_502, +2, -1, -2)
dm.fillH(10)
tps2 = time.time()
print(tps2 - tps1)

tps1 = time.time()
dm = DynamicMatrix(adn_501,adn_502, +2, -1, -2)
dm.fillH(20)
tps2 = time.time()
print(tps2 - tps1)

tps1 = time.time()
dm = DynamicMatrix(adn_501,adn_502, +2, -1, -2)
dm.fillH(40)
tps2 = time.time()
print(tps2 - tps1)


tps1 = time.time()
dm1000 = DynamicMatrix(adn_1001, adn_1002, +2, -1, -1)
dm1000.fill()
tps2 = time.time()
print(tps2 - tps1)


tps1 = time.time()
dm1000 = DynamicMatrix(adn_1001, adn_1002, +2, -1, -1)
dm1000.fillH(1)
tps2 = time.time()
print(tps2 - tps1)

tps1 = time.time()
dm1000 = DynamicMatrix(adn_1001, adn_1002, +2, -1, -1)
dm1000.fillH(5)
tps2 = time.time()
print(tps2 - tps1)

tps1 = time.time()
dm1000 = DynamicMatrix(adn_1001, adn_1002, +2, -1, -1)
dm1000.fillH(10)
tps2 = time.time()
print(tps2 - tps1)


tps1 = time.time()
dm1000 = DynamicMatrix(adn_1001, adn_1002, +2, -1, -1)
dm1000.fillH(20)
tps2 = time.time()
print(tps2 - tps1)

dm2000 = DynamicMatrix(adn_2001, adn_2002, +2, -1, -1)

tps1 = time.time()
dm2000.fill()
tps2 = time.time()
print(tps2 - tps1)


dm2000 = time.time()
dm2000 = DynamicMatrix(adn_2001, adn_2002, +2, -1, -1)
dm2000.fillH(1)
tps2 = time.time()
print(tps2 - tps1)

tps1 = time.time()
dm2000 = DynamicMatrix(adn_2001, adn_2002, +2, -1, -1)
dm2000.fillH(5)
tps2 = time.time()
print(tps2 - tps1)

tps1 = time.time()
dm2000 = DynamicMatrix(adn_2001, adn_2002, +2, -1, -1)
dm1000.fillH(10)
tps2 = time.time()
print(tps2 - tps1)


tps1 = time.time()
dm2000 = DynamicMatrix(adn_2001, adn_2002, +2, -1, -1)
dm2000.fillH(20)
tps2 = time.time()
print(tps2 - tps1)


tps1 = time.time()
dm2000 = DynamicMatrix(adn_2001, adn_2002, +2, -1, -1)
dm2000.fillH(40)
tps2 = time.time()
print(tps2 - tps1)


# On remarque que le temps diminue avec l'heuristique
# Plus la diagonale est large, plus ça se rapproche du score global
