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
        return self.matrix[len(self.S)][len(self.T)]

                


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
            
            

#dm = DynamicMatrix("CGAGCTGGTCCTAACCCGGAGACCGCAGGCTGCGCGCGTATCGCAGCATCTGGCATTACGCCGCATCGAGTGCATGCACGAGAGAAGGAAGGGCACTGTT", "CGAGCTGGTCCTAACCCGGAGACCGCAGGCTGCGCGCGTATCGCAGCATCTGGCATTACGCCGCATCGAGTGCATGCACGAGAGAAGGAAGGGCACTGTT", +1, 0, 0)
#dm.fillH(1)

#dm.printMatrix()
