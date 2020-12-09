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

    def __init__(self, S, T ):
        """
        Initialisation de la matrice d'alignement entre S et T
        :param S: Sequence
        :param T: Sequence
        """

        self.S = S
        self.T = T

        self.matrix = [0 for i in range(len(S) + 1)]
        for i in range(len(S) + 1):
            self.matrix[i] = [0 for j in range(len(T) + 1)]


    def score(self, a, b):
        """
        Comparaison de caractère, si le caractère sont égaux retourne 1 sinon retourne 0.

        :param a: Caractère
        :param b: Caractère
        :return: 1 ou 0
        """
        if a == b:
            return 1
        else:
            return 0

    def fillH(self):
        """
        Remplit la matrice de manière heuristique dans la diagonal central ( cas ou S et T sont de longueur
        égale ) avec les paramètres match = 1, mismatch = 0.

        :return: le score de l'alignement.
        """
        for i in range(1, len(self.S) + 1):
            for j in range(max(1, i), min(len(self.T) + 1, i + 1)):
                diagonal = self.matrix[i - 1][j - 1] + self.score(self.S[i - 1], self.T[j - 1])

                self.matrix[i][j] = diagonal
        return self.matrix[len(self.S)][len(self.T)]

    def printMatrix(self): # A SUPPRIMER JUSTE POUR VERIFIER
        ''' prints the matrix'''
        width = 4
        vide = " "
        line = f"{vide:>{2 * width}}"
        for j in range(0, len(self.T)):
            line += f"{self.T[j]:>{width}}"
        print(line)
        line = f"{vide:>{width}}"
        for j in range(0, len(self.T) + 1):
            line += f"{self.matrix[0][j]:>{width}}"
        print(line)
        for i in range(1, len(self.S) + 1):
            line = f"{self.S[i - 1]:>{width}}"
            for j in range(0, len(self.T) + 1):
                line += f"{self.matrix[i][j]:>{width}}"
            print(line)


mat = DynamicMatrix("ATCCCCCCTCCC","ATCCCCCCCCCC")
mat.printMatrix()
mat.fillH()
mat.printMatrix()