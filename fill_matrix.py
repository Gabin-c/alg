#!/usr/bin/env python3
# -*- coding: utf-8 -*-

class DynamicMatrix:
    """
    Stock une matrice |S|x|T| avec S et T deux séquences. Permet d'avoir le score d'alignement entre les 2 séquences.
    Ne prend en compte que les match et mismatch.
    """

    def __init__(self, S, T):
        """
        Initialisation de la matrice d'alignement entre S et T
        :param S: Sequence 1
        :param T: Sequence 2
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
        Remplit la matrice de manière heuristique dans la diagonal central (cas ou S et T sont de longueur égale)
        avec les paramètres match = 1, mismatch = 0.

        :return: le score de l'alignement.
        """
        for i in range(1, len(self.S) + 1):
            for j in range(max(1, i), min(len(self.T) + 1, i + 1)):
                diagonal = self.matrix[i - 1][j - 1] + self.score(self.S[i - 1], self.T[j - 1])

                self.matrix[i][j] = diagonal
        return self.matrix[len(self.S)][len(self.T)]

    def printMatrix(self):  # A SUPPRIMER JUSTE POUR VERIFIER
        """ prints the matrix
        """
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

