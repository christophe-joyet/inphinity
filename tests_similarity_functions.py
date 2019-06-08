#!/usr/bin/env python3
# File : tests_similarity_functions.py
# Author : Christophe Joyet
# Date : Juin 2019

import pandas as pd
import numpy as np
import os

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.ProteinJ import ProteinJson
import network_chart_couples_cl as network
import time

from Bio import Entrez
from Bio import Align 
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat.MatrixInfo import blosum62

from skbio.alignment import StripedSmithWaterman 
from skbio.alignment import AlignmentStructure

from Bio.SubsMat import MatrixInfo

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()
Entrez.email = "christophe.joyet@heig-vd.ch"


#==============================================================================
#2D dict for StripedSmithWaterman
#==============================================================================

striped_mx = {'A' : {'A':4, 'R':-1, 'N':-2, 'D':-2, 'C':0, 'Q':-1, 'E':-1, 'G':0, 'H':-2, 'I':-1, 'L':-1, 'K':-1, 'M':-1, 'F':-2, 'P':-1, 'S':1, 'T':0, 'W':-3, 'Y':-2, 'V':0, 'B':-2, 'Z':-1, 'X':0, '*':-4},
              'R' : {'A':-1, 'R':5, 'N':0, 'D':-2, 'C':-3, 'Q':1, 'E':0, 'G':-2, 'H':0, 'I':-3, 'L':-2, 'K':2, 'M':-1, 'F':-3, 'P':-2, 'S':-1, 'T':-1, 'W':-3, 'Y':-2, 'V':-3, 'B':-1, 'Z':0, 'X':-1, '*':-4},
              'N' : {'A':-2, 'R':0, 'N':6, 'D':1, 'C':-3, 'Q':0, 'E':0, 'G':0, 'H':1, 'I':-3, 'L':-3, 'K':0, 'M':-2, 'F':-3, 'P':-2, 'S':1, 'T':0, 'W':-4, 'Y':-2, 'V':-3, 'B':3, 'Z':0, 'X':-1, '*':-4},
              'D' : {'A':-2, 'R':-2, 'N':1, 'D':6, 'C':-3, 'Q':0, 'E':2, 'G':-1, 'H':-1, 'I':-3, 'L':-4, 'K':-1, 'M':-3, 'F':-3, 'P':-1, 'S':0, 'T':-1, 'W':-4, 'Y':-3, 'V':-3, 'B':4, 'Z':1, 'X':-1, '*':-4},
              'C' : {'A':0, 'R':-3, 'N':-3, 'D':-3, 'C':9, 'Q':-3, 'E':-4, 'G':-3, 'H':-3, 'I':-1, 'L':-1, 'K':-3, 'M':-1, 'F':-2, 'P':-3, 'S':-1, 'T':-1, 'W':-2, 'Y':-2, 'V':-1, 'B':-3, 'Z':-3, 'X':-2, '*':-4},
              'Q' : {'A':-1, 'R':1, 'N':0, 'D':0, 'C':-3, 'Q':5, 'E':2, 'G':-2, 'H':0, 'I':-3, 'L':-2, 'K':1, 'M':0, 'F':-3, 'P':-1, 'S':0, 'T':-1, 'W':-2, 'Y':-1, 'V':-2, 'B':0, 'Z':3, 'X':-1, '*':-4},
              'E' : {'A':-1, 'R':0, 'N':0, 'D':2, 'C':-4, 'Q':2, 'E':5, 'G':-2, 'H':0, 'I':-3, 'L':-3, 'K':1, 'M':-2, 'F':-3, 'P':-1, 'S':0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2, 'B':1, 'Z':4, 'X':-1, '*':-4},
              'G' : {'A':0, 'R':-2, 'N':0, 'D':-1, 'C':-3, 'Q':-2, 'E':-2, 'G':6, 'H':-2, 'I':-4, 'L':-4, 'K':-2, 'M':-3, 'F':-3, 'P':-2, 'S':0, 'T':-2, 'W':-2, 'Y':-3, 'V':-3, 'B':-1, 'Z':-2, 'X':-1, '*':-4},
              'H' : {'A':-2, 'R':0, 'N':1, 'D':-1, 'C':-3, 'Q':0, 'E':0, 'G':-2, 'H':8, 'I':-3, 'L':-3, 'K':-1, 'M':-2, 'F':-1, 'P':-2, 'S':-1, 'T':-2, 'W':-2, 'Y':2, 'V':-3, 'B':0, 'Z':0, 'X':-1, '*':-4},
              'I' : {'A':-1, 'R':-3, 'N':-3, 'D':-3, 'C':-1, 'Q':-3, 'E':-3, 'G':-4, 'H':-3, 'I':4, 'L':2, 'K':-3, 'M':1, 'F':0, 'P':-3, 'S':-2, 'T':-1, 'W':-3, 'Y':-1, 'V':3, 'B':-3, 'Z':-3, 'X':-1, '*':-4},
              'L' : {'A':-1, 'R':-2, 'N':-3, 'D':-4, 'C':-1, 'Q':-2, 'E':-3, 'G':-4, 'H':-3, 'I':2, 'L':4, 'K':-2, 'M':2, 'F':0, 'P':-3, 'S':-2, 'T':-1, 'W':-2, 'Y':-1, 'V':1, 'B':-4, 'Z':-3, 'X':-1, '*':-4},
              'K' : {'A':-1, 'R':2, 'N':0, 'D':-1, 'C':-3, 'Q':1, 'E':1, 'G':-2, 'H':-1, 'I':-3, 'L':-2, 'K':5, 'M':-1, 'F':-3, 'P':-1, 'S':0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2, 'B':0, 'Z':1, 'X':-1, '*':-4},
              'M' : {'A':-1, 'R':-1, 'N':-2, 'D':-3, 'C':-1, 'Q':0, 'E':-2, 'G':-3, 'H':-2, 'I':1, 'L':2, 'K':-1, 'M':5, 'F':0, 'P':-2, 'S':-1, 'T':-1, 'W':-1, 'Y':-1, 'V':1, 'B':-3, 'Z':-1, 'X':-1, '*':-4},
              'F' : {'A':-2, 'R':-3, 'N':-3, 'D':-3, 'C':-2, 'Q':-3, 'E':-3, 'G':-3, 'H':-1, 'I':0, 'L':0, 'K':-3, 'M':0, 'F':6, 'P':-4, 'S':-2, 'T':-2, 'W':1, 'Y':3, 'V':-1, 'B':-3, 'Z':-3, 'X':-1, '*':-4},
              'P' : {'A':-1, 'R':-2, 'N':-2, 'D':-1, 'C':-3, 'Q':-1, 'E':-1, 'G':-2, 'H':-2, 'I':-3, 'L':-3, 'K':-1, 'M':-2, 'F':-4, 'P':7, 'S':-1, 'T':-1, 'W':-4, 'Y':-3, 'V':-2, 'B':-2, 'Z':-1, 'X':-2, '*':-4},
              'S' : {'A':1, 'R':-1, 'N':1, 'D':0, 'C':-1, 'Q':0, 'E':0, 'G':0, 'H':-1, 'I':-2, 'L':-2, 'K':0, 'M':-1, 'F':-2, 'P':-1, 'S':4, 'T':1, 'W':-3, 'Y':-2, 'V':-2, 'B':0, 'Z':0, 'X':0, '*':-4},
              'T' : {'A':0, 'R':-1, 'N':0, 'D':-1, 'C':-1, 'Q':-1, 'E':-1, 'G':-2, 'H':-2, 'I':-1, 'L':-1, 'K':-1, 'M':-1, 'F':-2, 'P':-1, 'S':1, 'T':5, 'W':-2, 'Y':-2, 'V':0, 'B':-1, 'Z':-1, 'X':0, '*':-4},
              'W' : {'A':-3, 'R':-3, 'N':-4, 'D':-4, 'C':-2, 'Q':-2, 'E':-3, 'G':-2, 'H':-2, 'I':-3, 'L':-2, 'K':-3, 'M':-1, 'F':1, 'P':-4, 'S':-3, 'T':-2, 'W':11, 'Y':2, 'V':-3, 'B':-4, 'Z':-3, 'X':-2, '*':-4},
              'Y' : {'A':-2, 'R':-2, 'N':-2, 'D':-3, 'C':-2, 'Q':-1, 'E':-2, 'G':-3, 'H':2, 'I':-1, 'L':-1, 'K':-2, 'M':-1, 'F':3, 'P':-3, 'S':-2, 'T':-2, 'W':2, 'Y':7, 'V':-1, 'B':-3, 'Z':-2, 'X':-1, '*':-4},
              'V' : {'A':0, 'R':-3, 'N':-3, 'D':-3, 'C':-1, 'Q':-2, 'E':-2, 'G':-3, 'H':-3, 'I':3, 'L':1, 'K':-2, 'M':1, 'F':-1, 'P':-2, 'S':-2, 'T':0, 'W':-3, 'Y':-1, 'V':4, 'B':-3, 'Z':-2, 'X':-1, '*':-4},
              'B' : {'A':-2, 'R':-1, 'N':3, 'D':4, 'C':-3, 'Q':0, 'E':1, 'G':-1, 'H':0, 'I':-3, 'L':-4, 'K':0, 'M':-3, 'F':-3, 'P':-2, 'S':0, 'T':-1, 'W':-4, 'Y':-3, 'V':-3, 'B':4, 'Z':1, 'X':-1, '*':-4},
              'Z' : {'A':-1, 'R':0, 'N':0, 'D':1, 'C':-3, 'Q':3, 'E':4, 'G':-2, 'H':0, 'I':-3, 'L':-3, 'K':1, 'M':-1, 'F':-3, 'P':-1, 'S':0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2, 'B':1, 'Z':4, 'X':-1, '*':-4},
              'X' : {'A':0, 'R':-1, 'N':-1, 'D':-1, 'C':-2, 'Q':-1, 'E':-1, 'G':-1, 'H':-1, 'I':-1, 'L':-1, 'K':-1, 'M':-1, 'F':-1, 'P':-2, 'S':0, 'T':0, 'W':-2, 'Y':-1, 'V':-1, 'B':-1, 'Z':-1, 'X':-1, '*':-4},
              '*' : {'A':-4, 'R':-4, 'N':-4, 'D':-4, 'C':-4, 'Q':-4, 'E':-4, 'G':-4, 'H':-4, 'I':-4, 'L':-4, 'K':-4, 'M':-4, 'F':-4, 'P':-4, 'S':-4, 'T':-4, 'W':-4, 'Y':-4, 'V':-4, 'B':-4, 'Z':-4, 'X':-4, '*':1}}

#==============================================================================
#==============================================================================

def tests():
    """
    test some alignment algorithmes
    """
    list_prot = ProteinJson.getByOrganismID(430)
    print("Organisme testé : " + BacteriophageJson.getByID(430).designation)
    print("Nombres de protéines comparées : " + str(len(list_prot)))
    print("Penalty Gap Open = 1")
    print("Penalty Gap Extend = 0.5")
    

    ProteinA_sequence_AA_skbio = list_prot[0].sequence_AA
    # IUPAC provides some informations about the sequence AA
    ProteinA_sequence_AA = Seq(ProteinA_sequence_AA_skbio, IUPAC.protein)

    # display avg length of a protein sequence
    length = 0
    for prot in list_prot:
        length += len(prot.sequence_AA)

    print("Longeur moyenne d'une protéine : " + str(length/len(list_prot)))

    # tests with local alignment
    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.localxx(ProteinA_sequence_AA, ProteinB_sequence_AA)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print("PaireWise2 - Smith–Waterman : Aln local + paramètres auto : temps moyen de comparaison de 2 protéines = " + str(timeTot / len(list_prot)))

    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.localms(ProteinA_sequence_AA, ProteinB_sequence_AA, 1, 0, -1, -0.5)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print("PaireWise2 - Smith–Waterman : Aln local + paramètres arbitraires : temps moyen de comparaison de 2 protéines = " + str(timeTot / len(list_prot)))

    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner1 = pairwise2.align.localds(ProteinA_sequence_AA, ProteinB_sequence_AA, blosum62, -1, -0.5)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print("PaireWise2 - Smith–Waterman : Aln local + matrice similarité BLOSUM62 : temps moyen de comparaison de 2 protéines = " + str(timeTot / len(list_prot)))

    # tests with globals alignment
    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.globalxx(ProteinA_sequence_AA, ProteinB_sequence_AA)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print("PaireWise2 - Needleman–Wunsch : Aln global + paramètres auto : temps moyen de comparaison de 2 protéines = " + str(timeTot / len(list_prot)))

    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.globalms(ProteinA_sequence_AA, ProteinB_sequence_AA, 1, 0, -1, -0.5)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print("PaireWise2 - Needleman–Wunsch : Aln global + paramètres arbitraires : temps moyen de comparaison de 2 protéines = " + str(timeTot / len(list_prot)))

    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.globalds(ProteinA_sequence_AA, ProteinB_sequence_AA, blosum62, -1, -0.5)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print("PaireWise2 - Needleman–Wunsch : Aln global + matrice similarité BLOSUM62 : temps moyen de comparaison de 2 protéines = " + str(timeTot / len(list_prot)))

    # tests with optimised Smith-Waterman 
    query = StripedSmithWaterman(ProteinA_sequence_AA_skbio, protein=True, substitution_matrix=striped_mx, gap_open_penalty=1, gap_extend_penalty=0.5)

    timeTot = 0

    for i in list_prot:
        timeStart = time.time()
        alignment = query(i.sequence_AA)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print("skbio.alignment - StripedSmithWaterman : Aln local + matrice similarité BLOSUM62 : temps moyen de comparaison de 2 protéines = " + str(timeTot / len(list_prot)))

tests()