# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, '../inphinity')
import pandas as pd
import numpy as np
import os

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

from similarity_matrix import constants
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.ProteinJ import ProteinJson
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()
Entrez.email = "christophe.joyet@heig-vd.ch"

def tests():
    """
    Test alignment algorithmes
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

    print('Longeur moyenne d\'une protéine : ' + str(length/len(list_prot)))

    # tests with local alignment
    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.localxx(ProteinA_sequence_AA, ProteinB_sequence_AA)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print('PaireWise2 - Smith–Waterman : Aln local + paramètres auto : temps moyen de comparaison de 2 protéines = ' + str(timeTot / len(list_prot)))

    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.localms(ProteinA_sequence_AA, ProteinB_sequence_AA, 1, 0, -1, -0.5)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print('PaireWise2 - Smith–Waterman : Aln local + paramètres arbitraires : temps moyen de comparaison de 2 protéines = ' + str(timeTot / len(list_prot)))

    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner1 = pairwise2.align.localds(ProteinA_sequence_AA, ProteinB_sequence_AA, blosum62, -1, -0.5)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print('PaireWise2 - Smith–Waterman : Aln local + matrice similarité BLOSUM62 : temps moyen de comparaison de 2 protéines = ' + str(timeTot / len(list_prot)))

    # tests with globals alignment
    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.globalxx(ProteinA_sequence_AA, ProteinB_sequence_AA)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print('PaireWise2 - Needleman–Wunsch : Aln global + paramètres auto : temps moyen de comparaison de 2 protéines = ' + str(timeTot / len(list_prot)))

    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.globalms(ProteinA_sequence_AA, ProteinB_sequence_AA, 1, 0, -1, -0.5)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print('PaireWise2 - Needleman–Wunsch : Aln global + paramètres arbitraires : temps moyen de comparaison de 2 protéines = ' + str(timeTot / len(list_prot)))

    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.globalds(ProteinA_sequence_AA, ProteinB_sequence_AA, blosum62, -1, -0.5)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print('PaireWise2 - Needleman–Wunsch : Aln global + matrice similarité BLOSUM62 : temps moyen de comparaison de 2 protéines = ' + str(timeTot / len(list_prot)))

    # Tests with optimised Smith-Waterman 
    query = StripedSmithWaterman(ProteinA_sequence_AA_skbio, protein=True, substitution_matrix=constants.MATRIX_BLOSUM62, gap_open_penalty=1, gap_extend_penalty=0.5)

    timeTot = 0

    for i in list_prot:
        timeStart = time.time()
        alignment = query(i.sequence_AA)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print('skbio.alignment - StripedSmithWaterman : Aln local + matrice similarité BLOSUM62 : temps moyen de comparaison de 2 protéines = ' + str(timeTot / len(list_prot)))

tests()