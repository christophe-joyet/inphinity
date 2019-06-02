#!/usr/bin/env python3
# File : matrix_similarity.py
# Author : Christophe Joyet
# Date : Mai 2019

import pandas as pd
import numpy as np
import os

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.OrganismJ import OrganismJson
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

from skbio.alignment import local_pairwise_align_ssw 
from skbio.alignment import StripedSmithWaterman 
from skbio.alignment import AlignmentStructure

from Bio.SubsMat import MatrixInfo
from skbio.sequence import Protein




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

def getSimilarityScoreTwoProteinLocalAlignStripedSmithWaterman(proteinA:ProteinJson, proteinB:ProteinJson):
    """
  
    """
    # =======================================================================
    # Align the two sequences
    # =======================================================================

    # IUPAC provides some informations about the sequence AA
    ProteinA_sequence_AA = proteinA.sequence_AA
    ProteinB_sequence_AA = proteinB.sequence_AA

    query = StripedSmithWaterman(ProteinA_sequence_AA)
    aligner = query(ProteinB_sequence_AA)

    
    # get the score and the length of sequence
    print(aligner)
    
   
    #return (score / length)

#==============================================================================
#==============================================================================

def getSimilarityScoreTwoProteinLocalAlign(proteinA:ProteinJson, proteinB:ProteinJson):
    """
    get similarity score between two proteins with local alignment

    :param proteinA: protein A
    :param proteinB: protein B to compare to protein A

    :type proteinA: ProteinJson
    :type phage_B: ProteinJson

    :return: similarity score between two phage
    :rtype: float16   
    """
    # =======================================================================
    # Align the two sequences
    # =======================================================================

    # IUPAC provides some informations about the sequence AA
    ProteinA_sequence_AA = Seq(proteinA.sequence_AA, IUPAC.protein)
    ProteinB_sequence_AA = Seq(proteinB.sequence_AA, IUPAC.protein)

    # Match parameters :
    #   x     No parameters. Identical characters have score of 1, otherwise 0.
    #   m     A match score is the score of identical chars, otherwise mismatch
    #         score.
    #   d     A dictionary returns the score of any pair of characters.
    #   c     A callback function returns scores.
        
    # GAP penalty parameters :
    #   x     No gap penalties.
    #   s     Same open and extend gap penalties for both sequences.
    #   d     The sequences have different open and extend gap penalties.
    #   c     A callback function returns the gap penalties.
    
    # Do a global alignment. 
    # 1 point is given for identical characters
    # 0 point is deducted for each non-identical character. 
    # 1 point is deducted for each open gaps.
    # 0.5 point is deducted for each extended gaps.
    # aligner = pairwise2.align.localds(ProteinA_sequence_AA, ProteinB_sequence_AA, blosum62, -1, -0.5)
    # aligner = pairwise2.align.localxx(ProteinA_sequence_AA, 1, 0)
    # aligner = pairwise2.align.localms(ProteinA_sequence_AA, ProteinB_sequence_AA, 1, 0, -1, -0.5)
    query = StripedSmithWaterman(proteinA.sequence_AA, protein=True, substitution_matrix=striped_mx, gap_open_penalty=1, gap_extend_penalty=0.5, score_only=True)
    queryscoremax = StripedSmithWaterman(proteinA.sequence_AA, protein=True, substitution_matrix=striped_mx, gap_open_penalty=1, gap_extend_penalty=0.5, score_only=True)
    aligner_score1 = query(proteinB.sequence_AA)
    aligner_score2 = queryscoremax(proteinA.sequence_AA)

    score1 = int(str(aligner_score1).split(" ")[1])
    score2 = int(str(aligner_score2).split(" ")[1])
    # get the score and the length of sequence
    '''for a in aligner:
        al1, al2, score, begin, end = a'''

    # length al1 == al2
    # get the score with ponderation WITH BLOSUM62
    '''
    perfect_score = pairwise2.align.localds(ProteinA_sequence_AA, ProteinA_sequence_AA, blosum62, -10, -1, score_only=True)
    return ((score / perfect_score) * ((end - begin) / len(al1)))
    '''
    # get the score with ponderation WITHOUT BLOSUM62
    #return ((score / len(al1)) * ((end - begin) / len(al1)))

    return score1 / score2

#==============================================================================
#==============================================================================

def getSimilarityScoreTwoProteinGlobalAlign(proteinA:ProteinJson, proteinB:ProteinJson):
    """
    get similarity score between two proteins with global alignment

    :param proteinA: protein A
    :param proteinB: protein B to compare to protein A

    :type proteinA: ProteinJson
    :type phage_B: ProteinJson

    :return: similarity score between two phage
    :rtype: float16   
    """
    # =======================================================================
    # Align the two sequences
    # =======================================================================

    # IUPAC provides some informations about the sequence AA
    ProteinA_sequence_AA = Seq(proteinA.sequence_AA, IUPAC.protein)
    ProteinB_sequence_AA = Seq(proteinB.sequence_AA, IUPAC.protein)

    # Match parameters :
    #   x     No parameters. Identical characters have score of 1, otherwise 0.
    #   m     A match score is the score of identical chars, otherwise mismatch
    #         score.
    #   d     A dictionary returns the score of any pair of characters.
    #   c     A callback function returns scores.
        
    # GAP penalty parameters :
    #   x     No gap penalties.
    #   s     Same open and extend gap penalties for both sequences.
    #   d     The sequences have different open and extend gap penalties.
    #   c     A callback function returns the gap penalties.
    
    # Do a global alignment. 
    # 1 point is given for identical characters
    # 0 point is deducted for each non-identical character. 
    # 1 point is deducted for each open gaps.
    # 0.5 point is deducted for each extended gaps.
    aligner = pairwise2.align.globalms(ProteinA_sequence_AA, ProteinB_sequence_AA, 1, 0, -1, -0.5)
    
    # get the score and the length of sequence
    for a in aligner:
        al1, al2, score, begin, end = a

    # length al1 == al2
    # get the score
    print("MOTHERFUCKER")
    return score / len(al1)

    
#==============================================================================
#==============================================================================

def getSimilarityScoreTwoPhages(phage_A:BacteriophageJson, phage_B:BacteriophageJson, is_local = True):
    """
    get similarity score between two phages

    :param phage_A: phage A
    :param phage_B: phage to compare to phage A

    :type phage_A: BacteriophageJson
    :type phage_B: BacteriophageJson

    :return: similarity score between two phage
    :rtype: float16   
    """

    protein_list_phage_1 = ProteinJson.getByOrganismID(phage_A)
    protein_list_phage_2 = ProteinJson.getByOrganismID(phage_B)

    
    print("length list 1 " + str(len(protein_list_phage_1)))
    print("length list 2 " + str(len(protein_list_phage_2)))

    # contient les meilleurs match entre les paires de protéines
    list_score_best_match_protein = []
    i = 0
    for protein_phage_1 in protein_list_phage_1:
        # contient les scores des paires de protéines
        list_score_protein = []
        for protein_phage_2 in protein_list_phage_2:
            if is_local:
                list_score_protein.append(getSimilarityScoreTwoProteinLocalAlign(protein_phage_1, protein_phage_2))
            else:
                list_score_protein.append(getSimilarityScoreTwoProteinGlobalAlign(protein_phage_1, protein_phage_2))
        list_score_best_match_protein.append(np.max(list_score_protein))
        i += 1
        print(str(i) + "/" + str(len(protein_list_phage_1)) + " protein compared\n")

    '''
    print("Score max d'un match = " + str(np.amax(list_score_best_match_protein)))
    print("Length of list best match = " + str(len(list_score_best_match_protein)))'''

    # similarité basé sur la moyenne des meilleurs matchs des protéines 
    return np.mean(list_score_best_match_protein)


#==============================================================================
#==============================================================================

def getSimilarityMatrix(list_phages_to_compare:list, file_name:str, path:str):
    """
    get similarity matrix between a list of phages

    :param list_phages_to_compare: list of phages
    :param file_name: name of file where data will be stored
    :param path: where the matrix will be stored

    :type phage: list
    :type file_name: str
    :type path: str

    """

    print("Number of phages to analysed : " + (str(len(list_phages_to_compare))))

    matrice_similarity = []
    similarity = []
    i = 0

    matrix_size = len(list_phages_to_compare)
    for i in range(matrix_size):
        similarity = []
        for j in range(matrix_size):
            similarity.append(getSimilarityScoreTwoPhages(list_phages_to_compare[i], list_phages_to_compare[j]))
        matrice_similarity.append(similarity)
        i += 1
        print(str(i) + "/" + str(len(list_phages_to_compare)) + " phages compared\n")


    #get phages name
    phages_name = []
    for phage in list_phages_to_compare:
        phages_name.append(BacteriophageJson.getByID(phage).designation)

    df1 = pd.DataFrame(data=matrice_similarity, columns=phages_name, index=phages_name)
    df1.to_csv(os.path.join(path, file_name))
    ending_message = "file " + file_name + " saved in " + path
    
    print(ending_message)

def petitmodele():
    # IUPAC provides some informations about the sequence AA
    ProteinA_sequence_AA = Seq("MVARRKGERVVRKNEVENVQQRACANRRERQRTKELNDAFTLLRKLIPSMPSDKMSKIHTLRIATDYISFLDEMQKNGCKLYGHSIFDEKRGYNLQSAFNMWRGNNGYTPIAGPSQLPPLQSAHIPPPAPSSIPPHCLMPQPWYQTCPPPKQEFHELCPISTPNPNSNPNQLTPIHWQ", IUPAC.protein)
    ProteinB_sequence_AA = Seq("MVARRKGERVVRKNEVENVQQRACANRRERQRTKELNDAFTLLRKLIPSMPSDKMSKIHTLRIATDYISFLDEMQKNGCKLYGHSIFDEKRGYNLQSAFNMWRGNNGYTPIAGPSQLPPLQSAHIPPPAPSSIPPHCLMPQPWYQTCPPPKQEFHELCPISTPNPNSNPNQLTPIHWQ", IUPAC.protein)
    list_prot = ProteinJson.getByOrganismID(2050)
   
    print("MS")
    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.localds(ProteinA_sequence_AA, ProteinB_sequence_AA, blosum62, -1, -0.5)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    print("avg = " + str(timeTot / len(list_prot)))

    timeTot = 0

    for i in list_prot:
        ProteinB_sequence_AA = Seq(i.sequence_AA, IUPAC.protein)
        timeStart = time.time()
        aligner2 = pairwise2.align.globalds(ProteinA_sequence_AA, ProteinB_sequence_AA, blosum62, -1, -0.5)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart
    
    # get the score and the length of sequence
    for a in aligner2:
        al1, al2, score, begin, end = a

    '''print("Len 1 : " + str(len(ProteinA_sequence_AA)))
    print("Len 2 : " + str(len(ProteinB_sequence_AA)))
    print("al1 = " + str(len(al1)))
    print("al2 = " + str(len(al2)))
    print("Score = " + str(score))
    print("begin = " + str(begin))
    print("end = " + str(end))'''
    print("avg = " + str(timeTot / len(list_prot)))

    query = StripedSmithWaterman("MVARRKGERVVRKNEVENVQQRACANRRERQRTKELNDAFTLLRKLIPSMPSDKMSKIHTLRIATDYISFLDEMQKNGCKLYGHSIFDEKRGYNLQSAFNMWRGNNGYTPIAGPSQLPPLQSAHIPPPAPSSIPPHCLMPQPWYQTCPPPKQEFHELCPISTPNPNSNPNQLTPIHWQ", protein=True, substitution_matrix=striped_mx)
    alignment = query(i.sequence_AA)

    timeTot = 0

    for i in list_prot:
        timeStart = time.time()
        alignment = query(i.sequence_AA)
        timeEnd = time.time()
        timeTot += timeEnd - timeStart

    for i in alignment:
        print(str(i))
    print("avg = " + str(timeTot / len(list_prot)))

    # length al1 == al2
    # get the score'''

#petitmodele()



        
