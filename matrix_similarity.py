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

from Bio import Entrez
from Bio import Align 
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()
Entrez.email = "christophe.joyet@heig-vd.ch"

#==============================================================================
#==============================================================================

def getSimilarityScoreTwoProtein(proteinA:ProteinJson, proteinB:ProteinJson):
    """
    get similarity score between two proteins

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
    # 1 point is deducted for each gaps.
    aligner = pairwise2.align.globalms(ProteinA_sequence_AA, ProteinB_sequence_AA, 1, 0, -1, -1)
    
    # get the score and the length of sequence
    for a in aligner:
        al1, al2, score, begin, end = a

    # length al1 == al2
    # get the score
    return score / len(al1)

    '''
    print("Length PA : " + str(len(ProteinA_sequence_AA)))
    print("Length PB : " + str(len(ProteinB_sequence_AA)))
    print("Length al1 : " + str(len(al1)))
    print("Length al2 : " + str(len(al2)))
    '''
#==============================================================================
#==============================================================================

def getSimilarityScoreTwoPhages(phage_A:BacteriophageJson, phage_B:BacteriophageJson):
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

    '''
    print("length list 1 " + str(len(protein_list_phage_1)))
    print("length list 2 " + str(len(protein_list_phage_2)))'''

    # contient les meilleurs match entre les paires de protéines
    list_score_best_match_protein = []

    for protein_phage_1 in protein_list_phage_1:
        # contient les scores des paires de protéines
        list_score_protein = []
        for protein_phage_2 in protein_list_phage_2:
                list_score_protein.append(getSimilarityScoreTwoProtein(protein_phage_1, protein_phage_2))
        list_score_best_match_protein.append(np.max(list_score_protein))

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

    # devrait nous donner une matrice symétrique !
    matrix_size = len(list_phages_to_compare)
    for i in range(matrix_size):
        similarity = []
        for j in range(matrix_size):
            similarity.append(getSimilarityScoreTwoPhages(list_phages_to_compare[i], list_phages_to_compare[j]))
        matrice_similarity.append(similarity)

    #get phages name
    phages_name = []
    for phage in list_phages_to_compare:
        phages_name.append(BacteriophageJson.getByID(phage).designation)

    df1 = pd.DataFrame(data=matrice_similarity, columns=phages_name, index=phages_name)
    df1.to_csv(os.path.join(path, file_name))
    ending_message = "file " + file_name + " saved in " + path
    
    print(ending_message)



