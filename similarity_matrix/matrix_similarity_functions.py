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

import general_functions
from similarity_matrix import constants
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.CoupleJ import CoupleJson
from objects_API.ProteinJ import ProteinJson
from objects_API.ContigJ import ContigJson
from objects_API.SpecieJ import SpecieJson
from objects_API.StrainJ import StrainJson
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()
Entrez.email = "christophe.joyet@heig-vd.ch"

def matrixSimilarityScript(file_name:str, path:str, organism_id:int, is_phage=False):
    """
    | Create a matrix of similarity with many organisms.
    | In this version, get all the couples with organism id (according with a 
    | specific level of lysis). 
    | Calculates the similarity between all the different organisms in the couples

    :Remark: To compare other organisms, modifie the content of list_organism_to_compare

    :param file_name: name of the similarity matrix
    :param path: where save the matrix
    :param organism_id: id of the organism to find couples
    :param is_phage: to differenciate bacterium and phages

    :type file_name: str
    :type path: str
    :type organism_id: int
    :type is_phage: boolean
    """

    lysis_type = constants.ALL_CLEAR_LYSIS 
    file_name = file_name
    path = path
    organism_dict = {}  # Get couple from specific bacterie
    
    if is_phage == False:
        organism_dict['bacterium'] = organism_id
    else:
        organism_dict['bacteriophage'] = organism_id

    # Get all the couples with the organism gave in parameter
    liste_couple = (CoupleJson.getCouplesByFilterParameter(organism_dict))

    # Select couple according with the lysis attribute
    liste_couple_final = []
    for couple in liste_couple:
        if couple.lysis in lysis_type:
            liste_couple_final.append(couple)

    list_organism_to_compare = []

    for couple in liste_couple_final:
        if is_phage == False:
            # check if there is no duplicate phage in the list
            if couple.bacteriophage not in list_organism_to_compare:
                list_organism_to_compare.append(couple.bacteriophage)
        else:
            # check if there is no duplicate phage in the list
            if couple.bacterium not in list_organism_to_compare:
                list_organism_to_compare.append(couple.bacterium)

    # Change the contents of list_organism_to_compare to compare other organisms
    # list_organism_to_compare = general_functions.getCouplesLysis(lysis_type)

    createSimilarityMatrix(list_organism_to_compare, file_name, path, not is_phage)

def getSimilarityScoreTwoProteinLocalAlign(proteinA:ProteinJson, proteinB:ProteinJson):
    """
    | Get similarity score between two proteins with local alignment

    :param proteinA: protein A
    :param proteinB: protein B to compare to protein A

    :type proteinA: ProteinJson
    :type phage_B: ProteinJson

    :return: similarity score between two phage
    :rtype: float16   
    """

    # using optimized Smit-Waterman
    query = StripedSmithWaterman(proteinA.sequence_AA, protein=True, substitution_matrix=constants.MATRIX_BLOSUM62, gap_open_penalty=10, gap_extend_penalty=1, score_only=True)
    queryscoremax = StripedSmithWaterman(proteinA.sequence_AA, protein=True, substitution_matrix=constants.MATRIX_BLOSUM62, gap_open_penalty=10, gap_extend_penalty=1, score_only=True)
    aligner_score = query(proteinB.sequence_AA)
    aligner_score_max = queryscoremax(proteinA.sequence_AA)

    # get scores from the return of query functions
    score = int(str(aligner_score).split(" ")[1])
    score_max = int(str(aligner_score_max).split(" ")[1])
    return score / score_max

def getSimilarityScoreTwoProteinLocalAlignText(proteinA:str, proteinB:str):
    """
    | Get similarity score between two proteins with local alignment.

    :param proteinA: protein A
    :param proteinB: protein B to compare to protein A

    :type proteinA: ProteinJson
    :type phage_B: ProteinJson

    :return: similarity score between two phage
    :rtype: float16   
    """

    # using optimized Smit-Waterman
    query = StripedSmithWaterman(proteinA, protein=True, substitution_matrix=constants.MATRIX_BLOSUM62, gap_open_penalty=10, gap_extend_penalty=1, score_only=True)
    queryscoremax = StripedSmithWaterman(proteinA, protein=True, substitution_matrix=constants.MATRIX_BLOSUM62, gap_open_penalty=10, gap_extend_penalty=1, score_only=True)
    aligner_score = query(proteinB)
    aligner_score_max = queryscoremax(proteinA)

    # get scores from the return of query functions
    score = int(str(aligner_score).split(" ")[1])
    score_max = int(str(aligner_score_max).split(" ")[1])
    return score / score_max

def getSimilarityScoreTwoProteinGlobalAlign(proteinA:ProteinJson, proteinB:ProteinJson):
    """
    | Get similarity score between two proteins with global alignment

    :param proteinA: protein A
    :param proteinB: protein B to compare to protein A

    :type proteinA: ProteinJson
    :type proteinB: ProteinJson

    :return: similarity score between two phage
    :rtype: float16   
    """
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
    # aligner = pairwise2.align.globalds(ProteinA_sequence_AA, ProteinB_sequence_AA, blosum62, -1, -0.5)
    # aligner = pairwise2.align.globalxx(ProteinA_sequence_AA, ProteinB_sequence_AA)
    aligner = pairwise2.align.globalms(ProteinA_sequence_AA, ProteinB_sequence_AA, 1, 0, -1, -0.5)
    
    # get the score and the length of sequence
    for a in aligner:
        al1, al2, score, begin, end = a

    # get the score
    return score / len(al1)

def getSimilarityScoreTwoPhages(phage_A:int, phage_B:int, is_local = True):
    """
    | Get similarity score between two phages

    :param phage_A: id_phage A
    :param phage_B: id_phage to compare to phage A

    :type phage_A: int
    :type phage_B: int

    :return: similarity score between two phage
    :rtype: float16   
    """
    protein_list_phage_1 = ProteinJson.getByOrganismID(phage_A)
    protein_list_phage_2 = ProteinJson.getByOrganismID(phage_B)

    # Contains best match between protein
    list_score_best_match_protein = []
    i = 0
    for protein_phage_1 in protein_list_phage_1:
        # Contains score between protein
        list_score_protein = []
        for protein_phage_2 in protein_list_phage_2:
            if is_local:
                list_score_protein.append(getSimilarityScoreTwoProteinLocalAlign(protein_phage_1, protein_phage_2))
            else:
                list_score_protein.append(getSimilarityScoreTwoProteinGlobalAlign(protein_phage_1, protein_phage_2))
        # Vector contained the best scores between phage A and phage B
        list_score_best_match_protein.append(np.max(list_score_protein))
        i += 1

    # Similarity based on the mean of all best match
    return np.mean(list_score_best_match_protein)

def createSimilarityMatrix(list_organism_to_compare:list, file_name:str, path:str, is_phage=True):
    """
    | Create a similarity matrix with a list of organisms

    :param list_organism_to_compare: list of organisms
    :param file_name: name of file where data will be stored
    :param path: where the matrix will be stored

    :type phage: list
    :type file_name: str
    :type path: str

    """
    matrice_similarity = []
    similarity = []
    matrix_size = len(list_organism_to_compare)
    i = 0
    
    for i in range(matrix_size):
        similarity = []
        for j in range(matrix_size):
            # If this is the same phage, rating max -> 1.0
            if list_organism_to_compare[i] == list_organism_to_compare[j]:
                similarity.append(1.0)
                continue
            # Add similarity score in the list similarity
            similarity.append(getSimilarityScoreTwoPhages(list_organism_to_compare[i], list_organism_to_compare[j]))
        # Add the score in the matrix
        matrice_similarity.append(similarity)
        i += 1
        # print(str(i) + "/" + str(len(list_organism_to_compare)) + " organism compared\n")
     
    # Get organisms names
    organism_name = []
    if is_phage == True:
        for organism in list_organism_to_compare:
            organism_name.append(BacteriophageJson.getByID(organism).designation)
    else:
        for organism in list_organism_to_compare:
            strain = StrainJson.getByID(BacteriumJson.getByID(organism).strain)
            organism_name.append(SpecieJson.getByID(strain.specie).designation + " " + strain.designation)

    # Save file
    df1 = pd.DataFrame(data=matrice_similarity, columns=organism_name, index=organism_name)

    if is_phage == True:
        df1 = df1.rename_axis('Phages designation', axis='columns')
        df1.to_csv(os.path.join(path, file_name), index_label= 'Phages designation')
    else:
        df1 = df1.rename_axis('Bacterium designation', axis='columns')
        df1.to_csv(os.path.join(path, file_name), index_label= 'Bacterium designation')

    ending_message = "file " + file_name + " saved in " + path
    
    print(ending_message)

def createProteiqueAlignmentLocalSequenceInFileTxt(phage_A:int, phage_B:int, path:str, similarity_min=0.8, similarity_max=1.0):
    """
    | Compare protein sequences between two phages and create a txt file.
    | Sequences will be saved if the similarity score of proteins is 
    | between similarity_min and similarity_max.  

    :param phage_A: id_phage A
    :param phage_B: id_phage to compare to phage A
    :param similarity_min: similarity min
    :param similarity_max: similarity max

    :type phage_A: int
    :type phage_B: int
    :type similarity_min: float value between 0.0 and 1.0
    :type similarity_max: float value between 0.0 and 1.0
    """
    # For the ending message
    isSomethingCompared = 0

    # Get proteins from phages
    protein_list_phage_1 = ProteinJson.getByOrganismID(phage_A)
    protein_list_phage_2 = ProteinJson.getByOrganismID(phage_B)

    file_name = BacteriophageJson.getByID(phage_A).designation + " compare to " + BacteriophageJson.getByID(phage_B).designation
    
    for protein_phage_1 in protein_list_phage_1:
        for protein_phage_2 in protein_list_phage_2:
                # Caclul similarity score 
                similarity_score = getSimilarityScoreTwoProteinLocalAlign(protein_phage_1, protein_phage_2)
                if similarity_score >= similarity_min and similarity_score <= similarity_max:
                    isSomethingCompared = 1
                    query = StripedSmithWaterman(protein_phage_1.sequence_AA, protein=True, substitution_matrix=constants.MATRIX_BLOSUM62, gap_open_penalty=10, gap_extend_penalty=1, score_only=False)
                    alignement = query(protein_phage_2.sequence_AA)                  
                    # Save informations in the file txt  
                    # Open file in append mode (add txt at the end of file)
                    fichier = open(path + file_name, "a")
                    fichier.write("\n\n" + str(protein_phage_1.description) + " compare to " + str(protein_phage_2.description) + "\n" + str(alignement.aligned_query_sequence) + "\n" + str(alignement.aligned_target_sequence) + "\n" + "%.3f"%(similarity_score*100) + "%\n")
                    fichier.close()
                    
    if isSomethingCompared == 1:
        ending_message = 'file ' + file_name + ' saved in ' + path
        print(ending_message)
    else:
        print('no correspondance found')



