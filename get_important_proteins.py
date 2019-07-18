#!/usr/bin/env python
# encoding: utf-8

import os 

import constants
import matrix_similarity as ms
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from objects_API.ProteinJ import ProteinJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def getbestprot(organism_id_list:list, similarity_min=0.90, similarity_max=1.0):
    """
    Get the proteins who create the similarity inside a groupe of organisms

    :param organism_id_list: list of organisms' id
    :param similarity_min: similarity min between two proteins
    :param similarity_max: similarity max between two proteins

    :param organism_id_list: list of int
    :type similarity_min: float value between 0.0 and 1.0
    :type similarity_max: float value between 0.0 and 1.0
    
    """

    protein_list = []
    list_of_list_of_couples_grouped_by_organisms = []
    number_of_organisms = len(organism_id_list)
    # Get all the proteins of the organisms
    for organism in organism_id_list:
        protein_list.append(ProteinJson.getByOrganismID(organism))
    
    # This list will contain the paires of the proteins relative to organism
    # The first list represent the paires formed by protein of the first 
    # organism and the second organism. 
    # The second list represent the paires formed by protein of the first
    # organism and the third organism...
    list_of_couples_grouped_by_organisms = [[] for i in range (number_of_organisms - 1)]

    # Get the proteins of the first organism and make couples
    for protein_of_first_organism in protein_list[0]:
        for i in range (1, number_of_organisms):  # Begin at 1 to not compare the protein of first organism between them
            for protein_phage_n in protein_list[i]:
                # Compare proteins to get similarity
                similarity_score = ms.getSimilarityScoreTwoProteinLocalAlign(protein_of_first_organism, protein_phage_n)
                # Check if similarity is include between the thresholds
                if similarity_score >= similarity_min and similarity_score <= similarity_max:
                    couple = [protein_of_first_organism.sequence_AA, protein_phage_n.sequence_AA]  # Creation of couple
                    list_of_couples_grouped_by_organisms[i-1].append(couple) # Add the couple in the correct
    
    list_of_list_of_couples_grouped_by_organisms.append(list_of_couples_grouped_by_organisms)

    # Create the chain of similarity between proteins
    for i in range (number_of_organisms - 2):
        list_of_list_of_couples_grouped_by_organisms.append(getSimilarityProteinChain(list_of_list_of_couples_grouped_by_organisms[i], number_of_organisms, similarity_min, similarity_max))
    

def getSimilarityProteinChain(list_of_couples_grouped_by_organisms:list, number_of_organisms:int, similarity_min=0.70, similarity_max=1.0):
    """
    Get a list with couples of n proteique sequences resulting from
    comparison between the proteins from list_of_couples_grouped_by_organisms.
    (n = number_of_organisms - list_of_couples_grouped_by_organisms)


    :param list_of_couples_grouped_by_organisms: list of list of couples
    :param similarity_min: similarity min between two proteins
    :param similarity_max: similarity max between two proteins

    :param list_of_couples_grouped_by_organisms: list
    :type similarity_min: float value between 0.0 and 1.0
    :type similarity_max: float value between 0.0 and 1.0
    
    """
    number_of_organisms -= 1
    final_list_of_list_of_couples = [[] for i in range (len(list_of_couples_grouped_by_organisms) - 1)]

    for couple in list_of_couples_grouped_by_organisms[0]:
        for i in range (1, len(list_of_couples_grouped_by_organisms)):
            for j in range (len(list_of_couples_grouped_by_organisms[i])):
                    if couple[number_of_organisms - len(list_of_couples_grouped_by_organisms)] == list_of_couples_grouped_by_organisms[i][j][number_of_organisms - len(list_of_couples_grouped_by_organisms)]: # trouver les séquences du phages_1 identiques dans les deux couples des listes
                        similarity_score = ms.getSimilarityScoreTwoProteinLocalAlignText(couple[number_of_organisms - len(list_of_couples_grouped_by_organisms) + 1], list_of_couples_grouped_by_organisms[i][j][number_of_organisms - len(list_of_couples_grouped_by_organisms) + 1])
                        if similarity_score >= similarity_min and similarity_score <= similarity_max:
                            # creation d'un couple
                            couple2 = []
                            for k in range (number_of_organisms - len(list_of_couples_grouped_by_organisms) + 2):
                                couple2.append(couple[k])
                            couple2.append(list_of_couples_grouped_by_organisms[i][j][number_of_organisms - len(list_of_couples_grouped_by_organisms) + 1])
                            final_list_of_list_of_couples[i-1].append(couple2)

    return final_list_of_list_of_couples

getbestprot([5038,5030,5023,5045,5024])