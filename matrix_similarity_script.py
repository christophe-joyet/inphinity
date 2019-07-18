#!/usr/bin/env python
# encoding: utf-8

import constants
import network_chart_couples_cl as network
import matrix_similarity as ms
from objects_API.CoupleJ import CoupleJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.OrganismJ import OrganismJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.ProteinJ import ProteinJson
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def matrixSimilarityScript(file_name:str, path:str, organism_id:int, is_phage=False):
    """
    Create a matrix of similarity with many organisms
    In this version, we get all the couples we can make with the organism id (according with a 
    specific level of lysis). 
    We calculate the similarity between all the different organisms in the couples

    :param file_name: name of the similarity matrix
    :param path: where save the matrix
    :param organism_id: id of the organism to find couples
    :param is_phage: to differenciate bacterium and phages

    :type file_name: name of the similarity matrix
    :tpye path: where save the matrix
    :type organism_id: id of the organism to find couples
    :type is_phage: to differenciate bacterium and phages
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

    ms.getSimilarityMatrix(list_organism_to_compare, file_name, path, not is_phage)