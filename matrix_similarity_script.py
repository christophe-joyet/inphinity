#!/usr/bin/env python3
# File : matrix_similarity_script.py
# Author : Christophe Joyet
# Date : Mai 2019

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.OrganismJ import OrganismJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.ProteinJ import ProteinJson

import network_chart_couples_cl as network
import matrix_similarity as ms

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

import time

#=============================================================================================
#Constantes declarations
#=============================================================================================

#Lysis Type
CLEAR_LYSIS = 5
SEMI_CLEAR_LYSIS = 6
OPAQUE_LYSIS = 7
#dilution > 1e7
CLEAR_LYSIS_1E7PLUS = 8
SEMI_CLEAR_LYSIS_1E7PLUS = 10
#dilution < 1e7
CLEAR_LYSIS_1E7MINUS = 9
SEMI_CLEAR_LYSIS_1E7MINUS = 11

ALL_CLEAR_LYSIS = [CLEAR_LYSIS, CLEAR_LYSIS_1E7PLUS, CLEAR_LYSIS_1E7MINUS]
ALL_SEMI_CLEAR_LYSIS = [SEMI_CLEAR_LYSIS, SEMI_CLEAR_LYSIS_1E7PLUS, SEMI_CLEAR_LYSIS_1E7MINUS]

#=============================================================================================
#=============================================================================================

# choose what type of lysis we want
lysis_type = ALL_CLEAR_LYSIS

list_couples_lysis_type = []
list_couples_lysis_type = network.getCouplesLysis(lysis_type)

#=============================================================================================
#=============================================================================================
file_name = "similarity_Pseudomonas_aeruginosa_Muco16.csv"
path = "../../similarite/"

# bacterie to research by ID
bacterium_dict = {}
bacterium_dict['bacterium'] = 5190
liste_couple = (CoupleJson.getCouplesByFilterParameter(bacterium_dict))

# select only couples of a certain certain type
for couple in liste_couple:
    if not couple in list_couples_lysis_type:
        liste_couple.remove(couple)

list_phages_to_compare = []

for couple in liste_couple:
    # check if there is no duplicate phage in the list
    if not couple.bacteriophage in list_phages_to_compare:
        list_phages_to_compare.append(couple.bacteriophage)

ms.getSimilarityMatrix(list_phages_to_compare, file_name, path)
#=============================================================================================
#=============================================================================================