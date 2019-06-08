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

#==============================================================================
#==============================================================================

file_name = "similarity_Streptomyces_Venezuelae-ATCC_10712_gap_1_0.5.csv"
path = "../../similarite/"
start = time.time()
bacterium_dict = {}
bacterium_dict['bacterium'] = 126

list_couple_clear_lysis = CoupleJson.getCouplesByFilterParameter(bacterium_dict)
list_phages_to_compare = []

for couple in list_couple_clear_lysis:
    # vérifie qu'il n'y a pas déjà le même phage dans la liste
    if not couple.bacteriophage in list_phages_to_compare:
        list_phages_to_compare.append(couple.bacteriophage)

ms.getSimilarityMatrix(list_phages_to_compare, file_name, path)
end = time.time()
print('%.2f' % (end - start))
#==============================================================================
#==============================================================================