#!/usr/bin/env python3
# File : stats_phages_script.py
# Author : Christophe Joyet
# Date : Juin 2019

import csv
import os
import pandas as pd
import network_chart_couples_cl as network
import stats_phages as sp


from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from collections import OrderedDict
from objects_API.CoupleJ import CoupleJson
from objects_API.ProteinJ import ProteinJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.OrganismJ import OrganismJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

# =====================================================================================================
# =====================================================================================================

repository = '../../statistiques/CSV/'
file_name = 'Features_from_phages_in_Pseudomonas_aeruginosa_Muco16.csv'

# =====================================================================================================
# =====================================================================================================

CLEAR_LYSIS = 5
CLEAR_LYSIS_P1E7 = 8 
CLEAR_LYSIS_M1E7 = 9

# get the phages of all the couples with an interaction type
# list_couple = network.getCouplesInteraction(CoupleJson.getAllAPI(), interaction_type=True)

# get the phages of specific couples 
list_couples_lysis_type = network.getCouplesLysis([CLEAR_LYSIS, CLEAR_LYSIS_P1E7, CLEAR_LYSIS_M1E7])

# get couple from specific bacterie
bacterium_dict = {}
# research bact by ID
bacterium_dict['bacterium'] = 5190
list_couple = (CoupleJson.getCouplesByFilterParameter(bacterium_dict))

# select only couples of a certain certain type
for couple in list_couple:
    if not couple in list_couples_lysis_type:
        list_couple.remove(couple)

# =====================================================================================================
# =====================================================================================================

# phages to extract features will be here
list_phages = []

for couple in list_couple:
    if not couple.bacteriophage in list_phages:
        list_phages.append(couple.bacteriophage)

print("nombre de phages différents : " + str(len(list_phages)))

phages_amino_acids = OrderedDict()

for couple in list_couple:
    if not couple.bacteriophage in phages_amino_acids.keys():
        phage_amino_acid_dict = sp.getFeaturesForAPhage(BacteriophageJson.getByID(couple.bacteriophage), active_percentage=False, get_features=True)
        phages_amino_acids[couple.bacteriophage] = phage_amino_acid_dict

# Source : 
# ADAM SMITH [Pseudonyme], 2015. This is beyond simple if you can use pandas. 
# Stackoverflow [en ligne]. 16 Juillet 2015 à 16h38. 
# [Consulté le 3 Avril 2019]. Disponible à l'adresse : https://stackoverflow.com/questions/31436783/writing-dictionary-of-dictionaries-to-csv-file-in-a-particular-format
df1 = pd.DataFrame.from_dict(data=phages_amino_acids, orient="index")
df1.to_csv(os.path.join(repository, file_name))

ending_message = "file " + file_name + " saved in " + repository
print(ending_message)

# =====================================================================================================
# =====================================================================================================
