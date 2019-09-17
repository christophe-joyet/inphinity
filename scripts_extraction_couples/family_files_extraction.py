#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, './')

from objects_API.CoupleJ import CoupleJson
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

import matplotlib.pyplot as plt
import csv
import numpy as np

from general_graphs import constants
from objects_API.StrainJ import StrainJson
from objects_API.SpecieJ import SpecieJson
from objects_API.GenusJ import GenusJson
from objects_API.FamilyJ import FamilyJson
from objects_API.BacteriumJ import BacteriumJson


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def getCouplesOfGreg():
    list_couples_greg = CoupleJson.getCouplesByFilterParameter({'source_data_id':3})

    return list_couples_greg

list_phages = []
list_bacteria = []
list_couples_greg = getCouplesOfGreg()
list_couples_greg_id = []
list_couples_positive_interaction = []
list_couples_negative_interaction = []

for couple in list_couples_greg:
    if couple.bacteriophage not in list_phages:
        list_phages.append(couple.bacteriophage)
    if couple.bacterium not in list_bacteria:
        list_bacteria.append(couple.bacterium)
    if couple.interaction_type == True:
        list_couples_positive_interaction.append(couple.id)
    else:
        list_couples_negative_interaction.append(couple.id)
    
    list_couples_greg_id.append(str(couple.id))