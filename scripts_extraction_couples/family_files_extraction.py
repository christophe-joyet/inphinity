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
    #list_couples_greg = CoupleJson.getCouplesByFilterParameter({'source_data_id':3})
    list_couples_positive_interaction = CoupleJson.getCouplesByFilterParameter({'interaction_type':True})

    return list_couples_positive_interaction

list_phages = []
list_bacteria = []
list_couples_greg = getCouplesOfGreg()
list_couples_greg_id = []

for couple in list_couples_greg:
    if couple.bacteriophage not in list_phages:
        list_phages.append(couple.bacteriophage)
    if couple.bacterium not in list_bacteria:
        list_bacteria.append(couple.bacterium) 
    list_couples_greg_id.append(str(couple.id))

flag = 0

with open('/home/christophe/Documents/Phages4A/CH_dataset.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')

    for row in csv_reader:
        #transform row in str
        row_str = ""
        for el in row:
            row_str += el + ','
        row_str += '\n'
        
        #first line
        if flag == 0:
            with open('/home/christophe/Documents/Phages4A/CH_dataset_staphylococcaceae_positive.csv', mode='a') as csv_staph:
                    csv_staph.write(row_str)
            with open('/home/christophe/Documents/Phages4A/CH_dataset_moraxellaceae_positive.csv', mode='a') as csv_morax:
                        csv_morax.write(row_str)
            with open('/home/christophe/Documents/Phages4A/CH_dataset_enterobacteriaceae_positive.csv', mode='a') as csv_entero:
                        csv_entero.write(row_str)
            with open('/home/christophe/Documents/Phages4A/CH_dataset_pseudomonadaceae_positive.csv', mode='a') as csv_pseudo:
                        csv_pseudo.write(row_str)
            flag = 1

        #write in file according to family
        if row[0] in list_couples_greg_id:
            couple_id = int(row[0])
            couple = CoupleJson.getCouplesByFilterParameter({'id':couple_id})
            bacterium = BacteriumJson.getByID(couple[0].bacterium)
            strain = StrainJson.getByID(bacterium.strain)
            specie = SpecieJson.getByID(strain.specie)
            genus = GenusJson.getByID(specie.genus)
            family = FamilyJson.getByID(genus.family)

            if family.designation == 'Staphylococcaceae' and couple.source_data_id != 3:
                with open('/home/christophe/Documents/Phages4A/CH_dataset_staphylococcaceae_positive.csv', mode='a') as csv_staph:
                    csv_staph.write(row_str)
            if family.designation == 'Moraxellaceae' and couple.source_data_id != 3:
                with open('/home/christophe/Documents/Phages4A/CH_dataset_moraxellaceae_positive.csv', mode='a') as csv_morax:
                    csv_morax.write(row_str)
            if family.designation == 'Enterobacteriaceae' and couple.source_data_id != 3:
                with open('/home/christophe/Documents/Phages4A/CH_dataset_enterobacteriaceae_positive.csv', mode='a') as csv_entero:
                    csv_entero.write(row_str)
            if family.designation == 'Pseudomonadaceae' and couple.source_data_id != 3:
                with open('/home/christophe/Documents/Phages4A/CH_dataset_pseudomonadaceae_positive.csv', mode='a') as csv_pseudo:
                    csv_pseudo.write(row_str)
    