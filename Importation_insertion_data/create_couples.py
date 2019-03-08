# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 15:10:40 2018

@author: Diogo Leite
"""


import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import csv
import math

from objects_new.Organisms_new import Organism
from objects_new.Couples_new import Couple
 
df = pd.read_excel('couples_interactions_only_int.xlsx', sheetname='Feuil1')
 
print("Column headings:")
print(df.columns)

dict_values_print = {}

dict_lysis_type = {}

dict_lysis_type['CL'] = 1
dict_lysis_type['SCL'] = 2
dict_lysis_type['OL'] = 3
dict_lysis_type['CL&2M'] = 4
dict_lysis_type['SCL&2M'] = 5
dict_lysis_type['OL&2M'] = 6

list_of_couples = []

for line in df.iterrows():
    print("---------------")
    #print(line)
    bact_id = line[1]['bact_ident']
    bact_name = line[1]['bactphage']
    p68_lysis = line[1]['P68']
    p44AHJD_lysis = line[1]['44AHJD']
    p3A_lysis = line[1]['3A']
    p71_lysis = line[1]['71P']
    p77_lysis = line[1]['77P']
    pK_lysis = line[1]['K']
    pSb1_lysis = line[1]['Sb-1']
    print(bact_id)
    print(bact_name)
    print(p77_lysis)
    print(pK_lysis)
    dict_value = Organism.get_organism_id_by_strain_like(bact_name)
    print(dict_value)

    id_bact = bact_id
    id_type_interaction = 1
    type_interaction = True
    id_level_interaction = 4

    #def __init__(self, id_couple, interact_pn,  fk_bacteria, fk_phage, fk_type_inter, fk_level_interact, fk_lysis_inter = -1, fk_source_data = -1):

    #Check if any couple exists
    id_phage = 10107
    if p68_lysis != 'No':
        print("3333333")
        print(p68_lysis)
        id_lysis = dict_lysis_type[p68_lysis]
        couple_obj = Couple(id_couple = -1, interact_pn = True, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_lysis_inter = id_lysis, fk_source_data = 3)
        print(couple_obj)
    elif p68_lysis == 'No':
        couple_obj = Couple(id_couple = -1, interact_pn = False, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_source_data = 3)

    list_of_couples.append(couple_obj)

    id_phage = 9393
    if p44AHJD_lysis != 'No':
        print("3333333")
        print(p68_lysis)
        id_lysis = dict_lysis_type[p44AHJD_lysis]
        couple_obj = Couple(id_couple = -1, interact_pn = True, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_lysis_inter = id_lysis, fk_source_data = 3)
        print(couple_obj)
    elif p44AHJD_lysis == 'No':
        couple_obj = Couple(id_couple = -1, interact_pn = False, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_source_data = 3)

    list_of_couples.append(couple_obj)

    id_phage = 9999
    if p3A_lysis != 'No':
        print("3333333")
        print(p68_lysis)
        id_lysis = dict_lysis_type[p3A_lysis]
        couple_obj = Couple(id_couple = -1, interact_pn = True, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_lysis_inter = id_lysis, fk_source_data = 3)
        print(couple_obj)
    elif p3A_lysis == 'No':
        couple_obj = Couple(id_couple = -1, interact_pn = False, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_source_data = 3)

    list_of_couples.append(couple_obj)

    id_phage = 8656
    if p71_lysis != 'No':
        print("3333333")
        print(p68_lysis)
        id_lysis = dict_lysis_type[p71_lysis]
        couple_obj = Couple(id_couple = -1, interact_pn = True, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_lysis_inter = id_lysis, fk_source_data = 3)
        print(couple_obj)
    elif p71_lysis == 'No':
        couple_obj = Couple(id_couple = -1, interact_pn = False, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_source_data = 3)

    list_of_couples.append(couple_obj)

    id_phage = 10058
    if p77_lysis != 'No':
        print("3333333")
        print(p68_lysis)
        id_lysis = dict_lysis_type[p77_lysis]
        couple_obj = Couple(id_couple = -1, interact_pn = True, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_lysis_inter = id_lysis, fk_source_data = 3)
        print(couple_obj)
    elif p77_lysis == 'No':
        couple_obj = Couple(id_couple = -1, interact_pn = False, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_source_data = 3)

    list_of_couples.append(couple_obj)

    id_phage = 8981
    if pK_lysis != 'No':
        print("3333333")
        print(p68_lysis)
        id_lysis = dict_lysis_type[pK_lysis]
        couple_obj = Couple(id_couple = -1, interact_pn = True, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_lysis_inter = id_lysis, fk_source_data = 3)
        print(couple_obj)
    elif pK_lysis == 'No':
        couple_obj = Couple(id_couple = -1, interact_pn = False, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_source_data = 3)

    list_of_couples.append(couple_obj)

    id_phage = 9216
    if pSb1_lysis != 'No':
        print("3333333")
        print(p68_lysis)
        id_lysis = dict_lysis_type[pSb1_lysis]
        couple_obj = Couple(id_couple = -1, interact_pn = True, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_lysis_inter = id_lysis, fk_source_data = 3)
        print(couple_obj)
    elif pSb1_lysis == 'No':
        couple_obj = Couple(id_couple = -1, interact_pn = False, fk_bacteria = id_bact, fk_phage = id_phage, fk_type_inter = id_type_interaction, fk_level_interact = id_level_interaction, fk_source_data = 3)

    list_of_couples.append(couple_obj)

for couples_obj in list_of_couples:
    couples_obj.create_couple()

print(len(list_of_couples))



