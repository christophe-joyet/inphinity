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

flag = 0
with open('/home/christophe/Documents/Phages4A/CH_dataset.csv') as csv_file: 
    csv_reader = csv.reader(csv_file, delimiter=',')
    
    #add labels row
    for row in csv_reader:
        row_str = ""
        for el in row:
            row_str += el + ','
        row_str += '\n'
        break

    if flag == 0:
        with open('/home/christophe/Documents/Phages4A/CH_dataset_greg_couples.csv', mode='a') as greg_csv:
                greg_csv.write(row_str)
        #with open('/home/christophe/Documents/Phages4A/CH_dataset_greg_couples_positive_interaction.csv', mode='a') as positive_csv:
        #            positive_csv.write(row_str)
        #with open('/home/christophe/Documents/Phages4A/CH_dataset_greg_couples_negative_interaction.csv', mode='a') as negative_csv:
        #            negative_csv.write(row_str)
        flag = 1
    #end add labels row

    #file creation
    for row in csv_reader:
        if row[0] in list_couples_greg_id:
            #transform row in str
            row_str = ""
            for el in row:
                row_str += el + ','
            row_str += '\n'
            with open('/home/christophe/Documents/Phages4A/CH_dataset_greg_couples.csv', mode='a') as greg_csv:
                greg_csv.write(row_str)
            #if row[-1] == 'True':
            #    with open('/home/christophe/Documents/Phages4A/CH_dataset_greg_couples_positive_interaction.csv', mode='a') as positive_csv:
            #        positive_csv.write(row_str)
            #else:
            #    with open('/home/christophe/Documents/Phages4A/CH_dataset_greg_couples_negative_interaction.csv', mode='a') as negative_csv:
            #        negative_csv.write(row_str)


#Extract lysis value 
lysis_dictionnary = {constants.CLEAR_LYSIS: 0, 
                     constants.SEMI_CLEAR_LYSIS: 0,
                     constants.OPAQUE_LYSIS:0,
                     constants.CLEAR_LYSIS_1E7PLUS:0, 
                     constants.CLEAR_LYSIS_1E7MINUS:0,
                     constants.SEMI_CLEAR_LYSIS_1E7PLUS:0,
                     constants.SEMI_CLEAR_LYSIS_1E7MINUS:0}

for couple in list_couples_greg:
    if couple.lysis in lysis_dictionnary.keys():
        lysis_dictionnary[couple.lysis] += 1

lysis_name = ['Clear Lysis\nNumber of phages unknown', 'Semi Clear Lysis\nNumber of phages unknown', 'Opaque Lysis', 'Clear Lysis\nNumber of phages > 1E7', 'Clear Lysis\nNumber of phages < 1E7', 'Semi Clear Lysis\nNumber of phages > 1E7', 'Semi Clear Lysis\nNumber of phages < 1E7']
lysis_value = [0] * len(lysis_name)

i = 0
for value in lysis_dictionnary.values():
    lysis_value[i] = value
    i += 1


# Graph : number of phages and bacteria
labels = 'Phages\n(' + str(len(list_phages)) + ')', 'Bacteria\n(' + str(len(list_bacteria)) + ')'
size = len(list_phages), len(list_bacteria)
colors = ['gold', 'lightskyblue']
plt.pie(size, labels=labels, colors=colors, shadow=True, autopct='%1.1f%%')
plt.axis('equal')
plt.show()

# Graph : interactions positive negative
labels=['Couples with positive interaction\n(' + str(len(list_couples_positive_interaction)) + ')', 'Couples with negative interaction\n(' + str(len(list_couples_negative_interaction)) + ')']
size= len(list_couples_positive_interaction), len(list_couples_negative_interaction)
colors = ['gold', 'lightskyblue']
plt.pie(size, labels=labels, colors=colors, shadow=True, autopct='%1.1f%%')
plt.axis('equal')
plt.show()

#Graph : Lysis 
size = []
labels = []
for value in lysis_dictionnary.values():
    size.append(value)
labels = lysis_name
colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', '#FFEC00', '#52D726', '#007ED6']
plt.pie(size, colors=colors, labels=labels, shadow=True, autopct='%1.1f%%')
plt.axis('equal')
plt.show()


# Dictionnary {bacterie : nbr of bacteries}
bacterie_found = 0
bacteries_dictionnary = {}
for couple in list_couples_greg:
    bacterie_found = 0
    for bacterie in bacteries_dictionnary.keys():
            # Check if there is already the id of the bacterium in the dictionnary
            if couple.bacterium == bacterie.id:
                    # Increase de number of the bacterie whith this id
                    bacteries_dictionnary[bacterie] += 1
                    bacterie_found = 1
                    break

    # If we hadn't found the bacterie we add it in the dictionnary
    if bacterie_found != 1 :
        bacteries_dictionnary[BacteriumJson.getByID(couple.bacterium)] = 1

print(len(bacteries_dictionnary))

# Dictionnary {strain_name : number of bacteries in this strain}
strain_dictionnary = {}
strain_found = 0
for bacterie, bacterie_value in bacteries_dictionnary.items():
    # Get the strain_id from the bacterie
    strain_id = bacterie.strain
    for strain in strain_dictionnary.keys():
        strain_found = 0
        # Check if there is already the id of the strain in the dictionnary
        if strain_id == strain.id:
            # Increase de number of the bacterie of this strain
            strain_dictionnary[strain] += bacterie_value
            strain_found = 1
            break
    
    # If we hadn't found the strain we add it in the dictionnary
    if strain_found != 1 :
        # Add the key-value : strain
        strain_dictionnary[StrainJson.getByID(bacterie.strain)] = bacterie_value
        specie_dictionnary = {}

# Dictionnary {specie : number of strain belong to this specie}
specie_dictionnary = {}
specie_found = 0

for strain, strain_value in strain_dictionnary.items():
    # Get the specie_id from the strain
    specie_id = strain.specie
    for specie in specie_dictionnary.keys():
        specie_found = 0
        # Check if there is already the specie in the dictionnary
        if specie_id == specie.id:
            # Increase de number of the strain of this specie
            specie_dictionnary[specie] += strain_value
            specie_found = 1
            break
    # If we hadn't found the specie we add it in the dictionnary
    if specie_found != 1 :
        # Add the key-value : specie
        specie_dictionnary[SpecieJson.getByID(strain.specie)] = strain_value

print(len(specie_dictionnary))

# Dictionnary {genus : number of species belong to this genus}
genus_dictionnary = {}
genus_found = 0

for specie, species_value in specie_dictionnary.items():
    # Get the genus id from the specie
    genus_id = specie.genus
    for genus in genus_dictionnary.keys():
        genus_found = 0
        # Check if there is already the genus in the dictionnary
        if genus_id == genus.id:
            # Increase de number of the specie of this genus
            genus_dictionnary[genus] += species_value
            genus_found = 1
            break
    # If we hadn't found the genus we add it in the dictionnary
    if genus_found != 1 :
        # Add the key-value : genus
        genus_dictionnary[GenusJson.getByID(specie.genus)] = species_value

print(len(genus_dictionnary))

# Dictionnary {family : number of genus in this family}
family_dictionnary = {}
family_found = 0

for genus, genus_value in genus_dictionnary.items():
    # Get the family id from the genus
    family_id = genus.family
    for family in family_dictionnary.keys():
        family_found = 0
        # Check if there is already the family in the dictionnary
        if family_id == family.id:
            # Increase de number of the genus belong to this family
            family_dictionnary[family] += genus_value
            family_found = 1
            break
    # If we hadn't found the family we add it in the dictionnary
    if family_found != 1 :
        # Add the key-value : family
        family_dictionnary[FamilyJson.getByID(genus.family)] = genus_value


# Family repartition chart
number_of_family_tab = []
family_name = []
for family, number_of_family in family_dictionnary.items():
    # Take only specie with more than one souche
    number_of_family_tab.append(number_of_family)
    family_name.append(family.designation)

print(number_of_family_tab)
print(family_name)

# Display the plot
size = []
for el in family_dictionnary.values():
    size.append(el)
labels=family_name
colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', '#FFEC00', '#52D726', '#007ED6']
plt.pie(size, labels=labels, colors=colors, autopct='%1.1f%%', shadow=True)
plt.axis('equal')
plt.show()