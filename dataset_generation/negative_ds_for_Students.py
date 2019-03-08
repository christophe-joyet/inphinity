# -*- coding: utf-8 -*-
"""
Created on Wed May 1 12:02:32 2018

@author: Diogo
"""

from objects_new.Organisms_new import *
from objects_new.Couples_new import *
from objects_new.Species_new import *

from random import randint

list_id_bact = []
list_id_phage = []

def get_values_bact_ids():
    #get all values of bacterias
    list_id_bact = Organism.get_ids_all_organisms_by_params(type_organism = 1, source_data = 1)
    return list_id_bact

def get_values_phage_ids():
    #get all values of phages
    list_id_phage = Organism.get_ids_all_organisms_by_params(type_organism = 2, source_data = 1)
    list_id_phage = list_id_phage + Organism.get_ids_all_organisms_by_params(type_organism = 2, source_data = 2)
    return list_id_phage

def get_specie_by_bacterim_id(id_bacterium):
    #get organism Bact
    organism_obj = Organism.get_organism_by_id(id_bacterium)
    specie_bact_orga = Specie.get_specie_by_bacterium_id(organism_obj.id_organism)
    print(organism_obj)
    return specie_bact_orga


def get_couples_by_bacterium_id_strain_level_four(id_bacterium):
    #get couples of bacteria
    list_couples = Couple.get_all_couples_positiv_by_bacterium_level(id_bacterium, 4)
    return list_couples

list_id_bact = get_values_bact_ids()
list_id_phage = get_values_phage_ids()
list_couples = []
list_id_phages_negative = []

aux_quantities = []
print(len(list_id_bact))

for id_bacterium in list_id_bact:

    specie_obj_bact = get_specie_by_bacterim_id(id_bacterium)
    list_couples = get_couples_by_bacterium_id_strain_level_four(id_bacterium)
    quantity_couples = len(list_couples)

    count_couple = 0
    max_value = len(list_id_phage) -1
    list_id_phages_negative = []
    list_id_phages_positive = []

    while count_couple < quantity_couples:
        #random phage_id
        position_phage_list = randint(0, max_value)
        id_phage = list_id_phage[position_phage_list]
        list_couples_phage = Couple.get_all_positive_couples_by_phage_id_level_id(id_phage, 4)
        #Confirm that any phage attack a bacterium with the same ID of specie
        qty_couples = len(list_couples_phage)
        qty_couples_count = 0
        for couple_obj in list_couples_phage:
            specie_obj = Specie.get_specie_by_bacterium_id(couple_obj.fk_bacteria)
            if specie_obj.id_specie != specie_obj_bact.id_specie:
                qty_couples_count += 1   
        id_couple = Couple.verify_couple_exist_by_phage_bact(id_bacterium, id_phage)
        if qty_couples_count == qty_couples and id_couple == -1 and id_phage not in list_id_phages_negative:
            list_id_phages_negative.append(id_phage)
            count_couple += 1  
            print(len(list_id_phages_negative))
           
    for id_phage in list_id_phages_negative:
        couple_obj = Couple(id_couple = -1, interact_pn = 0, fk_bacteria = id_bacterium, fk_phage = id_phage, fk_type_inter = 4, fk_level_interact = 4, fk_source_data = 4)
        couple_obj.create_couple()
    print("--------------------")
    aux_quantities.append((id_bacterium, quantity_couples, len(list_id_phages_negative)))
    print(quantity_couples)
    print(len(list_id_phages_negative))
    list_id_phages_negative = []

aux = 123123
#print(list_id_bact[1])
#print(len(list_couples))
#quantity_couples = len(list_couples)
#count_couple = 0

#max_value = len(list_id_phage)

#list_id_phages_negative = []
#list_id_phages_positive = []
##Creation of the same quantity  of negative couples
#while count_couple < quantity_couples:
#    #random phage_id
#    position_phage_list = randint(0, max_value)
#    id_phage = list_id_phage[position_phage_list]
#    list_couples_phage = Couple.get_all_positive_couples_by_phage_id_level_id(id_phage, 4)
#    #Confirm that any phage attack a bacterium with the same ID
#    qty_couples = len(list_couples_phage)
#    qty_couples_count = 0
#    for couple_obj in list_couples_phage:
#        specie_obj = Specie.get_specie_by_bacterium_id(couple_obj.fk_bacteria)
#        if specie_obj.id_specie != specie_bact_orga.id_specie:
#            qty_couples_count += 1
#    id_couple = Couple.verify_couple_exist_by_phage_bact(organism_obj.id_organism, id_phage, 0)
#    if qty_couples_count == qty_couples and id_couple == -1:
#        list_id_phages_negative.append(id_phage)
#        count_couple += 1

##insert couples in the database
#couple_obj = Couple(id_couple = -1, interact_pn = 0, fk_bacteria = list_id_bact[1], fk_phage = list_id_phages_negative[0], fk_type_inter = 4, fk_level_interact = 4)

#print(len(list_id_phages_negative))

#print(list_id_phages_negative)

##delete
#for couple_obj in list_couples:
#    list_id_phages_positive.append(couple_obj.id_couple)
#print(list_id_phages_positive)
