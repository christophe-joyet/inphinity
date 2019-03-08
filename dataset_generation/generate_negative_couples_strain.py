# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 11:53:32 2018

@author: Diogo
"""
import csv

from objects_new.Organisms_new import *
from objects_new.Sources_data_new import *
from objects_new.Couples_new import *
from objects_new.Species_new import *

from configuration.configuration_data import *

from taxonomy_tree import *


#In the first part I create only with the organism which are obtained from NCBI and PhageBD

list_id_phages = []

list_id_phages_NCBI = Organism.get_ids_all_organisms_by_params(2, 1)
list_id_phages_PhageDB = Organism.get_ids_all_organisms_by_params(2, 2)

list_id_bacteria = Organism.get_ids_all_organisms_by_params(1, 1)

list_id_phages = list(set(list_id_phages_NCBI + list_id_phages_PhageDB)) 

print(len(list_id_phages))
print(len(list_id_bacteria))

def get_id_bacteria_positive_couple_by_phge_id(fk_phage):
    """
    Get ids of bacteira from positive couple by a phage id

    :param fk_phage: phage id

    :type fk_phage: int - required 

    :return: list of bacteria ids
    :rtype list(int)
    """
    list_id_bacteria_attacked = Couple.get_all_positive_couples_by_phage_id(id_phage)

    list_id_bacteria = []
    for couple in list_id_bacteria_attacked:
        if couple.fk_bacteria not in list_id_bacteria:
            list_id_bacteria.append(couple.fk_bacteria)
    return list_id_bacteria

def get_all_bacteria_id_same_specie(id_bacterium):
    """
    return ids of bacteira from the same specie of id_bacterium

    :param id_bacterium: phage id

    :type id_bacterium: int - required 

    :return: list of bacteria ids
    :rtype list(int)
    """
    list_id_bact = []

    specie_obj = Specie.get_specie_by_bacterium_id(id_bacterium)
    list_ob_bacteria = Organism.get_all_bacteria_by_id_specie(specie_obj.id_specie)
    for bacterium in list_id_bact:
        if bacterium.id_organism not in list_id_bact:
            list_id_bact.append(bacterium.id_organism)
    return list_id_bact

def generate_couples_contraing_id_bact(phage_id, list_id_bact_not_allow):
    """
    return list of couples based on the constrain used for the negative couple

    :param phage_id: phage used in these interactions
    :param list_id_bact_not_allow: list of bacteria ids that are in same taxonomy (given in "get_all_bacteria_id_same_specie" method) of positive couple

    :type phage_id: int . required
    :type list_id_bact_not_allow: list(int) - required 

    :return: list of couples
    :rtype list(Couple)
    """
    list_couples = []
    list_id_bacteria_allowed = list(set(list_id_bacteria) - set(list_id_bact_not_allow))
    for id_bact in list_id_bacteria_allowed:
        list_couples.append(Couple(id_couple = -1, interact_pn = 0, fk_bacteria = id_bact, fk_phage = phage_id, fk_level_interact = 4, fk_type_inter = 4))
    return list_couples


# main method
def create_couple():
    """
    This method is used to create all possible couples

    """
    count_interactions = 0
    count_interactions_interm = 0
    list_id_bact_specie = []
    count_phages = 0


    for id_phage in list_id_phages:
        print("----------------------")
        list_id_bacteria_attacked = get_id_bacteria_positive_couple_by_phge_id(id_phage)
        print(len(list_id_bacteria_attacked))
        count_interactions_interm = 0
        for id_bacterium in list_id_bacteria_attacked:
            list_id_bact_specie = get_all_bacteria_id_same_specie(id_bacterium)
            aux = generate_couples_contraing_id_bact(id_phage, list_id_bact_specie)
            #print(len(aux))
            count_interactions += len(aux)
            count_interactions_interm += len(aux)
        count_phages +=1 
        print(count_phages)

        print(count_interactions_interm)
        print(count_interactions)

################### new one #####################
#def get_species_frequency_by_specie(specie_obj):


taxo_obj = Taxonomy_tree()
taxo_obj.generate_tree()
##taxo_obj.tree_taxo.show(data_property='id_db')
taxo_obj.tree_taxo.show()


dict_val = Specie.get_all_quantity_couple_by_specie_by_phage_id_positive(9999)
list_id_species_positive = []
print(dict_val)
for key, value in dict_val.items() :
    print (key, " ¦ Frequence: {0}".format(value))
    list_id_species_positive.append((key.id_specie, key.fk_genus, value))
print(list_id_species_positive)

list_others_species = Specie.get_all_species_by_genus_id(list_id_species_positive[0][1])
print(list_others_species)

#max_species = 0
#list_id_specie = []

#array_fre = []
#dict_fre = {}
#for id_phage in list_id_phages:
#    list_id_specie = []
#    list_couple = Couple.get_all_couples_by_phage_id(id_phage)
#    for couple in list_couple:
#        id_specie = Specie.get_specie_id_by_organism_id(couple.fk_bacteria)
#        if id_specie not in list_id_specie and couple.interact_pn == 1:
#            list_id_specie.append(id_specie)
#    if len(list_id_specie) > max_species:
#        max_species = len(list_id_specie)
#    try:
#        dict_fre[len(list_id_specie)].append(id_phage)
#    except KeyError:
#        dict_fre[len(list_id_specie)] = [id_phage]

#print(dict_fre)


#with open('dict.csv', 'w') as csv_file:
#    writer = csv.writer(csv_file)
#    for key, value in dict_fre.items():
#       writer.writerow([key, value])
    

#print(max_species)