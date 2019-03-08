# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 13:44:54 2018

@author: Diogo
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from objects_new.Couples_new import *
from objects_new.Protein_dom_new import *
from objects_new.Organisms_new import *
from objects_new.Proteins_new import *
from objects_new.Ddi_interactions_new import *
from objects_new.DDI_interactions_DB_new import *

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed

import multiprocessing

import numpy as np


def get_couples_NCBI_phagesDB():
    #Get all couples from NCBI and Phages DB
    list_couples = Couple.get_all_couples_by_type_level_source(1, 4, 1)
    list_couples = list_couples + Couple.get_all_couples_by_type_level_source(1, 4, 2)
    list_couples = list_couples + Couple.get_all_couples_by_type_level_source(0, 4, 4)
    return list_couples

def serach_domains(protein_id):
    global dict_prot_dom
    values = dict_prot_dom.get(protein_id, [])
    return protein_id, values


def get_list_domains(organism_id):
    max_proce = multiprocessing.cpu_count()
    max_proce = 1
    dict_domains = {}
    list_proteins = Protein.get_all_Proteins_by_organism_id(organism_id)
    with ThreadPoolExecutor(max_workers=max_proce) as executor:
        futures = [executor.submit(serach_domains, protein_obj.id_protein) for protein_obj in list_proteins]
        for future in as_completed(futures):
            protein_id, values = future.result()
            if len(values) > 0:
                dict_domains[protein_id] = values
    return dict_domains

def calculate_scores(list_db_domains):
    if 1 in list_db_domains or 2 in list_db_domains or 3 in list_db_domains:
        return 9
    else:
        return len(list_db_domains)

def combine_domaines(list_dom_bact, list_dom_phage):
    list_scores_ddi = []
    id_interaction = -1
    for dom_bact in list_dom_bact:
        for dom_phage in list_dom_phage:
            id_interaction = dict_ddi_id_inter.get((dom_bact, dom_phage), -1)
            if id_interaction == -1:
                id_interaction = dict_ddi_id_inter.get((dom_phage, dom_bact), -1)
            if id_interaction != -1:
                # An error here TypeError("unhashable type: 'list'",)
                list_db_domains = dict_interactions_ddi_db.get(id_interaction, [])
                score_ddi = calculate_scores(list_db_domains)
                list_scores_ddi.append(score_ddi)
    return list_scores_ddi


def create_bins_scores_SB(list_scores_values, size_bins):
    #Size of bins
    global max_value
    number_bins = int(max_value/ size_bins)
    range_bins = (max_value + 0.5) / number_bins
    range_bins_vec = np.arange(0, (max_value + range_bins), range_bins)
    list_of_quantites, bins_intervals, bins_size_draw  = plt.hist(list_scores_values, bins=range_bins_vec) 
    plt.close('all')
    return list_of_quantites, bins_intervals

def set_id_label_couple(id_couple, label_couple, vec_values):
    vec_scores = np.concatenate(([id_couple], vec_values, [label_couple]), axis = 0)
    return vec_scores

def save_DS(matrix_values, file_name, config_bins = False):
    nb_columns = len(matrix_values[0]) -2 
    string_head = "Id_interactions,"
    format_values = "%i,"
    aux_column_number = 0
    type_info = ""
    if config_bins == False:
        type_info = "%i,"
    else:
        type_info = "%1.4f,"

    while nb_columns > aux_column_number:
        string_head += "B" + str(aux_column_number) + ","
        format_values += type_info
        aux_column_number += 1

    format_values += "%i"
    string_head += "label"

    print(format_values)
    #https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.savetxt.html
    np.savetxt(file_name, matrix_values, delimiter=',', fmt=format_values, header=string_head)

list_couples = get_couples_NCBI_phagesDB()
print('End of the couples collection')
dict_ddi_id_inter = DDI_interaction.get_all_DDI_interaction_to_dictionary()
print('End of the DDI collection')
dict_interactions_ddi_db = DDI_interaction_DB.get_all_DDI_interactionDB_to_dictionary()
print('End of the interaction DDI collection')
dict_prot_dom = ProteinDom.get_all_protein_domain_dict()
print('End of the Prot_DOM collection')

list_domains_bact = []
list_domains_phage = []
list_scores_ppi = []


size_bins = 5
size_matrix = int(190 / 5)
max_value = 190

matrix_results_size_bins = np.empty([0, (size_matrix + 2)])


for couple_obj in list_couples:
    bacteria_id = couple_obj.fk_bacteria
    phage_id = couple_obj.fk_phage
    dict_bacterium = get_list_domains(bacteria_id)
    dict_phage = get_list_domains(phage_id)
    aux = 0
    for key_bact in dict_bacterium:
        list_domains_bact = dict_bacterium[key_bact]
        for key_phage in dict_phage:
            list_domains_phage = dict_phage[key_phage]
            list_scores = combine_domaines(list_domains_bact, list_domains_phage)
            list_scores_ppi = list_scores_ppi + list_scores
    print(len(list_scores_ppi))
    print('End of the couple {0}'.format(couple_obj.id_couple))

    bins_scores_SB, bins_size_SB = create_bins_scores_SB(list_scores_ppi, size_bins)


    bins_scores_SB = set_id_label_couple(couple_obj.id_couple, couple_obj.interact_pn, bins_scores_SB)

    matrix_results_size_bins = np.append(matrix_results_size_bins, [bins_scores_SB], axis = 0)
    list_scores_ppi = []


save_DS(matrix_results_size_bins, 'SB_' + str(size_bins) + '_ML.csv', False)


