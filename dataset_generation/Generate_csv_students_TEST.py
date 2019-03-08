# -*- coding: utf-8 -*-
"""
Created on Fri May 4 14:51:54 2018

@author: Diogo
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from objects_new.PPI_preview_new import *
from objects_new.Couples_new import *
from objects_new.Proteins_new import *
import time
import copy


def get_couples_NCBI_phagesDB():
    #Get all couples from NCBI and Phages DB
    list_couples = Couple.get_all_couples_by_type_level_source(1, 4, 1)
    list_couples = list_couples + Couple.get_all_couples_by_type_level_source(1, 4, 2)
    list_couples = list_couples + Couple.get_all_couples_by_type_level_source(0, 4, 4)
    return list_couples

def get_couples_Greg():
    #Get all couples from NCBI and Phages DB
    list_couples = Couple.get_all_couples_by_type_level_source(1, 4, 3)
    list_couples = list_couples + Couple.get_all_couples_by_type_level_source(0, 4, 3)
    return list_couples

def create_bins_scores_NB(list_scores_values, qty_bins, only_zeros = False):
    #Number of bins
    list_of_quantites, bins_intervals, bins_size_draw  = plt.hist(list_scores_values, bins=qty_bins) 
    plt.close('all')
    if only_zeros == True:
        position = np.where(list_of_quantites != 0)
        position = position[0][0]
        list_of_quantites[0], list_of_quantites[position] = list_of_quantites[position], list_of_quantites[0]
        qty_bins = len(bins_intervals)
        bins_intervals = np.zeros(qty_bins)

    return list_of_quantites, bins_intervals

def create_bins_scores_SB(list_scores_values, size_bins):
    #Size of bins
    global max_value
    number_bins = int(max_value/ size_bins)
    range_bins = (max_value + 0.5) / number_bins
    range_bins_vec = np.arange(0, (max_value + range_bins), range_bins)
    list_of_quantites, bins_intervals, bins_size_draw  = plt.hist(list_scores_values, bins=range_bins_vec) 
    plt.close('all')
    return list_of_quantites, bins_intervals



def set_zeros_values_vec(list_scores_values, n_times):
    #To add v, n times, to l:
    #l += n * [v]
    list_with_zeros = copy.copy(list_scores_values)
    list_with_zeros += n_times * [0]
    return list_with_zeros

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

def set_id_label_couple(id_couple, label_couple, vec_values):
    vec_scores = np.concatenate(([id_couple], vec_values, [label_couple]), axis = 0)
    return vec_scores

def save_DS_append(matrix_values, file_name, config_bins = False, bool_header = False):
    f_handle = open(file_name, 'a')

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
    if bool_header == True:
        np.savetxt(file_name, matrix_values, delimiter=',', fmt=format_values, header=string_head)
    else:
        np.savetxt(file_name, matrix_values, delimiter=',', fmt=format_values)

def simple_save_file(file_name, array_values, config_bins = False, bool_header = False):
    f_handle = open(file_name, 'a+')

    nb_columns = len(array_values) -2 
    string_head = "Id_interactions,"
    format_values = "%i,"
    aux_column_number = 0
    type_info = ""
    if config_bins == False:
        type_info = "%i,"
    else:
        type_info = "%1.4f, "

    while nb_columns > aux_column_number:
        string_head += "B" + str(aux_column_number) + ","
        format_values += type_info
        aux_column_number += 1

    format_values += "%i"
    string_head += "label"

    #https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.savetxt.html
    strin_value = ''
    if bool_header == True:
        np.savetxt(file_name, [array_values], delimiter=',', fmt=format_values, header=string_head)
    else:
        string_value = np.array2string(array_values, separator=',', formatter = {'float_kind':lambda array_values: "%i" % array_values})
        string_value = string_value.replace(" ", "")
        string_value = string_value.replace("\n", "")[1:-1]
        np.savetxt(file_name, [array_values], delimiter=',',fmt=format_values )

    with open('aaaa' + file_name, 'a+') as file:
        file.write(string_value + "\n")
    f_handle.close()


def reset_matrix():
    global matrix_results_qty_bins
    global matrix_results_size_bins

    global matrix_results_qty_bins_BV
    global matrix_results_size_bins_BV

    global matrix_results_qty_bins_ZERO
    global matrix_results_size_bins_ZERO

    global matrix_results_qty_bins_BV_ZERO
    global matrix_results_size_bins_BV_ZERO

    matrix_results_qty_bins = np.empty([0, (n_bins + 2)])
    matrix_results_size_bins = np.empty([0, (size_bins + 2)])

    matrix_results_qty_bins_BV = np.empty([0,  (n_bins + 3)])
    matrix_results_size_bins_BV = np.empty([0, (size_bins + 3)])

    matrix_results_qty_bins_ZERO = np.empty([0, (n_bins + 2)])
    matrix_results_size_bins_ZERO = np.empty([0, (size_bins + 2)])

    matrix_results_qty_bins_BV_ZERO = np.empty([0, (n_bins + 3)])
    matrix_results_size_bins_BV_ZERO = np.empty([0, (size_bins + 3)])

#Finir le script pour générer toutes les intéractions....

#list_of_scores = PPI_preview.get_ppi_preview_scores_grouped_by_couple_id(11949)

list_possible_values = PPI_preview.get_all_ppi_preview_couple()
print(len(list_possible_values))
#print(list_possible_values)
########### A voir mais max = 108
#max_value = PPI_preview.get_max_ppi_score()
max_value = 162

qty_columns = max_value + 2
list_couples = get_couples_Greg()

print(len(list_couples))

list_axu = list(set(list_couples))
print(len(list_axu))

n_bins = 108
size_bins = 54

qty_bins = int(max_value/size_bins)

matrix_results_qty_bins = np.empty([0, (n_bins + 2)])
matrix_results_size_bins = np.empty([0, (qty_bins + 2)])

matrix_results_qty_bins_BV = np.empty([0,  (n_bins + 3)])
matrix_results_size_bins_BV = np.empty([0, (qty_bins + 3)])

matrix_results_qty_bins_ZERO = np.empty([0, (n_bins + 2)])
matrix_results_size_bins_ZERO = np.empty([0, (qty_bins + 2)])

matrix_results_qty_bins_BV_ZERO = np.empty([0, (n_bins + 3)])
matrix_results_size_bins_BV_ZERO = np.empty([0, (qty_bins + 3)])

id_previou_bact = -1
id_previou_phage = -1

qty_proteins_bact = 0
qty_proteins_phage = 0

aux_count = 0
only_zeros = False

first_copy = True


aux = 0
start_time = time.time()
for couple in list_couples:
    time_start_execution = time.time()
    only_zeros = False
    if couple.fk_bacteria != id_previou_bact:
        qty_proteins_bact = Protein.get_qty_proteins_by_organism_id(couple.fk_bacteria)
        id_previou_bact = couple.fk_bacteria
        print("Change bacterium")

    if couple.fk_phage != id_previou_phage:
        qty_proteins_phage = Protein.get_qty_proteins_by_organism_id(couple.fk_phage)
        id_previou_phage = couple.fk_phage

    max_protein = qty_proteins_bact * qty_proteins_phage

    #print("max proteins %d" %max_protein)






    list_of_scores = PPI_preview.get_ppi_preview_scores_grouped_by_couple_id(couple.id_couple)
    list_of_scores_zeros = []
    ###############
    ppi_score_zero = max_protein - len(list_of_scores)
    if(len(list_of_scores) == 0):
        only_zeros = True
        list_of_scores_zeros = set_zeros_values_vec(list_of_scores_zeros, ppi_score_zero)
    else:
        list_of_scores_zeros = set_zeros_values_vec(list_of_scores, ppi_score_zero)

    ##Without zeros

    bins_scores_NB, bins_size_NB = create_bins_scores_NB(list_of_scores, n_bins, False)
    bins_scores_SB, bins_size_SB = create_bins_scores_SB(list_of_scores, size_bins)
    

    bins_scores_NB = set_id_label_couple(couple.id_couple, couple.interact_pn, bins_scores_NB)
    bins_scores_SB = set_id_label_couple(couple.id_couple, couple.interact_pn, bins_scores_SB)

    matrix_results_qty_bins = np.append(matrix_results_qty_bins, [bins_scores_NB], axis = 0)
    matrix_results_size_bins = np.append(matrix_results_size_bins, [bins_scores_SB], axis = 0)

    ##With zeros
    bins_scores_NB_ZERO, bins_size_NB_ZERO = create_bins_scores_NB(list_of_scores_zeros, n_bins, only_zeros)
    bins_socres_SB_ZERO, bins_size_SB_ZERO = create_bins_scores_SB(list_of_scores_zeros, size_bins)

    bins_scores_NB_ZERO = set_id_label_couple(couple.id_couple, couple.interact_pn, bins_scores_NB_ZERO)
    bins_socres_SB_ZERO = set_id_label_couple(couple.id_couple, couple.interact_pn, bins_socres_SB_ZERO)


    matrix_results_qty_bins_ZERO = np.append(matrix_results_qty_bins_ZERO, [bins_scores_NB_ZERO], axis = 0)
    matrix_results_size_bins_ZERO = np.append(matrix_results_size_bins_ZERO, [bins_socres_SB_ZERO], axis = 0)




#    simple_save_file('NB_ZERO_' + str(n_bins) + '_ML.csv', bins_scores_NB_ZERO, False)

    #############
    aux_count += 1

    del bins_scores_NB_ZERO
    del bins_size_NB_ZERO
    time_end_execution = time.time()
    time_execution = time_end_execution - time_start_execution
    print("%5d: takes: %5.4f seconds to execute" % (aux, time_execution))
    aux += 1
end_time = time.time()
total_time = end_time - start_time
print("how many time? : %8.6f" % total_time)

print("Saving files")

save_DS(matrix_results_qty_bins, 'GREG_NB_' + str(n_bins) + '_ML.csv', False)
save_DS(matrix_results_size_bins, 'GREG_SB_' + str(size_bins) + '_ML.csv', False)
#save_DS(matrix_results_qty_bins_BV, 'NB_bins_config_' + str(n_bins) + '_ML.csv', True)
#save_DS(matrix_results_size_bins_BV, 'SB_bins_config_' + str(size_bins) + '_ML.csv', True)

save_DS(matrix_results_qty_bins_ZERO, 'GREG_NB_ZERO_' + str(n_bins) + '_ML.csv', False)
save_DS(matrix_results_size_bins_ZERO, 'GREG_SB_ZERO_' + str(size_bins) + '_ML.csv', False)
#save_DS(matrix_results_qty_bins_BV_ZERO, 'NB_bins_config_ZERO_' + str(n_bins) + '_ML.csv', True)
#save_DS(matrix_results_size_bins_BV_ZERO, 'SB_bins_config_ZERO_' + str(size_bins) + '_ML.csv', True)




print(len(list_of_scores))
print(list_of_scores)
print(type(list_of_scores))

list_of_quantites, bins_intervals, bins_size_draw  = plt.hist(list_of_scores, bins=6) 
print(list_of_quantites)
print(bins_intervals)

