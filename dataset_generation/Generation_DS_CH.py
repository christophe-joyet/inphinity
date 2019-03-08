# -*- coding: utf-8 -*-
"""
Created on Sun May 6 09:32:54 2018

@author: Diogo
"""

from objects_new.Couples_new import *
from objects_new.Proteins_new import *
import numpy as np

def calculatePercenteAAChWeight(sequence):
    isCaract = False;
    Weight = 0
    PercCaract = np.zeros(21)
    #0 = C
    #1 = H
    #2 = O
    #3 = N
    #4 = S
    PercChimical = np.zeros(5)
    for caracter in sequence:
        caracter = caracter.upper()
        isCaract = False;
        if caracter == 'A':
            PercCaract[0] = PercCaract[0] + 1
            PercChimical[0] = PercChimical[0] + 3
            PercChimical[1] = PercChimical[1] + 7
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 89
            isCaract = True
        if caracter == 'C':
            PercCaract[1] = PercCaract[1] + 1
            PercChimical[0] = PercChimical[0] + 3
            PercChimical[1] = PercChimical[1] + 7
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 1
            Weight = Weight + 133
            isCaract = True
        if caracter == 'D':
            PercCaract[2] = PercCaract[2] + 1
            PercChimical[0] = PercChimical[0] + 4
            PercChimical[1] = PercChimical[1] + 6
            PercChimical[2] = PercChimical[2] + 4
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 133
            isCaract = True
        if caracter == 'E':
            PercCaract[3] = PercCaract[3] + 1
            PercChimical[0] = PercChimical[0] + 5
            PercChimical[1] = PercChimical[1] + 8
            PercChimical[2] = PercChimical[2] + 4
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 147
            isCaract = True
        if caracter == 'F':
            PercCaract[4] = PercCaract[4] + 1
            PercChimical[0] = PercChimical[0] + 3
            PercChimical[1] = PercChimical[1] + 6
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 165
            isCaract = True
        if caracter == 'G':
            PercCaract[5] = PercCaract[5] + 1
            PercChimical[0] = PercChimical[0] + 2
            PercChimical[1] = PercChimical[1] + 5
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 75
            isCaract = True
        if caracter == 'H':
            PercCaract[6] = PercCaract[6] + 1
            PercChimical[0] = PercChimical[0] + 6
            PercChimical[1] = PercChimical[1] + 9
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 3
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 155
            isCaract = True
        if caracter == 'I':
            PercCaract[7] = PercCaract[7] + 1
            PercChimical[0] = PercChimical[0] + 6
            PercChimical[1] = PercChimical[1] + 13
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 131
            isCaract = True
        if caracter == 'K':
            PercCaract[8] = PercCaract[8] + 1
            PercChimical[0] = PercChimical[0] + 6
            PercChimical[1] = PercChimical[1] + 15
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 2
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 146
            isCaract = True
        if caracter == 'L':
            PercCaract[9] = PercCaract[9] + 1
            PercChimical[0] = PercChimical[0] + 6
            PercChimical[1] = PercChimical[1] + 13
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 131
            isCaract = True
        if caracter == 'M':
            PercCaract[10] = PercCaract[10] + 1
            PercChimical[0] = PercChimical[0] + 5
            PercChimical[1] = PercChimical[1] + 11
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 1
            Weight = Weight + 149
            isCaract = True
        if caracter == 'N':
            PercCaract[11] = PercCaract[11] + 1
            PercChimical[0] = PercChimical[0] + 4
            PercChimical[1] = PercChimical[1] + 8
            PercChimical[2] = PercChimical[2] + 3
            PercChimical[3] = PercChimical[3] + 2
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 132
            isCaract = True
        if caracter == 'P':
            PercCaract[12] = PercCaract[12] + 1
            PercChimical[0] = PercChimical[0] + 5
            PercChimical[1] = PercChimical[1] + 9
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 115
            isCaract = True
        if caracter == 'Q':
            PercCaract[13] = PercCaract[13] + 1
            PercChimical[0] = PercChimical[0] + 5
            PercChimical[1] = PercChimical[1] + 10
            PercChimical[2] = PercChimical[2] + 3
            PercChimical[3] = PercChimical[3] + 2
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 146
            isCaract = True
        if caracter == 'R':
            PercCaract[14] = PercCaract[14] + 1
            PercChimical[0] = PercChimical[0] + 6
            PercChimical[1] = PercChimical[1] + 15
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 4
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 174
            isCaract = True
        if caracter == 'S':
            PercCaract[15] = PercCaract[15] + 1
            PercChimical[0] = PercChimical[0] + 3
            PercChimical[1] = PercChimical[1] + 7
            PercChimical[2] = PercChimical[2] + 3
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 105
            isCaract = True
        if caracter == 'T':
            PercCaract[16] = PercCaract[16] + 1
            PercChimical[0] = PercChimical[0] + 4
            PercChimical[1] = PercChimical[1] + 9
            PercChimical[2] = PercChimical[2] + 3
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 119
            isCaract = True
        if caracter == 'V':
            PercCaract[17] = PercCaract[17] + 1
            PercChimical[0] = PercChimical[0] + 5
            PercChimical[1] = PercChimical[1] + 11
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 117
            isCaract = True
        if caracter == 'W':
            PercCaract[18] = PercCaract[18] + 1
            PercChimical[0] = PercChimical[0] + 5
            PercChimical[1] = PercChimical[1] + 8
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 2
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 204
            isCaract = True
        if caracter == 'Y':
            PercCaract[19] = PercCaract[19] + 1
            PercChimical[0] = PercChimical[0] + 3
            PercChimical[1] = PercChimical[1] + 7
            PercChimical[2] = PercChimical[2] + 2
            PercChimical[3] = PercChimical[3] + 1
            PercChimical[4] = PercChimical[4] + 0
            Weight = Weight + 181
            isCaract = True
        if isCaract == False:
            PercCaract[20] = PercCaract[20] + 1
        isCaract = False
    PercCaract = np.divide(PercCaract, np.sum(PercCaract))
    PercChimical = np.divide(PercChimical, np.sum(PercChimical))
    results = np.concatenate(([Weight], PercCaract, PercChimical), axis=0)
    #results = np.concatenate(([Weight], results), axis = 0)
    #return Weight, PercCaract, PercChimical
    return results

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

def get_values_for_list_prots(list_proteins):
    qty_proteins = len(list_proteins)
    matrix_values = np.empty([1, 27])
    protein_number = 0

    for protein in list_proteins:
        results = calculatePercenteAAChWeight(protein.sequence_prot)
        matrix_values = np.append(matrix_values, [results], axis=0)
        protein_number += 1
    matrix_values = np.delete(matrix_values, 0, 0)
    return matrix_values

def creat_total_matrix_bact_phage(matrix_values_bact, matrix_values_phage):
    matrix_values_results_temp = np.empty([0, 54])
    matrix_values_results_final = np.empty([0, 54])
    qty_row_bact = matrix_values_bact.shape[0]
    qty_row_phage = matrix_values_phage.shape[0]
    aux_row_bact = 0
    aux_row_phage = 0
    matrix_PPI_bact = np.tile(matrix_values_bact, (qty_row_phage,1))
    matrix_PPI_Phage = np.tile(matrix_values_phage, (qty_row_bact,1))

    print(matrix_values_bact.shape)
    print(matrix_values_phage.shape)

    print(matrix_PPI_bact.shape)
    print(matrix_PPI_Phage.shape)

    #np.savetxt('testA.csv', matrix_values_bact, delimiter=',') 
    #np.savetxt('testB.csv', matrix_values_phage, delimiter=',') 
    #np.savetxt('testC.csv', matrix_PPI_bact, delimiter=',') 
    #np.savetxt('testD.csv', matrix_PPI_Phage, delimiter=',') 

    #while qty_row_bact > aux_row_bact:

    #    while qty_row_phage > aux_row_phage:
    #        matrix_values_results_temp = np.concatenate((matrix_values_bact[aux_row_bact], matrix_values_phage[aux_row_phage]), axis=0)
    #        matrix_values_results_final = np.append(matrix_values_results_final, [matrix_values_results_temp], axis = 0)
    #        aux_row_phage += 1
    #    aux_row_phage = 0
    #    aux_row_bact += 1
    matrix_values_results_final = np.concatenate((matrix_PPI_bact, matrix_PPI_Phage), axis = 1)
    print(matrix_values_results_final.shape)

    #np.savetxt('testE.csv', matrix_values_results_final, delimiter=',') 

    vector_mean = np.mean(matrix_values_results_final, axis=0)
    vector_std = np.std(matrix_values_results_final, axis=0)
    return vector_mean, vector_std

def save_DS(matrix_values, file_name):
    #https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.savetxt.html
    np.savetxt(file_name, matrix_values, delimiter=',', fmt='%1.14f')

#list_couples = get_couples_NCBI_phagesDB()
#list_couples = get_couples_Greg()


#
couple_obj = Couple.get_couples_by_list_id([35290, 35343 , 35343, 35458, 35876 , 36041, 36348, 36397])

#

list_couples = couple_obj
print(len(list_couples))
print("qwqwqwqwqwqw")

previous_phage = -1
previous_bact = -1
list_prot_phage = []
list_prot_bacterium = []
matrix_values_prot_Phage = None
matrix_values_prot_Bact = None
qty_couples = len(list_couples)

matrix_results = np.empty([0, 110])

for couple in list_couples:
    print(couple)
    if couple.fk_phage != previous_phage:
        previous_phage = couple.fk_phage
        list_prot_phage = Protein.get_all_Proteins_by_organism_id(couple.fk_phage)
        matrix_values_prot_Phage = get_values_for_list_prots(list_prot_phage)
    print(matrix_values_prot_Phage.shape)

    if couple.fk_bacteria != previous_bact:
        previous_bact = couple.fk_bacteria
        list_prot_bacterium = Protein.get_all_Proteins_by_organism_id(couple.fk_bacteria)
        matrix_values_prot_Bact = get_values_for_list_prots(list_prot_bacterium)
    print(matrix_values_prot_Bact.shape)

    vec_mean, vec_std = creat_total_matrix_bact_phage(matrix_values_prot_Bact, matrix_values_prot_Phage)

    id_interaction = couple.id_couple
    lable_interaction = couple.interact_pn
    vec_features_value = np.concatenate(([id_interaction], vec_mean, vec_std, [lable_interaction]), axis = 0)

    matrix_results = np.append(matrix_results, [vec_features_value], axis = 0)


    save_DS(matrix_results, 'DS_CH_student_version_error_ds.csv')