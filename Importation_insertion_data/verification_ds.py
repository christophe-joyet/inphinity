# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 10:57:10 2017

@author: Stage
"""

import csv

from objects_new.Couples_new import Couple

def read_csv(file_name):
    aux = 0
    list_couple_id = []
    with open(file_name, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='\n', quoting= csv.QUOTE_NONNUMERIC)
        for row in spamreader:
            row = [int(i) for i in row[:-1]]
            qty_element = row.count(0)
            if qty_element > 106:
                aux += 1
                list_couple_id.append(row[0])
    return list_couple_id

def write_ids_organism(list_organism):
    with open('list_organisms_verify.csv', 'w') as file_csv:
        wr = csv.writer(file_csv, quoting=csv.QUOTE_ALL)
        wr.writerow(list_organism)


file_name = 'aaaaNB_ZERO_108_ML.csv'
list_couple_id = read_csv(file_name)
print(list_couple_id)

list_couples_obj = Couple.get_couples_by_list_id(list_couple_id)

list_organisms = []
for couple in list_couples_obj:
    if couple.fk_bacteria not in list_organisms:
        list_organisms.append(couple.fk_bacteria)
    if couple.fk_phage not in list_organisms:
        list_organisms.append(couple.fk_phage)

print(len(list_organisms))
write_ids_organism(list_organisms)

