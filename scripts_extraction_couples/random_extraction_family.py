#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, './')

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

import matplotlib.pyplot as plt
import csv
import numpy as np

from random import randint

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

negatif = 3250 # TOCOMPLETE
positif = 721 # TOCOMPLETE
files = int(negatif/positif) #change if negatif < positif
tab = np.zeros(files)
flag = 0
negative = 721 #TOCOMPLETE
range_min = 1
range_max = files + 1
t = 0
with open('/home/christophe/Documents/Phages4A/Family_data_split/greg/CH_dataset_greg_couples_2.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')

    for row in csv_reader:
        first = 1
        isnotok = 1

        #transform row in str
        row_str = ""
        for el in row:
            if first == 1:
                row_str += el
                first = 0
            else:
                row_str += ',' + el
        row_str += '\n'
        
        #first line
        if flag == 0:
            for i in range(range_min,range_max):
                with open('/home/christophe/Documents/Phages4A/Family_data_split/greg/CH_dataset_greg_F' + str(i) + '.csv', mode='a') as csv_entero:
                            csv_entero.write(row_str)
            flag = 1

        # si négatif < que positif changer '1' -> '0'
        if row[-1] == '1':
            for i in range(range_min, range_max):
                with open('/home/christophe/Documents/Phages4A/Family_data_split/greg/CH_dataset_greg_F' + str(i) + '.csv', mode='a') as csv_entero:
                            csv_entero.write(row_str)
            t += 1
            print(t)

        # si négatif < que positif changer '0' -> '1'
        elif row[-1] == '0':
            while isnotok:
                random_number = randint(range_min,range_max - 1)
                if tab[random_number - 1] < negative:
                    with open('/home/christophe/Documents/Phages4A/Family_data_split/greg/CH_dataset_greg_F' + str(random_number) + '.csv', mode='a') as csv_entero:
                        csv_entero.write(row_str)
                        tab[random_number - 1] += 1
                        isnotok = 0
        
            
