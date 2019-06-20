
# File : read_csv.py
# Author : Christophe Joyet
# Date : Juin 2019

import csv

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from objects_API.BacteriophageJ import BacteriophageJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def read_csv_with_mean_features(csv_file, coordinates, phage_designation):
    """
    read csv file with mean features
    fill coordinates
    get phage designation from csv and add them into phage_designation

    :param csv_file: csv file with features
    :param coordinates: list of features
    :param phage_designation: phages' designation

    :type csv_file: str
    :type coordinates: list
    :type phage_designation: list

    :return: None
    """
    with open(csv_file, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            if 'MEAN_AA_M' in row:
                mean_aa_m = row.index('MEAN_AA_M')
            if 'MEAN_AA_E' in row:
                mean_aa_e = row.index('MEAN_AA_E')
            if 'MEAN_AA_N' in row:
                mean_aa_n = row.index('MEAN_AA_N')
            if 'MEAN_AA_Y' in row:
                mean_aa_y = row.index('MEAN_AA_Y')
            if 'MEAN_AA_K' in row:
                mean_aa_k = row.index('MEAN_AA_K')
            if 'MEAN_AA_F' in row:
                mean_aa_f = row.index('MEAN_AA_F')
            if 'MEAN_AA_I' in row:
                mean_aa_i = row.index('MEAN_AA_I')
            if 'MEAN_AA_A' in row:
                mean_aa_a = row.index('MEAN_AA_A')
            if 'MEAN_AA_H' in row:
                mean_aa_h = row.index('MEAN_AA_H')
            if 'MEAN_AA_L' in row:
                mean_aa_l = row.index('MEAN_AA_L')
            if 'MEAN_AA_V' in row:
                mean_aa_v = row.index('MEAN_AA_V')
            if 'MEAN_AA_Q' in row:
                mean_aa_q = row.index('MEAN_AA_Q')
            if 'MEAN_AA_R' in row:
                mean_aa_r = row.index('MEAN_AA_R')
            if 'MEAN_AA_T' in row:
                mean_aa_t = row.index('MEAN_AA_T')
            if 'MEAN_AA_D' in row:
                mean_aa_d = row.index('MEAN_AA_D')
            if 'MEAN_AA_S' in row:
                mean_aa_s = row.index('MEAN_AA_S')
            if 'MEAN_AA_W' in row:
                mean_aa_w = row.index('MEAN_AA_W')
            if 'MEAN_AA_C' in row:
                mean_aa_c = row.index('MEAN_AA_C')
            if 'MEAN_AA_G' in row:
                mean_aa_g = row.index('MEAN_AA_G')
            if 'MEAN_AA_P' in row:
                mean_aa_p = row.index('MEAN_AA_P')   
            if 'MEAN_AA_X' in row:
                mean_aa_x = row.index('MEAN_AA_X')
                continue

            #handle empty case for X
            if row[mean_aa_x] == "":
                row[mean_aa_x] = 0.0

            #the number of features is the dimension of our final matrice
            couple = [float(row[mean_aa_a]),float(row[mean_aa_c]),float(row[mean_aa_d]),float(row[mean_aa_e]),
                        float(row[mean_aa_f]),float(row[mean_aa_g]),float(row[mean_aa_y]),float(row[mean_aa_h]),
                        float(row[mean_aa_i]),float(row[mean_aa_k]),float(row[mean_aa_l]),float(row[mean_aa_m]),
                        float(row[mean_aa_n]),float(row[mean_aa_p]),float(row[mean_aa_q]),float(row[mean_aa_r]),
                        float(row[mean_aa_s]),float(row[mean_aa_t]),float(row[mean_aa_v]),float(row[mean_aa_w]),
                        float(row[mean_aa_x])]

            coordinates.append(couple)
            #add phages designation
            phage_designation.append(BacteriophageJson.getByID(row[0]).designation)
    f.close()