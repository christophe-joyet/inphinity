# File : get_features.py
# Author : Christophe Joyet
# Date : Juin 2019

import csv
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from objects_API.BacteriophageJ import BacteriophageJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def get_features(csv_file:str, matrix_of_features:list, organisms_designation:list):
    """
    Read a csv file who contained features and fill matrix_of_features list with them
    Get organisms label from a csv and add them into organisms_designation

    :param csv_file: csv file with features
    :param matrix_of_features: matrix of features
    :param organisms_designation: list of organisms' designation

    :type csv_file: str
    :type matrix_of_features: list
    :type organisms_designation: list

    :return: None
    """
    # /!\ If a row description is not in the csv_file, add it /!\ 
    with open(csv_file, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            # Get the indexes of columns of features on the first iteration
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
            if 'MEAN_CE_H' in row:
                mean_ce_h = row.index('MEAN_CE_H')
            if 'MEAN_CE_N' in row:
                mean_ce_n = row.index('MEAN_CE_N')
            if 'MEAN_CE_O' in row:
                mean_ce_o = row.index('MEAN_CE_O')
            if 'MEAN_CE_C' in row:
                mean_ce_c = row.index('MEAN_CE_C')
            if 'MEAN_CE_S' in row:
                mean_ce_s = row.index('MEAN_CE_S')
            if 'WEIGHT (Da)' in row:
                weight = row.index('WEIGHT (Da)')
            if 'ISO POINT' in row:
                iso_point = row.index('ISO POINT')
            if 'AROMATICITY' in row:
                aromaticity = row.index('AROMATICITY')
                continue # On the first iteration we are on the first row and there is no data only labels

            # Handle empty case for X
            if row[mean_aa_x] == "":
                row[mean_aa_x] = 0.0

            #the number of features is the dimension of our final matrice (29 features)
            features = [float(row[mean_aa_a]), float(row[mean_aa_c]), float(row[mean_aa_d]), float(row[mean_aa_e]),
                        float(row[mean_aa_f]), float(row[mean_aa_g]), float(row[mean_aa_y]), float(row[mean_aa_h]),
                        float(row[mean_aa_i]), float(row[mean_aa_k]), float(row[mean_aa_l]), float(row[mean_aa_m]),
                        float(row[mean_aa_n]), float(row[mean_aa_p]), float(row[mean_aa_q]), float(row[mean_aa_r]),
                        float(row[mean_aa_s]), float(row[mean_aa_t]), float(row[mean_aa_v]), float(row[mean_aa_w]),
                        float(row[mean_ce_h]), float(row[mean_ce_n]), float(row[mean_ce_o]), float(row[mean_ce_c]),
                        float(row[mean_aa_x]), float(row[mean_ce_s]), float(row[weight]), float(row[iso_point]),
                        float(row[aromaticity])]

            matrix_of_features.append(features) # Add the features in the matrix
            organisms_designation.append(row[0]) # Add organisms designation
    f.close()