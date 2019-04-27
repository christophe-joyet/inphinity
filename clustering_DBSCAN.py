# File : clustering_DBSCAN.py
# Author : Christophe Joyet
# Date : April 2019

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA #Features Exctraction


from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from objects_API.BacteriophageJ import BacteriophageJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

# ==================================================================
# take data from csv file
# ==================================================================

coordinates = []
phage_designation = []

with open('../../statistiques/CSV/phages_in_clear_lysis_couple.csv', 'r') as f:
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
        #ajout des designations des phages
        phage_designation.append(BacteriophageJson.getByID(row[0]).designation)
f.close()

#reduction of all features
pca = PCA(n_components=2)
coordinates = pca.fit_transform(coordinates)

# ==================================================================
# clustering
# ==================================================================

scaler = StandardScaler()
X_scaled = scaler.fit_transform(coordinates)
# cluster the data
dbscan = DBSCAN(eps=0.123, min_samples = 2)
clusters = dbscan.fit_predict(X_scaled)
# plot the cluster assignments
plt.scatter(coordinates[:, 0], coordinates[:, 1], c=clusters, cmap="plasma")
plt.xticks([])
plt.yticks([])
plt.title("Graph - DBSCAN")

# ==================================================================
# display phages'name
# ==================================================================
for i, txt in enumerate(phage_designation):
    plt.annotate(txt, (coordinates[i,0],coordinates[i,1]))

plt.show()