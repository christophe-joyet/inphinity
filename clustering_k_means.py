# File : clustering_k_means.py
# Author : Christophe Joyet
# Date : April 2019

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA #Features Exctraction

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from objects_API.BacteriophageJ import BacteriophageJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

coordinates = []
phage_designation = []

# ==================================================================
# take data from csv file
# ==================================================================

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
        #add phages designation
        phage_designation.append(BacteriophageJson.getByID(row[0]).designation)
f.close()

# =================================================================================================
# 3D Graph
# =================================================================================================
'''
#we need a 3D representation -> we use PCA to reduce dimension
pca = PCA(n_components=3)
coordinates = pca.fit_transform(coordinates)
coordinates = np.array(coordinates)
#n_cluster => nombre de clusters que l'on veut.
#n_init => nombre de fois que l'on effectue l'algorithme depuis le départ
#max_iter => nombre d'itération pour trouver le meilleur centre
kmeans = KMeans(n_clusters=3, init='random', n_init=10, max_iter=100)
kmeans.fit(coordinates)
y_kmeans = kmeans.predict(coordinates)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], c=y_kmeans, cmap='viridis')
centers = kmeans.cluster_centers_
ax.scatter(centers[:, 0], centers[:, 1], c='green', s=200, alpha=0.5)

#affiche les noms
for i, txt in enumerate(phage_designation):
    ax.text(coordinates[i,0],coordinates[i,1],coordinates[i,2],  phage_designation[i], size=6) 

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.set_title("Graph - K-Mean Algorithme with K = 3")

'''
# =================================================================================================
# 2D Graph
# =================================================================================================
# we need a 2D representation -> we use PCA to reduce dimension of our matrix
pca = PCA(n_components=2)
coordinates = pca.fit_transform(coordinates)

coordinates = np.array(coordinates)
#n_cluster => number of clusters we want.
#n_init    => number of execution of the algorithm
#max_iter  => number of iteration to find the center
kmeans = KMeans(n_clusters=4, init='random', n_init=10, max_iter=100)
kmeans.fit(coordinates)
y_kmeans = kmeans.predict(coordinates)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(coordinates[:, 0], coordinates[:, 1], c=y_kmeans, s=50, cmap='viridis')
centers = kmeans.cluster_centers_
plt.scatter(centers[:, 0], centers[:, 1], c='black', s=200, alpha=0.5)

#display phages' name
for i, txt in enumerate(phage_designation):
    plt.annotate(txt, (coordinates[i,0],coordinates[i,1]))

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_title("Graph - K-Mean Algorithme with K = 4")

plt.show()

# Source : 
# VANDERPLAS, Jake, N.D. In Depth: k-Means Clustering. 
# Python Data Science HandBook [en ligne]. 
# [Consulté le 10 Avril 2019]. Disponible à l'adresse : https://jakevdp.github.io/PythonDataScienceHandbook/05.11-k-means.html
