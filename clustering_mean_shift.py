# File : clustering_mean_shift.py
# Author : Christophe Joyet
# Date : April 2019

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from sklearn.cluster import *
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from objects_API.BacteriophageJ import BacteriophageJson
from sklearn.decomposition import PCA #Features Exctraction


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()


list_mean_coord_X = []
list_std_coord_Y = []
X = []
phage_id = []

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

        X.append(couple)
        #ajout des id
        phage_id.append(BacteriophageJson.getByID(row[0]).designation)
f.close()

# we need a 3D representation -> we use PCA to reduce dimension
pca = PCA(n_components=2)
X = pca.fit_transform(X)
X = np.array(X)

# ==================================================================
# Compute clustering with MeanShift
# ==================================================================

# The following bandwidth can be automatically detected using
bandwidth = estimate_bandwidth(X, quantile=0.5, n_samples=500)

ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)

ms.fit(X)
labels = ms.labels_
cluster_centers = ms.cluster_centers_

labels_unique = np.unique(labels)
n_clusters_ = len(labels_unique)

print("number of estimated clusters : %d" % n_clusters_)

# ==================================================================
# Plot result
# ==================================================================
import matplotlib.pyplot as plt
from itertools import cycle

fig = plt.figure()
ax = fig.add_subplot(111)

plt.clf()

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(X[my_members, 0], X[my_members, 1], col + '.')
    plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)
plt.xticks([])
plt.yticks([])
plt.title("Graph - K-MeanShift -  %d cluster(s) estimated" % n_clusters_)

# ==================================================================
# display name of phages
# ==================================================================

for i, txt in enumerate(phage_id):
    plt.annotate(txt, (X[:, 0][i], X[:, 1][i]), size = 12)

plt.show()

# Source : 
# 2018. A demo of the mean-shift clustering algorithm. 
# scikit-learn [en ligne]. 
# [Consulté le 27 Avril 2019]. Disponible à l'adresse : 
# https://scikit-learn.org/stable/auto_examples/cluster/plot_mean_shift.html#sphx-glr-auto-examples-cluster-plot-mean-shift-py
