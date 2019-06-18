# File : clustering_mean_shift.py
# Author : Christophe Joyet
# Date : April 2019

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import read_csv as rcsv
from sklearn.cluster import MeanShift
from sklearn.cluster import estimate_bandwidth
from sklearn.decomposition import PCA #Features Exctraction

# =================================================================================================
# /!\ WARNING : IF NO "MEAN_AA_X" IN FILE -> ADD IT /!\
# take data from csv file
csv_file = '../../statistiques/CSV/Features_from_phages_in_Pseudomonas_aeruginosa_Muco16.csv'

list_mean_coord_X = []
list_std_coord_Y = []
coordinates = []
phage_designation = []

rcsv.read_csv_with_mean_features(csv_file, coordinates, phage_designation)

# =================================================================================================
# Compute clustering with MeanShift
# =================================================================================================

# we need a 3D representation -> we use PCA to reduce dimension
pca = PCA(n_components=2)
coordinates = pca.fit_transform(coordinates)
coordinates = np.array(coordinates)

# The following bandwidth can be automatically detected using
bandwidth = estimate_bandwidth(coordinates, quantile=0.5, n_samples=500)

ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)

ms.fit(coordinates)
labels = ms.labels_
cluster_centers = ms.cluster_centers_

labels_unique = np.unique(labels)
n_clusters_ = len(labels_unique)

print("number of estimated clusters : %d" % n_clusters_)

# =================================================================================================
# Plot result
# =================================================================================================
import matplotlib.pyplot as plt
from itertools import cycle

fig = plt.figure()
ax = fig.add_subplot(111)

plt.clf()

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(coordinates[my_members, 0], coordinates[my_members, 1], col + '.')
    plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)
plt.xticks([])
plt.yticks([])
plt.title("Graph - K-MeanShift -  %d cluster(s) estimated" % n_clusters_)

# =================================================================================================
# display name of phages
# =================================================================================================

for i, txt in enumerate(phage_designation):
    plt.annotate(txt, (coordinates[:, 0][i], coordinates[:, 1][i]), size = 12)

plt.show()

# Source : 
# 2018. A demo of the mean-shift clustering algorithm. 
# scikit-learn [en ligne]. 
# [Consulté le 27 Avril 2019]. Disponible à l'adresse : 
# https://scikit-learn.org/stable/auto_examples/cluster/plot_mean_shift.html#sphx-glr-auto-examples-cluster-plot-mean-shift-py
