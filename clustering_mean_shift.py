# File : clustering_mean_shift.py
# Author : Christophe Joyet
# Date : April 2019
# Last modification : July 2019

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import get_features_from_csv as gfcsv
from itertools import cycle
from sklearn.cluster import MeanShift
from sklearn.cluster import estimate_bandwidth
from sklearn.decomposition import PCA # Using for features exctraction

csv_file = '../../extraction_features/AA_CE_WEIGHT_ARO_ISO/features_70_phages_clear_lysis.csv' # Take data from csv file
matrix_of_features = []
organisms_designation = []

gfcsv.get_features(csv_file=csv_file, 
                   matrix_of_features=matrix_of_features, 
                   organisms_designation=organisms_designation) # Fill the matrix with data and add phage designation
pca = PCA(n_components=2) # We want a 2D representation so we'll reduce our matrice to 2 components
matrix_of_features = pca.fit_transform(matrix_of_features) 
matrix_of_features = np.array(matrix_of_features) # Transform the matrix in a 2D array

# Compute clustering with MeanShift
bandwidth = estimate_bandwidth(matrix_of_features, quantile=0.15, n_samples=500) 
mean_shift_clustering = MeanShift(bandwidth=bandwidth, bin_seeding=True) 
mean_shift_clustering.fit(matrix_of_features)
labels = mean_shift_clustering.labels_ # Labels of each phage
cluster_centers = mean_shift_clustering.cluster_centers_

# Defining how much label we have and so the number of clusters
labels_unique = np.unique(labels)
n_clusters_ = len(labels_unique) 

fig = plt.figure()
ax = fig.add_subplot(111) # The first 1 = number of row, second 1 = number of columns, third 1 = index of subplot
plt.clf() # Clear current figure

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(matrix_of_features[my_members, 0], matrix_of_features[my_members, 1], col + '.')
    plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)
plt.xticks([]) # Remove axe X
plt.yticks([]) # Remove axe Y
plt.title("Mean Shift -  %d cluster(s)" % n_clusters_)

# Display phages' names
for i, txt in enumerate(organisms_designation):
    plt.annotate(txt, (matrix_of_features[:, 0][i], matrix_of_features[:, 1][i]), size = 12)

plt.show() # Display the graphic

# Source : 
# 2018. A demo of the mean-shift clustering algorithm. 
# scikit-learn [en ligne]. 
# [Consulté le 27 Avril 2019]. Disponible à l'adresse : 
# https://scikit-learn.org/stable/auto_examples/cluster/plot_mean_shift.html#sphx-glr-auto-examples-cluster-plot-mean-shift-py
