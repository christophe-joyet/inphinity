#!/usr/bin/env python
# encoding: utf-8

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA #Features Exctraction
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import get_features_from_csv as gfcsv

csv_file = '../../extraction_features/AA_CE_WEIGHT_ARO_ISO/features_70_phages_clear_lysis.csv'
matrix_of_features = []
organisms_designation = []
n_clusters = 9  # Number of clusters

# Fill the matrix with data and add phage designation
gfcsv.get_features(csv_file=csv_file, 
                   matrix_of_features=matrix_of_features, 
                   organisms_designation=organisms_designation) 

pca = PCA(n_components=2)  # We want a 2D representation so we'll reduce our matrice to 2 components
matrix_of_features = pca.fit_transform(matrix_of_features)
matrix_of_features = np.array(matrix_of_features)  # Transform the matrix in a 2D array

# n_init = Number of executions of the algorithm
# max_iter = Number of iterations to find the center
k_means_clustering = KMeans(n_clusters=n_clusters, init='random', n_init=10, max_iter=100)
k_means_clustering.fit(matrix_of_features)
labels = k_means_clustering.predict(matrix_of_features)  # Labels of each phage

fig = plt.figure()
ax = fig.add_subplot(111)  # The first 1 = number of row, second 1 = number of columns, third 1 = index of subplot
plt.clf()  # Clear current figure
plt.scatter(matrix_of_features[:, 0], matrix_of_features[:, 1], c=labels, s=50, cmap='viridis')
cluster_centers = k_means_clustering.cluster_centers_
plt.scatter(cluster_centers[:, 0], cluster_centers[:, 1], c='black', s=200, alpha=0.5)

# Display organisms' names
for i, txt in enumerate(organisms_designation):
    plt.annotate(txt, (matrix_of_features[i,0],matrix_of_features[i,1]), size=10)

ax.set_xticklabels([])  # Clean axe X
ax.set_yticklabels([])  # Clean axe Y
ax.set_title('K-Means Algorithm with K = ' + str(n_clusters))

# To have a 3D reprensatation of the data decomment
'''
pca = PCA(n_components=3)
matrix_of_features = pca.fit_transform(matrix_of_features)
matrix_of_features = np.array(matrix_of_features)
k_means_clustering = KMeans(n_clusters=3, init='random', n_init=10, max_iter=100)
k_means_clustering.fit(matrix_of_features)
labels = k_means_clustering.predict(matrix_of_features)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(matrix_of_features[:, 0], matrix_of_features[:, 1], matrix_of_features[:, 2], c=y_k_means_clustering, cmap='viridis')
centers = k_means_clustering.cluster_centers_
ax.scatter(centers[:, 0], centers[:, 1], c='green', s=200, alpha=0.5)

# Display phages' names
for i, txt in enumerate(organisms_designation):
    ax.text(matrix_of_features[i,0],matrix_of_features[i,1],matrix_of_features[i,2],  organisms_designation[i], size=12) 

ax.set_xticklabels([])  # Clean axe X
ax.set_yticklabels([])  # Clean axe Y
ax.set_zticklabels([])  # Clean axe Z
ax.set_title('Graph - K-Means Algorithme')
'''

plt.show()  # Display the graphic

# Source : 
# VANDERPLAS, Jake, N.D. In Depth: k-Means Clustering. 
# Python Data Science HandBook [en ligne]. 
# [Consulté le 10 Avril 2019]. Disponible à l'adresse : https://jakevdp.github.io/PythonDataScienceHandbook/05.11-k-means.html
