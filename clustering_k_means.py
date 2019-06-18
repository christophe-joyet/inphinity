# File : clustering_k_means.py
# Author : Christophe Joyet
# Date : April 2019
# Modification : June 2019 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import read_csv as rcsv
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA #Features Exctraction
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

# =================================================================================================
# /!\ WARNING : IF NO "MEAN_AA_X" IN FILE -> ADD IT /!\
# take data from csv file
csv_file = '../../statistiques/CSV/Features_from_phages_in_Pseudomonas_aeruginosa_Muco16.csv'

coordinates = []
phage_designation = []

rcsv.read_csv_with_mean_features(csv_file, coordinates, phage_designation)

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
    ax.text(coordinates[i,0],coordinates[i,1],coordinates[i,2],  phage_designation[i], size=12) 

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
kmeans = KMeans(n_clusters=3, init='random', n_init=10, max_iter=100)
kmeans.fit(coordinates)
y_kmeans = kmeans.predict(coordinates)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(coordinates[:, 0], coordinates[:, 1], c=y_kmeans, s=50, cmap='viridis')
centers = kmeans.cluster_centers_
plt.scatter(centers[:, 0], centers[:, 1], c='black', s=200, alpha=0.5)

#display phages' name
for i, txt in enumerate(phage_designation):
    plt.annotate(txt, (coordinates[i,0],coordinates[i,1]), size=12)

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_title("Graph - K-Mean Algorithme with K = 3")


plt.show()
# Source : 
# VANDERPLAS, Jake, N.D. In Depth: k-Means Clustering. 
# Python Data Science HandBook [en ligne]. 
# [Consulté le 10 Avril 2019]. Disponible à l'adresse : https://jakevdp.github.io/PythonDataScienceHandbook/05.11-k-means.html
