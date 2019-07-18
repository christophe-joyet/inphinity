#!/usr/bin/env python
# encoding: utf-8

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv

from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA #Features Exctraction

import get_features_from_csv as gfcsv

csv_file = '../../extraction_features/AA_CE_WEIGHT_ARO_ISO/features_70_phages_clear_lysis.csv'
matrix_of_features = []
organisms_designation = []

# Fill the matrix with data and add phage designation
gfcsv.get_features(csv_file=csv_file,
                   matrix_of_features=matrix_of_features, 
                   organisms_designation=organisms_designation)
                   
pca = PCA(n_components=2)  # We want a 2D representation so we'll reduce our matrice to 2 components
matrix_of_features = pca.fit_transform(matrix_of_features)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(matrix_of_features)

dbscan_clustering = DBSCAN(eps=0.5, min_samples = 2)
clusters = dbscan_clustering.fit_predict(X_scaled)

# plot the cluster assignments
plt.scatter(matrix_of_features[:, 0], matrix_of_features[:, 1], c=clusters, cmap="tab20")
plt.xticks([])  # Remove axe X
plt.yticks([])  # Remove axe Y
plt.title('Graph - DBSCAN')

# Display phages' names
for i, txt in enumerate(organisms_designation):
    plt.annotate(txt, (matrix_of_features[i,0],matrix_of_features[i,1]))

plt.show()  # Display the graphic

# Source : 
# PIEROBON, Gabriel, 2018. DBSCAN clustering for data shapes k-means can’t handle well (in Python)
# Medium [en ligne]. 
# [Consulté le 10 Avril 2019]. 
# Disponible à l'adresse : 
# https://towardsdatascience.com/dbscan-clustering-for-data-shapes-k-means-cant-handle-well-in-python-6be89af4e6ea
