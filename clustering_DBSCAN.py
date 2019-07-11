# File : clustering_DBSCAN.py
# Author : Christophe Joyet
# Date : April 2019
# Modification : June 2019 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import read_csv as rcsv
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA #Features Exctraction

# =================================================================================================
# /!\ WARNING : IF COLUMN "MEAN_AA_X" DOESN'T EXIST IN FILE -> ADD IT /!\
# take data from csv file
csv_file = '../../extraction_features/AA_CE_WEIGHT_ARO_ISO/features_70_phages_clear_lysis.csv'

coordinates = []
phage_designation = []

rcsv.read_csv_with_mean_features(csv_file, coordinates, phage_designation)

# ==================================================================
# clustering
# ==================================================================

#reduction of all features
pca = PCA(n_components=2)
coordinates = pca.fit_transform(coordinates)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(coordinates)
# cluster the data
dbscan = DBSCAN(eps=0.250, min_samples = 2)
clusters = dbscan.fit_predict(X_scaled)
# plot the cluster assignments
plt.scatter(coordinates[:, 0], coordinates[:, 1], c=clusters, cmap="tab20")
plt.xticks([])
plt.yticks([])
plt.title("Graph - DBSCAN")

# ==================================================================
# display phages'name
# ==================================================================
for i, txt in enumerate(phage_designation):
    plt.annotate(txt, (coordinates[i,0],coordinates[i,1]))

plt.show()