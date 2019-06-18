# File : clustering_t_SNE.py
# Author : Christophe Joyet
# Date : April 2019

from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import read_csv as rcsv

# =================================================================================================
# /!\ WARNING : IF NO "MEAN_AA_X" IN FILE -> ADD IT /!\
# take data from csv file
csv_file = '../../statistiques/CSV/Features_from_phages_in_Pseudomonas_aeruginosa_Muco16.csv'

coordinates = []
phage_designation = []

rcsv.read_csv_with_mean_features(csv_file, coordinates, phage_designation)

# =================================================================================================
# clustering
# =================================================================================================

X = np.array(coordinates)

#reduce the matrix
coordinates_embedded = TSNE(n_components=2, perplexity=5.0, learning_rate=500).fit_transform(coordinates)
y_data = coordinates_embedded[:,1]
plt.scatter(coordinates_embedded[:, 0], coordinates_embedded[:, 1], c=y_data)
plt.colorbar(orientation='horizontal', ticks=range(10))
plt.clim(0, 10)
plt.title("Graph - T-SNE")
plt.xticks([])
plt.yticks([])


# =================================================================================================
# display phages'name
# =================================================================================================
for i, txt in enumerate(phage_designation):
    plt.annotate(txt, (coordinates_embedded[i,0],coordinates_embedded[i,1]))
plt.show()
