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
# /!\ WARNING : IF COLUMN "MEAN_AA_X" DOESN'T EXIST IN FILE -> ADD IT /!\
# take data from csv file
csv_file = '../../extraction_features/AA_CE_WEIGHT_ARO_ISO/features_70_phages_clear_lysis.csv'

coordinates = []
phage_designation = []

rcsv.read_csv_with_mean_features(csv_file, coordinates, phage_designation)

# =================================================================================================
# clustering
# =================================================================================================

X = np.array(coordinates)
#reduce the matrix
coordinates_embedded = TSNE(n_components=2, perplexity=25.0, learning_rate=50).fit_transform(coordinates)
y_data = coordinates_embedded[:,1]
plt.scatter(coordinates_embedded[:, 0], coordinates_embedded[:, 1], c=y_data, cmap="tab10")

plt.colorbar(orientation='horizontal', ticks=range(10))
plt.clim(0, 10)
plt.title("Graph - T-SNE")

# =================================================================================================
# display phages'name
# =================================================================================================
for i, txt in enumerate(phage_designation):
    plt.annotate(txt, (coordinates_embedded[i,0],coordinates_embedded[i,1]))
plt.show()
