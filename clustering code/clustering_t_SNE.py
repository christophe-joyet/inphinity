#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, './')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv

from sklearn.manifold import TSNE

from features import features_functions

csv_file = './features/CSV_files/features_70_phages_clear_lysis.csv'  # Take data from csv file
matrix_of_features = []
organisms_designation = []

# Fill the matrix with data and add phage designation
features_functions.createMatrixFeatures(csv_file=csv_file, 
                                        matrix_of_features=matrix_of_features, 
                                        organisms_designation=organisms_designation) 

matrix_of_features = np.array(matrix_of_features)
matrix_of_features_embedded = TSNE(n_components=2, perplexity=25.0, learning_rate=50).fit_transform(matrix_of_features)
y_data = matrix_of_features_embedded[:,1]
plt.scatter(matrix_of_features_embedded[:, 0], matrix_of_features_embedded[:, 1], c=y_data, cmap="tab10")
plt.colorbar(orientation='horizontal', ticks=range(10))
plt.clim(0, 10)
plt.title('Graph - T-SNE')

# Display organisms' names
for i, txt in enumerate(organisms_designation):
    plt.annotate(txt, (matrix_of_features_embedded[i,0],matrix_of_features_embedded[i,1]))

plt.show()  # Display the graphic
