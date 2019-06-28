# File : dendogram.py
# Author : Christophe Joyet
# Date : Mai 2019

# Libraries
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage

import pandas as pd
import numpy as np
import seaborn as sns

def DisplayDendrogramFromMatrix(matrix_path:str, graph:str):
    """
    Display Dendrogram From csv file

    :param matrix_path: path to the matrix

    :type matrix_path: str   
    """

    # Data set
    df = pd.read_csv(matrix_path)

    # delete the row with the labels
    df = df.set_index('Phages designation')
    del df.index.name
    df

    # Calculate the distance between each sample
    Z = hierarchy.linkage(df, 'ward')
    
    if graph == "hierarchique":
        # Plot with Custom leaves
        hierarchy.dendrogram(Z, orientation="right", leaf_rotation=0, leaf_font_size=8, labels=df.index)
    elif graph == "heatmap":
        # Standardize or Normalize every column in the figure
        sns.heatmap(df, cmap="RdYlGn", square=True, xticklabels=True, yticklabels=True)
    elif graph == "clustermap":
        sns.clustermap(df, metric="euclidean", standard_scale=1, method="ward")

    # Display
    plt.show()
