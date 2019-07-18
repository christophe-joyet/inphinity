#!/usr/bin/env python
# encoding: utf-8

import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage

def DisplayDendrogramFromMatrix(matrix_path:str, graph:str):
    """
    Display graphic from csv file

    :param matrix_path: path to the matrix

    :type matrix_path: str   
    """

    df = pd.read_csv(matrix_path)  # Data set

    df = df.set_index('Phages designation')  # delete the row with the labels
    del df.index.name
    df

    Z = hierarchy.linkage(df, 'ward')  # Calculate the distance between each sample
    
    if graph == "hierarchique":
        hierarchy.dendrogram(Z, color_threshold=1250000, orientation="right", leaf_rotation=0, leaf_font_size=8, labels=df.index)  # Plot with Custom leaves
    elif graph == "heatmap":
        sns.heatmap(df, cmap="RdYlGn", square=True, xticklabels=True, yticklabels=True)
    elif graph == "clustermap":
        sns.set(font_scale=0.5)
        sns.clustermap(df, metric="euclidean", standard_scale=1, method="ward", yticklabels=True, xticklabels=True, annot=True)
    
    # Display grapphic
    plt.show()