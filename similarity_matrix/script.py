#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, './')

from similarity_matrix import display_graph
from similarity_matrix import matrix_similarity_functions

# Create similatry matrix
matrix_similarity_functions.matrixSimilarityScript('Bacteries', './similarity_matrix/CSV_files', 5190)

# Display a heatmap of similatry matrix
display_graph.displayGraphFromMatrix('./similarity_matrix/CSV_files/Bacteries', 'heatmap')