#import seaborn as sns; sns.set()  # for plot styling

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from sklearn.cluster import *
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from objects_API.BacteriophageJ import BacteriophageJson
from sklearn.decomposition import PCA #Features Exctraction


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()


list_mean_coord_X = []
list_std_coord_Y = []
X = []
phage_id = []

with open('../../statistiques/CSV/fichier_test.csv', 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        if 'P_AA_MEAN' in row:
            mean_index = row.index('P_AA_MEAN')
        if 'P_AA_STD' in row: 
            std_index = row.index('P_AA_STD')
        if 'P_AA_MED' in row:
            med_index = row.index('P_AA_MED')
        if 'P_AA_VAR' in row:
            var_index = row.index('P_AA_VAR')
        if 'P_AA_SKEWNESS' in row:
            skewness_index = row.index('P_AA_SKEWNESS')
        if 'P_AA_KURTOSIS' in row:
            kurt_index = row.index('P_AA_KURTOSIS')
            continue

        #the number of features is the dimension of our final matrice
        couple = [float(row[mean_index]),float(row[std_index]),float(row[med_index]),float(row[var_index]),float(row[skewness_index]),float(row[kurt_index])]
        X.append(couple)
        #ajout des id
        phage_id.append(BacteriophageJson.getByID(row[0]).designation)
f.close()

print(X)
#we need a 3D representation -> we use PCA to reduce dimension
pca = PCA(n_components=2)
X = pca.fit_transform(X)
print(X)
X = np.array(X)
#Source : https://scikit-learn.org/stable/auto_examples/cluster/plot_mean_shift.html#sphx-glr-auto-examples-cluster-plot-mean-shift-py
# #############################################################################
# Compute clustering with MeanShift

# The following bandwidth can be automatically detected using
bandwidth = estimate_bandwidth(X, quantile=0.2, n_samples=500)

ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
ms.fit(X)
labels = ms.labels_
cluster_centers = ms.cluster_centers_

labels_unique = np.unique(labels)
n_clusters_ = len(labels_unique)

print("number of estimated clusters : %d" % n_clusters_)

# #############################################################################
# Plot result
import matplotlib.pyplot as plt
from itertools import cycle

plt.figure(1)
plt.clf()

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(X[my_members, 0], X[my_members, 1], col + '.')
    plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)
plt.title('Estimated number of clusters: %d' % n_clusters_)
#afficher les noms
for i, txt in enumerate(phage_id):
    plt.annotate(txt, (X[:, 0][i], X[:, 1][i]))
plt.show()
