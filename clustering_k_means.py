#import seaborn as sns; sns.set()  # for plot styling

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA #Features Exctraction

#source : https://matplotlib.org/gallery/mplot3d/scatter3d.html
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import


from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from objects_API.BacteriophageJ import BacteriophageJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

coordinates = []
phage_designation = []

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
        coordinates.append(couple)
        #ajout des designations des phages
        phage_designation.append(BacteriophageJson.getByID(row[0]).designation)
f.close()

#=========================================3D========================================================
#=========================================3D========================================================
'''
#we need a 3D representation -> we use PCA to reduce dimension
pca = PCA(n_components=3)
coordinates = pca.fit_transform(coordinates)
#source : https://jakevdp.github.io/PythonDataScienceHandbook/05.11-k-means.html
coordinates = np.array(coordinates)
#n_cluster => nombre de clusters que l'on veut.
#n_init => nombre de fois que l'on effectue l'algorithme depuis le départ
#max_iter => nombre d'itération pour trouver le meilleur centre
kmeans = KMeans(n_clusters=4, init='random', n_init=10, max_iter=100)
kmeans.fit(coordinates)
y_kmeans = kmeans.predict(coordinates)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], c=y_kmeans, cmap='viridis')
centers = kmeans.cluster_centers_
ax.scatter(centers[:, 0], centers[:, 1], c='green', s=200, alpha=0.5)

#affiche les noms
for i, txt in enumerate(phage_designation):
    ax.text(coordinates[i,0],coordinates[i,1],coordinates[i,2],  phage_designation[i], size=3) 

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
'''
#=========================================2D========================================================
#=========================================2D========================================================
#we need a 3D representation -> we use PCA to reduce dimension
pca = PCA(n_components=2)
coordinates = pca.fit_transform(coordinates)
#source : https://jakevdp.github.io/PythonDataScienceHandbook/05.11-k-means.html
coordinates = np.array(coordinates)
#n_cluster => nombre de clusters que l'on veut.
#n_init => nombre de fois que l'on effectue l'algorithme depuis le départ
#max_iter => nombre d'itération pour trouver le meilleur centre
kmeans = KMeans(n_clusters=4, init='random', n_init=10, max_iter=100)
kmeans.fit(coordinates)
y_kmeans = kmeans.predict(coordinates)

fig = plt.figure()
plt.scatter(coordinates[:, 0], coordinates[:, 1], c=y_kmeans, s=50, cmap='viridis')
centers = kmeans.cluster_centers_
plt.scatter(centers[:, 0], centers[:, 1], c='black', s=200, alpha=0.5)

#affiche les noms
for i, txt in enumerate(phage_designation):
    plt.annotate(txt, (coordinates[i,0],coordinates[i,1]))

plt.show()