from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv

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

X = np.array(coordinates)
print(X.shape)
#reduire notre matrice 
coordinates_embedded = TSNE(n_components=2).fit_transform(coordinates)
plt.scatter(coordinates_embedded[:, 0], coordinates_embedded[:, 1], label=phage_designation)
plt.legend()
plt.show()
