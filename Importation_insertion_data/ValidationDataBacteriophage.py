from objects_API.BacteriophageJ import BacteriophageJson
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_new.Organisms_new import Organism
from objects_new.Proteins_new import Protein
import pandas as pd 
import numpy as np
import csv

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

file_write = 'qty_proteins_validation.csv'
data = pd.read_csv("./verification_bacteriophages.csv")

def writeCSV(file_name:str, id_bacteriophage_API:int, id_bacteriophage_DB:int, qty_prots_API:int, qty_prots_DB:int, status_validation:bool):

    array_values = np.array([id_bacteriophage_API, id_bacteriophage_DB, qty_prots_API, qty_prots_DB, status_validation])

    with open(file_name, 'a') as csvFile:
        writer = csv.writer(csvFile, delimiter = ',', lineterminator='\n')
        writer.writerow(array_values)
    csvFile.close()



def getProteinsAPI(id_bacteriophage):
    bacteriophage_json = BacteriophageJson.getByID(id_bacteriophage)
    qty_prots_API = len(bacteriophage_json.protein_organism)

    return qty_prots_API

def getProteinsDB(id_bacteriophage):
    bacteriophage = Organism.get_organism_by_id(id_bacteriophage)
    list_proteins = Protein.get_all_Proteins_by_organism_id(id_bacteriophage)
    qty_prots_DB = len(list_proteins)

    return qty_prots_DB



for index, row in data.iterrows():
    id_old_db = int(row['id_old_db'])
    id_api = int(row['id_api'])
    qty_prots_old_db = getProteinsDB(id_old_db)
    qty_prots_api = getProteinsAPI(id_api)
    validator = False
    print(id_api)
    if ((qty_prots_old_db == qty_prots_api) or 
        (qty_prots_old_db == (qty_prots_api+1)) or
        (qty_prots_old_db == (qty_prots_api-1))):
        validator = True
    writeCSV(file_write, id_api, id_old_db, qty_prots_api, qty_prots_old_db, validator)
print('Fini')