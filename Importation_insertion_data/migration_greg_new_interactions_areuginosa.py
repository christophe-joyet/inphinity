import pandas as pd
import csv
import math

from objects_API.BacteriumJ import BacteriumJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.StrainJ import StrainJson
from objects_API.SpecieJ import SpecieJson
from objects_API.CoupleJ import CoupleJson

from objects_new.Couples_new import Couple
from objects_new.Organisms_new import Organism
from objects_new.Strains_new import Strain
from objects_new.Species_new import Specie

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI



def getIdOrganismContainName(datafram_data, name_organism:str):
    value = datafram_data[datafram_data['designation'].str.contains(name_organism)]

    quantity_results = value.shape[0]
    if quantity_results == 1:
        str_write = name_organism + ',' + str(value['organism_id'].item()) + '\n'
        with open(file_results_ids,'a') as fd:
            fd.write(str_write)
    else:
        str_write = name_organism + '\n'
        with open(file_results_ids_error,'a') as fd:
            fd.write(str_write)


def getDictionaryFromCSV(path_csv:str):
    reader = csv.reader(open(path_csv, 'r'))
    dict_data = {}
    for row in reader:
       k, v = row
       dict_data[k] = v
    return dict_data


def addInteractionsNewDB(interaction_type:bool, bacterium_id:int, phage_id:int, level_id:int, lysis_id:int, persone_responsible:int, source_data_id:int, validity_id:bool):
    couple_obj_json = None 
    couple_obj = CoupleJson(interaction_type = interaction_type, bacteriophage = phage_id, bacterium = bacterium_id, level = level_id, lysis = lysis_id, person_responsible = persone_responsible, source_data = source_data_id, validity = validity_id)

    if lysis_id != None:
        couple_obj_json = couple_obj.setCoupleWithLysis()
    else:
        couple_obj_json = couple_obj.setCouple()
    return couple_obj_json



def getBacteriumStrainSpecieDesignationById(id_bacterium):


    bacterium_json_obj = BacteriumJson.getByID(id_bacterium)
    strain_id = bacterium_json_obj.strain
    strain_obj = StrainJson.getByID(strain_id)
    strain_designation = strain_obj.designation

    specie_id = strain_obj.specie
    specie_obj = SpecieJson.getByID(specie_id)
    specie_designation = specie_obj.designation

    taxonomy_bacterium = 'Specie: ' + specie_designation + ' Strain: ' + strain_designation

    return taxonomy_bacterium

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def validateDataRow(list_phages_names:list, bacterium_id:int, dict_phages_id:dict, dict_types_interactions_id:dict, row_dataframe):
    couple_interaction_type = True
    type_lysis = None
    for phage_name in list_phages_names:
        phage_id = dict_phages_id[str(phage_name)]
        interaction_type = row_dataframe[phage_name]
        print(type(interaction_type))
        if RepresentsInt(interaction_type) == False:
            couple_interaction_type = False
            type_lysis = None
            print(couple_interaction_type)
        elif interaction_type in dict_types_interactions_id:
            couple_interaction_type = True
            type_lysis = dict_types_interactions_id[interaction_type]

        addInteractionsNewDB(couple_interaction_type, bacterium_id, phage_id, 4, type_lysis, persone_responsible = 5, source_data_id = 3, validity_id = 1)

        print(phage_id)
    print('Hello')

def getBacteriumDesignationNewByID(id_new_bacterium:int):
    print(type(id_new_bacterium))
    bacterium_obj = BacteriumJson.getByID(id_new_bacterium)
    strain_id = bacterium_obj.strain
    strain_obj = StrainJson.getByID(strain_id)
    strain_designation = strain_obj.designation

    specie_id = strain_obj.specie
    specie_obj = SpecieJson.getByID(specie_id)
    specie_designation = specie_obj.designation

    taxonomy_bacterium = 'Specie: ' + specie_designation + ' Strain: ' + strain_designation

    return taxonomy_bacterium

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()


path_excel_file_interaction = 'data_csv_import_greg/e_coli_interaction_rest.xlsx'
path_excel_file_bacterium_list = 'bacteria_greg.csv'

file_bacterium_ids = 'data_csv_import_greg/e_coli_ids.csv'
dict_ids_bacterium = getDictionaryFromCSV(file_bacterium_ids)

file_phages_ids = 'data_csv_import_greg/e_coli_ids_phage.csv'
dict_ids_bacteriophage = getDictionaryFromCSV(file_phages_ids)



dataframe_file_interaction = pd.read_excel(path_excel_file_interaction)
#dataframe_file_bacterium_list = pd.read_csv(path_excel_file_bacterium_list)

print('Hello?')



dict_types_interactions = {}
dict_types_interactions[1] = 8
dict_types_interactions[2] = 9
dict_types_interactions[3] = 10
dict_types_interactions[4] = 11


list_heads = list(dataframe_file_interaction.head(0))
list_heads_phages = list_heads[1:]

for index, row in dataframe_file_interaction.iterrows():
    bacterium_name = row['bac/phage']

    bacterium_name = str(bacterium_name)

    id_bacterium = dict_ids_bacterium[bacterium_name]
    id_bacterium = int(id_bacterium)
    taxonomy_bacterium = getBacteriumDesignationNewByID(id_bacterium)
    print(taxonomy_bacterium)
    #ask the user if we can insert the data
    validateDataRow(list_heads_phages, id_bacterium, dict_ids_bacteriophage, dict_types_interactions, row)
    print(bacterium_name)
    #getIdOrganismContainName(dataframe_file_bacterium_list, bacterium_name)
