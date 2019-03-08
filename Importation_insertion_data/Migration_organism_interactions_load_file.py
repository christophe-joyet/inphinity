import pandas as pd
import numpy as np
from objects_API.BacteriumJ import BacteriumJson
from objects_API.BacteriophageJ import BacteriophageJson

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

def get_data_interactions_by_csv(path_file):
    """
    load the interactions according the path file


    :param path_file: file path name

    :type path_file: str - required   

    :return dataframe with the information
    :rtype pandas dataframe
    """
    dataframe_interaction = pd.read_csv(filepath_or_buffer = path_file)
    return dataframe_interaction


def get_id_bacterium_by_acc(bacterium_acc):
    """
    get the id of a bacterium according to its acc number


    :param bacterium_acc: acc number

    :type path_file: str - required   

    :return id of the bacterium
    :rtype int
    """

    bacterium_obj = BacteriumJson.getByAccnumber(bacterium_acc)
    bacterium_id = bacterium_obj.id
    return bacterium_id

def get_id_bacteriophage_by_designation(phage_designation:str):
    """
    get the id of a bacteriophage according to its designation

    :param phage_designation: designation of the bacteriophage

    :type phage_designation: str - required   

    :return id of the bacteriophage
    :rtype int
    """
    bacteriophage_obj = BacteriophageJson.getBydesignation(phage_designation)
    bacteriophage_id = bacteriophage_obj.id
    return bacteriophage_id


def get_dictionary_index_bacterium(dataframe_interactions):
    """
    get a dictionary with the ids of the bacteria in the dataframe

    :param dataframe_interactions: dataframe with all the data

    :type dataframe_interactions: dataframe - required   

    :return dictionary with the bacteriums ids
    :rtype dict{[acc_value]:id}
    """

    dataframe_interactions = dataframe_interactions.sort_values(by=['acc_bacterium'], ascending=[1])
    old_acc = ''
    dict_ids_bact = {}
    for index, row in dataframe_data.iterrows():
        bacterium_acc = row["acc_bacterium"]
        if bacterium_acc != old_acc:
            try:
                bacterium_id = get_id_bacterium_by_acc(bacterium_acc)
                old_acc = bacterium_acc
                dict_ids_bact[bacterium_acc] = bacterium_id
            except:
                dict_ids_bact[bacterium_acc] = -1

    return dict_ids_bact

def getIdsBacteriophagesByDataframe(dataframe_interactions):
    """
    get a dictionary with the ids of the bacteriophages in the dataframe

    :param dataframe_interactions: dataframe with all the data

    :type dataframe_interactions: dataframe - required   

    :return dictionary with the bacteriophage ids
    :rtype dict{[acc_value]:id}
    """
    dataframe_interactions = dataframe_interactions.sort_values(by=['phage_name'], ascending=[1])
    old_design = ''
    dict_designation_id = {}
    for index, row in dataframe_interactions.iterrows():
        bacteriophage_design = row["phage_name"]
        if bacteriophage_design != old_design:
            try:
                bacteriophage_id = get_id_bacteriophage_by_designation(bacteriophage_design)
                old_design = bacteriophage_design
                dict_designation_id[bacteriophage_design] = bacteriophage_id
            except:
                dict_designation_id[bacteriophage_design] = -1
    return dict_designation_id


def createIdFileInteractions(dict_bact, dict_phage, dataframe_data_interaction, path_file_name):
    array_information = []
    for index, row in dataframe_data_interaction.iterrows():
        id_bacterium = dict_bact[row['acc_bacterium']]
        id_phage = dict_phage[row['phage_name']]
        array_information.append([row['acc_bacterium'], row['phage_name'], id_bacterium, id_phage])
    datafram_data = pd.DataFrame(array_information, columns = ['acc_bact','design_phage', 'id_bact', 'id_phage'] )
    datafram_data.to_csv(path_or_buf = path_file_name)
    

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

path_csv_file = 'data_interactions.csv'

dataframe_data = get_data_interactions_by_csv(path_csv_file)


dict_ids_bacteria = get_dictionary_index_bacterium(dataframe_data)
dict_ids_bacteriophage = getIdsBacteriophagesByDataframe(dataframe_data)

path_save_file = 'data_ids_interaction.csv'
createIdFileInteractions(dict_ids_bacteria, dict_ids_bacteriophage, dataframe_data, path_save_file)




print(dict_designation_id)