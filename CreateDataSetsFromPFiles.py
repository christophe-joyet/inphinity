
import pickle
import os
import glob

import numpy as np
from numpy import genfromtxt
import csv
import re

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from objects_API.CoupleJ import CoupleJson

def getFilesFromPath(path_directory:str):
    """
    get all the files *.p in a give path. These files contains the ppi scores

    :param path_directory: path of the directory

    :type path_directory: str 

    :return: array with all the files name
    :rtype: array[str]

    """
    list_files = []
    os.chdir(path_directory)
    for file in glob.glob("*.p"):
        list_files.append(file)
    return list_files

def loadDictFromPIckle(path_file:str):
    """
    load the dictionary from a file

    :param path_file: path of the picke file

    :type path_file: str 

    :return: dictionary with the ddis sources ids
    :rtype: dict[(domain_id_a, domain_id_b)] : id sources

    """
    pickle_in = open(path_file,"rb")
    dict_ddi_source_ids = pickle.load(pickle_in)
    return dict_ddi_source_ids


def getMaxKeyValuePPI(dict_values_freq:dict):
    """
    get the maximum key value in a dictionary

    :param dict_values_freq: dictionary with the PPI scores

    :type dict_values_freq: dictionary[int:array] 

    :return: max value of PPI score
    :rtype: int

    """
    max_value_ppi = max(dict_values_freq, key=int)
    return max_value_ppi

def getMaxPPIValuesListOfFiles(name_files:list, path_directory:str):
    """
    get the maximum PPI score in all interactions files

    :param name_files: list with all the pickle interactions files
    :param path_directory: path of the directory where are the *.p files

    :type name_files: list[str]
    :type path_directory: str

    :return: max value found in all the PPI
    :rtype: int

    """
    max_score_found = -1
    for file_name in name_files:
        #complete_path = path_directory + '/' +  file_name
        complete_path = file_name
        dict_values = loadDictFromPIckle(complete_path)
        max_values = getMaxKeyValuePPI(dict_values)
        if max_values > max_score_found:
            max_score_found = max_values
    return max_score_found

def createArrayScoresPPI(dict_PPI_Scores:dict):
    """
    get an array with the PPI in scores in order to create the CSV

    :param dict_PPI_Scores: dict with the frequencies PPI scores

    :type dict_PPI_Scores: dict[int]:int

    :return: the array with the scores
    :rtype: array[int]

    """
    array_PPI_scores = []
    for key, value in dict_PPI_Scores.items():
        array_PPI_scores += value * [key]
    return array_PPI_scores

def createHistogram(array_scores_PPI:list, max_Score:int):
    """
    get an array with the histogram value for the given array_Score based on the max score and bin size = 1

    :note: you can find more information here: https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html

    :param array_scores_PPI: array with the scores frequency
    :param max_Score: max score PPI

    :type array_scores_PPI: array[int]
    :type max_Score: int

    :return: array with the frequency AND array with the bins
    :rtype: array[int], array[int]

    """
    max_Score = max_Score + 1
    bins_array = np.arange(max_Score)
    histo_array, bins_array = np.histogram(array_scores_PPI, bins=bins_array)
    return histo_array, bins_array

def getCoupleTyById(id_couple:int):
    """
    return a interaction type of a given couple ID

    :param id_couple: id of the couple

    :type id_couple: int

    :return: type of the couple
    :rtype: int

    """
    dict_parameters_couple = {}
    dict_parameters_couple['id'] = id_couple
    couple_obj = CoupleJson.getCouplesByFilterParameter(dict_parameters_couple)
    type_interaction = couple_obj.interaction_type
    return type_interaction

def createHistogramAllInteractions(list_files:list, path_directory:str):
    """
    create the list of arrays with the PPI scores for all the interactions

    :param list_files: list of the *.p files that contain the interactions
    :param path_directory: max score PPI

    :type array_scores_PPI: array[str,...]
    :type max_Score: int

    :return: List of array that contains the histograms
    :rtype: list[array]

    """

    p = re.compile(r'\d+')
    array2d_scores_interaction = []
    max_Score_found = getMaxPPIValuesListOfFiles(list_files, path_directory)
    complete_path = ''
    for list_p in list_files:
        #complete_path = path_directory + '/' +  list_p
        complete_path = list_p
        dict_ddis_ids = loadDictFromPIckle(complete_path)
        arrayScores = createArrayScoresPPI(dict_ddis_ids)
        array_freq, array_hist = createHistogram(arrayScores, max_Score_found)
        id_interaction = p.findall(list_p)
        assert len(id_interaction) == 1
        interaction_type = getCoupleTyById(id_interaction[0])
        array_freq = np.insert(array_freq, 0, id_interaction)
        array_freq = np.append(array_freq, interaction_type)

        array2d_scores_interaction.append(array_freq)
    return array2d_scores_interaction


def writeCSV(path_save_file:str, list_values_bins:list):
    """
    write CSV histogram data

    :param path_save_file: path where you want to save the dataset
    :param list_values_bins: list of arrays with the scores

    :type path_save_file: str
    :type list_values_bins: list[array]


    """
    with open(path_save_file, "w", newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(list_values_bins)


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

dict_ddis_ids = loadDictFromPIckle('files_data/couple_1000.p')
arrayScores = createArrayScoresPPI(dict_ddis_ids)
array_freq, array_hist = createHistogram(arrayScores, 579)

path_directory = 'files_data'
list_files = getFilesFromPath(path_directory)

max_Score_found = getMaxPPIValuesListOfFiles(list_files, path_directory)

array_scores = createHistogramAllInteractions(list_files, path_directory)
fil_csv_name ='bins_base.csv'
writeCSV(fil_csv_name, array_scores)
print('Hello')