from pathlib import Path

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from objects_API.ProteinPFAMJ import ProteinPFAMJson
from objects_API.ProteinJ import ProteinJson
from objects_API.CoupleJ import CoupleJson

import pandas as pd
import numpy as np 
import itertools


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()


def getAllCouplesByParameters(parameter_dic:dict):
    list_couples = CoupleJson.getCouplesByFilterParameter(parameter_dic)
    return list_couples

def getAllProteinsByOrganism(organism_id:int):
    dict_params_protein = {}
    dict_params_protein['organism_id'] = organism_id
    list_proteins = ProteinJson.getProteinByFilterParameter(dict_params_protein)
    return list_proteins

def countChemicalComponents(proteic_sequence:str):
    is_aminoacid = False
    perc_charact = {}
    perc_charact['C_ch'] = 0
    perc_charact['H_ch'] = 0
    perc_charact['O_ch'] = 0
    perc_charact['N_ch'] = 0
    perc_charact['S_ch'] = 0


    for aminoacid in proteic_sequence:
        if aminoacid == 'A':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 3
            perc_charact['H_ch'] = perc_charact['H_ch'] + 7
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'C':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 3
            perc_charact['H_ch'] = perc_charact['H_ch'] + 7
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 1
        if aminoacid == 'D':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 4
            perc_charact['H_ch'] = perc_charact['H_ch'] + 6
            perc_charact['O_ch'] = perc_charact['O_ch'] + 4
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'E':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 5
            perc_charact['H_ch'] = perc_charact['H_ch'] + 8
            perc_charact['O_ch'] = perc_charact['O_ch'] + 4
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'F':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 3
            perc_charact['H_ch'] = perc_charact['H_ch'] + 6
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'G':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 2
            perc_charact['H_ch'] = perc_charact['H_ch'] + 5
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'H':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 6
            perc_charact['H_ch'] = perc_charact['H_ch'] + 9
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 3
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'I':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 6
            perc_charact['H_ch'] = perc_charact['H_ch'] + 13
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'K':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 6
            perc_charact['H_ch'] = perc_charact['H_ch'] + 15
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 2
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'L':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 6
            perc_charact['H_ch'] = perc_charact['H_ch'] + 13
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'M':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 5
            perc_charact['H_ch'] = perc_charact['H_ch'] + 11
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 1
        if aminoacid == 'N':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 4
            perc_charact['H_ch'] = perc_charact['H_ch'] + 8
            perc_charact['O_ch'] = perc_charact['O_ch'] + 3
            perc_charact['N_ch'] = perc_charact['N_ch'] + 2
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'P':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 5
            perc_charact['H_ch'] = perc_charact['H_ch'] + 9
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'Q':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 5
            perc_charact['H_ch'] = perc_charact['H_ch'] + 10
            perc_charact['O_ch'] = perc_charact['O_ch'] + 3
            perc_charact['N_ch'] = perc_charact['N_ch'] + 2
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'R':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 6
            perc_charact['H_ch'] = perc_charact['H_ch'] + 15
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 4
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'S':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 3
            perc_charact['H_ch'] = perc_charact['H_ch'] + 7
            perc_charact['O_ch'] = perc_charact['O_ch'] + 3
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'T':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 4
            perc_charact['H_ch'] = perc_charact['H_ch'] + 9
            perc_charact['O_ch'] = perc_charact['O_ch'] + 3
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'V':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 5
            perc_charact['H_ch'] = perc_charact['H_ch'] + 11
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'W':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 5
            perc_charact['H_ch'] = perc_charact['H_ch'] + 8
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 2
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0
        if aminoacid == 'Y':
            perc_charact['C_ch'] = perc_charact['C_ch'] + 3
            perc_charact['H_ch'] = perc_charact['H_ch'] + 7
            perc_charact['O_ch'] = perc_charact['O_ch'] + 2
            perc_charact['N_ch'] = perc_charact['N_ch'] + 1
            perc_charact['S_ch'] = perc_charact['S_ch'] + 0

    results_ch_percent = {}
    nb_chemical_componenents = sum(perc_charact.values())
    for key, value in perc_charact.items():
        results_ch_percent[key] = value / nb_chemical_componenents

    dataframe_chemica_components = pd.DataFrame.from_dict([results_ch_percent])
    return dataframe_chemica_components

def cartesian(df1:pd.DataFrame, df2:pd.DataFrame):
    """
    perform a cartesian product between two dataframe

    :param df1: dataframe 1
    :param df2: dataframe 2

    :type df1: dataframe 
    :type df2: dataframe

    :return: Dataframe with the cartesian product of both dataframes
    :rtype: pd.DataFrame

    """

    df_final = df1.assign(foo=1).merge(df2.assign(foo=1)).drop('foo', 1)
    return df_final
    #rows = itertools.product(df1.iterrows(), df2.iterrows())

    #df = pd.DataFrame(left.append(right) for (_, left), (_, right) in rows)
    #return df.reset_index(drop=True)

def computeMeanStd(dataframe_data:pd.DataFrame):
    """
    perform the mean and std of a given Dataframe, concatenate the data and return a new DataFrame with the results

    :param dataframe_data: dataframe with the results

    :type dataframe_data: dataframe 

    :return: Dataframe with the mean and STD of all columns
    :rtype: pd.DataFrame

    """

    dataframe_results = pd.DataFrame()
    results_mean = dataframe_data.mean().to_frame().T
    results_mean.columns = [str(col) + '_mean' for col in results_mean.columns]

    results_std = dataframe_data.std().to_frame().T
    results_std.columns = [str(col) + '_std' for col in results_std.columns]

    dataframe_results = pd.concat([results_mean, results_std],  axis = 1, sort=False)
    return dataframe_results

def validateProteinSequence(proteic_sequence:str):
    if '*' in proteic_sequence:
        proteic_sequence = proteic_sequence.replace('*','')
    if '"' in proteic_sequence:
        proteic_sequence = proteic_sequence.replace('"','')
    if 'X' in proteic_sequence or 'x' in proteic_sequence:
        proteic_sequence = proteic_sequence.replace('X','')
    if 'J' in proteic_sequence:
        proteic_sequence = proteic_sequence.replace('J','')
    if 'B' in proteic_sequence:
        proteic_sequence = proteic_sequence.replace('B','')
    if 'Z' in proteic_sequence:
        proteic_sequence = proteic_sequence.replace('Z','')


    return proteic_sequence

def calculatePercentAAMolecularWeightByListProteins(list_proteins:list, is_bacterium:bool):
    """

    calculate the percentage of all AA and these molecular weight for a list of proteins. Complete the columns name according a bacterium or phage and return a dataframe

    :param list_proteins: list of proteins
    :param is_bacterium: True if it is a list bacterium proteins

    :type list_proteins: list[ProteinJ] 
    :type is_bacterium: bool

    :Note pay attention that in this version the X sequence WAS NOT considered at all

    :return: Dataframe with the cartesian product of both dataframes
    :rtype: pd.DataFrame

    """

    dataframe_percents = pd.DataFrame()
    index_df = 0
    assert len(list_proteins) >= 3, "Are you sure that you may have proteins in this organism? "
    for protein in list_proteins:
        #print(protein.id)

        #prot_seq = protein.sequence_AA + 'XXXXXXXXXXXXXXXXXXXX'
        #sequence = Seq(prot_seq, IUPAC.extended_protein)
        #X.sequence.alphabet = IUPACData.protein_letters
        protein_sequence = protein.sequence_AA

        protein_sequence_molecular_weight = validateProteinSequence(protein_sequence)
        protein_sequence_molecular_weight = protein_sequence_molecular_weight.upper()
        prot_seq_treated = ProteinAnalysis(protein_sequence_molecular_weight)
        molecular_weight = prot_seq_treated.molecular_weight()
        chemical_component_dataframe = countChemicalComponents(protein_sequence)



        X = ProteinAnalysis(protein_sequence)
        dict_percent = X.get_amino_acids_percent()
        dict_percent['id_prot'] = protein.id

        dataframe_aux = pd.DataFrame.from_dict([dict_percent])
        dataframe_aux['moelcular_weight'] = molecular_weight
        dataframe_aux = pd.concat([dataframe_aux,chemical_component_dataframe], axis = 1, sort = False)

        dataframe_percents = dataframe_percents.append(dataframe_aux, ignore_index=True)

        index_df += 1
    dataframe_percents = dataframe_percents.set_index('id_prot')
    if is_bacterium == True:
        dataframe_percents.columns = [str(col) + '_bact' for col in dataframe_percents.columns]
    else:
        dataframe_percents.columns = [str(col) + '_phage' for col in dataframe_percents.columns]
    return dataframe_percents

def write_indexs_couples(list_index:list, path_save:str):
    """

    write the index of the couples already analysed 

    :param list_index: list of couples already treated
    :param path_save: path where you want to save the information

    :type list_index: list[ProteinJ] 
    :type path_save: str

    """
    np.savetxt(path_save, list_index, delimiter=',') 

def load_index_couple(path_load:str):
    """

    read the index of the couples treated

    :param path_load: path where you want to save the information

    :type path_load: str

    :return: array with the indexs
    :rtype: numPy array

    """

    array_index = []
    my_file = Path(path_load)
    if my_file.is_file():
        array_index = np.loadtxt(path_load, delimiter=',').astype(int)
    return array_index

def load_dataset_couples(path_load:str):
    """

    read the dataset already existant

    :param path_load: path where you want to save the information

    :type path_load: str

    :return: the dataframe that contain the dataset
    :rtype: pandas dataframe

    """

    dataframe_csv_existant = pd.DataFrame()
    my_file = Path(path_load)
    if my_file.is_file():
        dataframe_csv_existant = pd.read_csv(path_load) 
    return dataframe_csv_existant

dict_param_couple = {}
dict_param_couple['level_id'] = 1
list_couples = getAllCouplesByParameters(dict_param_couple)

id_old_bacterium = -1
list_prots_phage = []
list_prots_bact = []
dataframe_percents_bacterium = pd.DataFrame()


path_or_buffile_to_save = 'dataset_CH.csv'
path_index_couples_treated = 'index_couples.csv'
path_dataset_existant = 'dataset_CH.csv'

dataframe_results_CH =load_dataset_couples(path_dataset_existant)


array_id_couple_treatment = load_index_couple(path_index_couples_treated)


for couple_obj in list_couples:
    print('start couple {0}'.format(str(couple_obj.id)))
    id_couple = couple_obj.id

    if id_couple not in array_id_couple_treatment:
        id_new_phage = couple_obj.bacteriophage
        id_new_bacterium = couple_obj.bacterium


        list_prots_phage = getAllProteinsByOrganism(couple_obj.bacteriophage)
        dataframe_percents_bacteriophage = calculatePercentAAMolecularWeightByListProteins(list_prots_phage, False)

        if id_new_bacterium != id_old_bacterium:
            id_old_bacterium = id_new_bacterium
            list_prots_bact = getAllProteinsByOrganism(couple_obj.bacterium)
            dataframe_percents_bacterium = calculatePercentAAMolecularWeightByListProteins(list_prots_bact, True)

        dataFrame_Cartezian = cartesian(dataframe_percents_bacterium, dataframe_percents_bacteriophage)
        dataframe_resume_mean_std = computeMeanStd(dataFrame_Cartezian)

        dataframe_resume_mean_std['id_couple'] = couple_obj.id
        dataframe_resume_mean_std['label'] = couple_obj.interaction_type

        dataframe_results_CH = dataframe_results_CH.append(dataframe_resume_mean_std, ignore_index=True)

        array_id_couple_treatment = np.append(array_id_couple_treatment, id_couple)

        write_indexs_couples(array_id_couple_treatment, path_index_couples_treated)
        
        
        dataframe_results_CH.to_csv(path_or_buffile_to_save, index=False)
    print('End couple {0}'.format(str(couple_obj.id)))


dataframe_results_CH.to_csv(path_or_buffile_to_save, index=False)

list_prots = getAllProteinsByOrganism(list_couples[0].bacteriophage)

datafram_percents_bacteriophage = calculatePercentAAMolecularWeightByListProteins(list_prots, False)

list_prots = getAllProteinsByOrganism(list_couples[0].bacterium)
datafram_percents_bacterium = calculatePercentAAByListProteins(list_prots, True)


dataFrame_Cartezian = cartesian(datafram_percents_bacterium, datafram_percents_bacteriophage)
print(dataFrame_Cartezian)

print(datafram_percents_bacteriophage.shape)
print(datafram_percents_bacterium.shape)
print(dataFrame_Cartezian.shape)

computeMeanStd(dataFrame_Cartezian)

print(datafram_percents_bacterium)
print('Hello')
