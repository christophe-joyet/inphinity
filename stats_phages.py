
# File : stats_phages.py
# Author : Christophe Joyet
# Date : Avril 2019

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import os
import network_chart_couples_cl as network
from collections import Counter
from collections import OrderedDict

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.ProteinJ import ProteinJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.OrganismJ import OrganismJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def getAllChemicalStructureOfAminoAcids(amino_acid:str):
    """
    get chemical strucutre of a given amino acid

    :param amino_acid: amino acid

    :type amino_acid: str

    :return: chemical strucutre of the given amino acid
    :rtype: dictionnary {atomium : occurence}
    """
    #source : https://www.britannica.com/science/amino-acid
    #dict contains 21 proteins with their molecules
    molecule_amino_acid = { 'A':{
                                'H' : 7,
                                'N' : 1,
                                'C' : 3,
                                'O' : 2},
                            'V':{
                                'H' : 11,
                                'N' : 1,
                                'C' : 5,
                                'O' : 2},
                            'L':{
                                'H' : 13,
                                'N' : 1,
                                'C' : 6,
                                'O' : 2},
                            'G':{
                                'H' : 5,
                                'N' : 1,
                                'C' : 2,
                                'O' : 2},
                            'I':{
                                'H' : 13,
                                'N' : 1,
                                'C' : 6,
                                'O' : 2},
                            'M':{
                                'H' : 11,
                                'N' : 1,
                                'C' : 5,
                                'O' : 2,
                                'S' : 1},
                            'W':{
                                'H' : 12,
                                'N' : 2 ,
                                'C' : 11,
                                'O' : 2},
                            'F':{
                                'H' : 11,
                                'N' : 1,
                                'C' : 9,
                                'O' : 2},
                            'P':{
                                'H' : 9,
                                'N' : 1,
                                'C' : 5,
                                'O' : 2},
                            'S':{
                                'H' : 7,
                                'N' : 1,
                                'C' : 3,
                                'O' : 3},
                            'C':{
                                'H' : 7,
                                'N' : 1,
                                'C' : 3,
                                'O' : 2,
                                'S' : 1},
                            'N':{
                                'H' : 8,
                                'N' : 2,
                                'C' : 4,
                                'O' : 3},
                            'Q':{
                                'H' : 10,
                                'N' : 2,
                                'C' : 5,
                                'O' : 3},
                            'T':{
                                'H' : 9,
                                'N' : 1,
                                'C' : 4,
                                'O' : 3},
                            'Y':{
                                'H' : 11,
                                'N' : 1,
                                'C' : 9,
                                'O' : 3}, 
                            'D':{
                                'H' : 7,
                                'N' : 1,
                                'C' : 4,
                                'O' : 4},
                            'E':{
                                'H' : 9,
                                'N' : 1,
                                'C' : 5,
                                'O' : 4},
                            'K':{
                                'H' : 14,
                                'N' : 2,
                                'C' : 6,
                                'O' : 2},
                            'R':{
                                'H' : 14,
                                'N' : 4,
                                'C' : 6,
                                'O' : 2},
                            'H':{
                                'H' : 9,
                                'N' : 3,
                                'C' : 6,
                                'O' : 2},
                            'U':{
                                'H' : 7,
                                'N' : 1,
                                'C' : 3,
                                'O' : 2,
                                'Se': 1}}  
    
    if amino_acid in molecule_amino_acid.keys():
        return molecule_amino_acid[amino_acid]
    else:
        return {'Unknow' : 1}
        
def getFrequenceAllAminoAcidForAProtein(protein:ProteinJson):
    """
    get the frequence of each amino acid of a given protein

    :param protein: protein object

    :type protein: ProteinJson

    :return: all the frequences of each amino acids of the given protein
    :rtype: dictionnary {name amino acid : frequence}
    """
    #dictionnary {amino acid : occurence}
    protein_amino_acid_dict = OrderedDict()
    #get the sequence with the amino acid of the protein
    sequence_of_amino_acid = protein.sequence_AA
    #count all the amino acid in the sequence
    for amino_acid in sequence_of_amino_acid:
        if not amino_acid in protein_amino_acid_dict.keys():
            protein_amino_acid_dict[amino_acid] = 1
        else:
            protein_amino_acid_dict[amino_acid] += 1

    #transform the occurence of each amino acids in frequence
    for amino_acid_key, amino_acid_value in protein_amino_acid_dict.items():
        #divide the occurence by the length of the sequence
        protein_amino_acid_dict[amino_acid_key] = amino_acid_value / len(protein.sequence_AA)

    return protein_amino_acid_dict

def getFeaturesForAPhage(phage:BacteriophageJson, active_percentage:bool=False, get_features:bool=False):
    """
    get features for a given phage

    :param phage: phage object
    :param active_percentage: provide result in percentage
    :param get_features: add mean and std
    
    :type phage: BacteriophageJson
    :type active_percentage: bool
    :type get_features: bool

    :return: features of the given phages (mean and std of 21 amino acids of the phage)
    :rtype: dictionnary {Amino Acid Description : features}
    """
 
    list_prot = ProteinJson.getByOrganismID(phage.id)
    phage_amino_acid_dict = OrderedDict()

    #get all amino acids of the phage
    for protein in list_prot:
        protein_amino_acids = getFrequenceAllAminoAcidForAProtein(protein)

        #add values from two dictionnaries in one
        # Source : 
        # MSEIFERT [Pseudonyme], 2017. You can use collections.Counter which implements addition + that way. 
        # Stackoverflow [en ligne]. 16 août 2017 à 2h29. 
        # [Consulté le 3 Avril 2019]. Disponible à l'adresse : https://stackoverflow.com/questions/45713887/add-values-from-two-dictionaries
        phage_amino_acid_dict = Counter(Counter(phage_amino_acid_dict) + Counter(protein_amino_acids))

    #calcul mean
    tot_protein_of_the_phage = len(ProteinJson.getByOrganismID(phage.id))
    dict_of_features = {}
    if get_features == True:
        #calcul mean of each amino acid
        for amino_acid_key in phage_amino_acid_dict.keys():
            dict_of_features['MEAN_AA_' + str(amino_acid_key)] = phage_amino_acid_dict[amino_acid_key] / tot_protein_of_the_phage
    
    return dict_of_features

def converstionValuesInPercentage(dictionnary:dict):
    """
    convert all value of a dictionnary in percent

    :param dictionnary: dictionnary

    :type dictionnary: dict

    :return: --
    :rtype: --
    """
    total_value = sum(dictionnary.values())
    for key, value in dictionnary.items():
        dictionnary[key] = (value * 100) / total_value
        
#==============================================================================
#====================================NOT WORKING===============================

def getAllChemicalStructureForAPhage(phage:BacteriophageJson, active_percentage:bool=False):
    """
    get all chemical structure for a phage

    :param phage: phage object
    :param active_percentage: provide result in percentage
    
    :type phage: BacteriophageJson
    :type active_percentage: bool

    :return: all the chemical strucure of the given phage
    :rtype: dictionnary {name amino acid : occurence}
    """
    phage_amino_acid_dict = getFeaturesForAPhage(phage)
    phage_chemical_struct_dict = OrderedDict()
    #get all amino acids of the phage
    for amino_acid, number_of_amino_acid in phage_amino_acid_dict.items():
        #get each molecule of the amino acid
        molecule_dict = getAllChemicalStructureOfAminoAcids(amino_acid)
        #add number of molecule to phage_chemical_struct_dict
        for molecule, number_of_molecules in molecule_dict.items():
            if not molecule in phage_chemical_struct_dict.keys():
                phage_chemical_struct_dict[molecule] = number_of_amino_acid * number_of_molecules
            else:
                phage_chemical_struct_dict[molecule] += number_of_amino_acid * number_of_molecules

    return phage_chemical_struct_dict

#==============================================================================
#==============================================================================
repository = '../../statistiques/CSV/'
file_name = 'fichier_test_5.csv'

CLEAR_LYSIS = 5
CLEAR_LYSIS_M1E7 = 8 
CLEAR_LYSIS_P1E7 = 9

#get the phages of all the couples with an interaction type
#list_couple = network.getCouplesInteraction(CoupleJson.getAllAPI(), interaction_type=True)
#get the phages of specific couples 
list_couple = network.getCouplesLysis([CLEAR_LYSIS])
list_phages = []
for couple in list_couple:
    if not couple.bacteriophage in list_phages:
        list_phages.append(couple.bacteriophage)

print("nombre de phages différents : " + str(len(list_phages)))

phages_amino_acids = OrderedDict()
for couple in list_couple:
    if not couple.bacteriophage in phages_amino_acids.keys():
        phage_amino_acid_dict = getFeaturesForAPhage(BacteriophageJson.getByID(couple.bacteriophage), active_percentage=False, get_features=True)
        phages_amino_acids[couple.bacteriophage] = phage_amino_acid_dict

# Source : 
# ADAM SMITH [Pseudonyme], 2015. This is beyond simple if you can use pandas. 
# Stackoverflow [en ligne]. 16 Juillet 2015 à 16h38. 
# [Consulté le 3 Avril 2019]. Disponible à l'adresse : https://stackoverflow.com/questions/31436783/writing-dictionary-of-dictionaries-to-csv-file-in-a-particular-format
df1 = pd.DataFrame.from_dict(data=phages_amino_acids, orient="index")
df1.to_csv(os.path.join(repository, file_name))

ending_message = "file " + file_name + " saved in " + repository
print(ending_message)

#==============================================================================
#==============================================================================