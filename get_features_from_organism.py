
# File : stats_phages.py
# Author : Christophe Joyet
# Date : Avril 2019
# Last Modification : July 2019
# - get the weight of protein

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import os
from collections import Counter
from collections import OrderedDict

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.ProteinJ import ProteinJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.OrganismJ import OrganismJson
from Bio.SeqUtils.ProtParam import ProteinAnalysis

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()
        
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
    :rtype: dictionnary {Amino Acid Description : MEAN}
    """
 
    list_prot = ProteinJson.getByOrganismID(phage.id)
    phage_amino_acid_dict = OrderedDict()
    protein_weight = 0

    #get all amino acids of the phage
    for protein in list_prot:
        protein_amino_acids = getFrequenceAllAminoAcidForAProtein(protein)

        #add values from two dictionnaries in one
        # Source : 
        # MSEIFERT [Pseudonyme], 2017. You can use collections.Counter which implements addition + that way. 
        # Stackoverflow [en ligne]. 16 août 2017 à 2h29. 
        # [Consulté le 3 Avril 2019]. Disponible à l'adresse : https://stackoverflow.com/questions/45713887/add-values-from-two-dictionaries
        phage_amino_acid_dict = Counter(Counter(phage_amino_acid_dict) + Counter(protein_amino_acids))

        # get the weight
        seq = protein.sequence_AA
        # put off X to count the weight
        seq = seq.replace("X", "")
        analysed_seq = ProteinAnalysis(seq)
        protein_weight += analysed_seq.molecular_weight()

    # calcul mean
    tot_protein_of_the_phage = len(ProteinJson.getByOrganismID(phage.id))
    
    # dict contained mean std and weight for the phage
    dict_of_features = {}
    if get_features == True:
        # calcul mean of each amino acid
        for amino_acid_key in phage_amino_acid_dict.keys():
            dict_of_features['MEAN_AA_' + str(amino_acid_key)] = phage_amino_acid_dict[amino_acid_key] / tot_protein_of_the_phage

        # put weight in dict
        dict_of_features['WEIGHT (Da)'] = protein_weight

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
        