
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


def getAllChemicalStructureOfAminoAcids(amino_acid:str):
    """
    get chemical strucutre of a given amino acid
    :param amino_acid: amino acid
    :type amino_acid: str
    :return: chemical strucutre of the given amino acid
    :rtype: dictionnary {atomium : occurence}
    """
    # Source : 
    # K.REDDY, Michael, 1998. "Amino acid CHEMICAL COMPOUND". 
    # ENCYCLOPAEDIA BRITANNICA [en ligne]. 20 Juillet 1998. 
    # [Consulté le 27 février 2019]. Disponible à l'adresse : https://www.britannica.com/science/amino-acid
    
    #dict contains 20 proteins with their molecules
    molecule_amino_acid = { 'A':{
                                'H' : 7,
                                'N' : 1,
                                'C' : 3,
                                'O' : 2,
                                'S' : 0},
                            'R':{
                                'H' : 14,
                                'N' : 4,
                                'C' : 6,
                                'O' : 2,
                                'S' : 0},
                            'N':{
                                'H' : 8,
                                'N' : 2,
                                'C' : 4,
                                'O' : 3,
                                'S' : 0},
                            'D':{
                                'H' : 7,
                                'N' : 1,
                                'C' : 4,
                                'O' : 4,
                                'S' : 0},
                            'C':{
                                'H' : 7,
                                'N' : 1,
                                'C' : 3,
                                'O' : 2,
                                'S' : 1},
                            'E':{
                                'H' : 9,
                                'N' : 1,
                                'C' : 5,
                                'O' : 4,
                                'S' : 0},
                            'Q':{
                                'H' : 10,
                                'N' : 2,
                                'C' : 5,
                                'O' : 3,
                                'S' : 0},
                            'G':{
                                'H' : 5,
                                'N' : 1,
                                'C' : 2,
                                'O' : 2,
                                'S' : 0},
                            'H':{
                                'H' : 9,
                                'N' : 3,
                                'C' : 6,
                                'O' : 2,
                                'S' : 0},
                            'I':{
                                'H' : 13,
                                'N' : 1,
                                'C' : 6,
                                'O' : 2,
                                'S' : 0},
                            'L':{
                                'H' : 13,
                                'N' : 1,
                                'C' : 6,
                                'O' : 2,
                                'S' : 0},
                            'K':{
                                'H' : 14,
                                'N' : 2,
                                'C' : 6,
                                'O' : 2,
                                'S' : 0},
                            'M':{
                                'H' : 11,
                                'N' : 1,
                                'C' : 5,
                                'O' : 2,
                                'S' : 1},
                            'F':{
                                'H' : 11,
                                'N' : 1,
                                'C' : 9,
                                'O' : 2,
                                'S' : 0},
                            'P':{
                                'H' : 9,
                                'N' : 1,
                                'C' : 5,
                                'O' : 2,
                                'S' : 0},
                            'S':{
                                'H' : 7,
                                'N' : 1,
                                'C' : 3,
                                'O' : 3,
                                'S' : 0},
                            'T':{
                                'H' : 9,
                                'N' : 1,
                                'C' : 4,
                                'O' : 3,
                                'S' : 0},
                            'W':{
                                'H' : 12,
                                'N' : 2 ,
                                'C' : 11,
                                'O' : 2,
                                'S' : 0},
                            'Y':{
                                'H' : 11,
                                'N' : 1,
                                'C' : 9,
                                'O' : 3,
                                'S' : 0},
                            'V':{
                                'H' : 11,
                                'N' : 1,
                                'C' : 5,
                                'O' : 2,
                                'S' : 0}}  
    
    if amino_acid in molecule_amino_acid.keys():
        return molecule_amino_acid[amino_acid]
    else:
        return {'Unknown' : 1}

def getFrequenceAllAminoAcidForAProtein(protein:ProteinJson):
    """
    get the frequence of each amino acid of a given protein

    :param protein: protein object

    :type protein: ProteinJson

    :return: all the frequences of each amino acids and chemicals elements of the given protein
    :rtype: dictionnary {name amino acid : frequence}
    """
    #dictionnary {amino acid : occurence} + {chimical element : occurence}
    protein_amino_acid_dict = OrderedDict()
    protein_chemical_element_dict = OrderedDict()

    #get the proteique sequence of current prot
    sequence_of_amino_acid = protein.sequence_AA

    #count all the amino acid in the sequence
    for amino_acid in sequence_of_amino_acid:
        if not amino_acid in protein_amino_acid_dict.keys():
            protein_amino_acid_dict[amino_acid] = 1
        else:
            protein_amino_acid_dict[amino_acid] += 1

        #get chemical element of amino acid
        protein_chemical_element_dict = Counter(Counter(getAllChemicalStructureOfAminoAcids(amino_acid)) + Counter(protein_chemical_element_dict))

    #transform the occurence of each amino acids in frequence
    for amino_acid_key, amino_acid_value in protein_amino_acid_dict.items():
        #divide the occurence by the length of the sequence
        protein_amino_acid_dict[amino_acid_key] = amino_acid_value / len(protein.sequence_AA)
    
    #transform the occurence of each chemical element in frequence
    total_chemical_elements = protein_chemical_element_dict['H'] + protein_chemical_element_dict['O'] + protein_chemical_element_dict['N'] + protein_chemical_element_dict['C'] + protein_chemical_element_dict['S']

    for chemical_element_key, chemical_element_value in protein_chemical_element_dict.items():
        #divide the occurence by the length of the sequence
        protein_chemical_element_dict[chemical_element_key] = chemical_element_value / total_chemical_elements 

    return protein_amino_acid_dict, protein_chemical_element_dict

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
    phage_chemical_element_dict = OrderedDict()
    protein_weight = 0
    protein_aromatiticy = 0
    protein_isoeletric_point = 0

    #get all amino acids of the phage
    for protein in list_prot:
        protein_amino_acids, protein_chemical_element = getFrequenceAllAminoAcidForAProtein(protein)
         
        #add values from two dictionnaries in one
        # Source : 
        # MSEIFERT [Pseudonyme], 2017. You can use collections.Counter which implements addition + that way. 
        # Stackoverflow [en ligne]. 16 août 2017 à 2h29. 
        # [Consulté le 3 Avril 2019]. Disponible à l'adresse : https://stackoverflow.com/questions/45713887/add-values-from-two-dictionaries
        phage_amino_acid_dict = Counter(Counter(phage_amino_acid_dict) + Counter(protein_amino_acids))
        phage_chemical_element_dict = Counter(Counter(phage_chemical_element_dict) + Counter(protein_chemical_element))

        # get the weight
        seq = protein.sequence_AA
        # put off X to count the weight
        seq = seq.replace("X", "")
        analysed_seq = ProteinAnalysis(seq)
        protein_weight += analysed_seq.molecular_weight()
        # get the aromaticity
        protein_aromatiticy += analysed_seq.aromaticity()
        # get isoeletric_point
        protein_isoeletric_point += analysed_seq.isoelectric_point()


    # calcul mean
    tot_protein_of_the_phage = len(ProteinJson.getByOrganismID(phage.id))
    
    # dict contained mean std and weight for the phage
    dict_of_features = {}
    if get_features == True:
        # calcul mean of each amino acid and chemical element
        for amino_acid_key in phage_amino_acid_dict.keys():
            dict_of_features['MEAN_AA_' + str(amino_acid_key)] = phage_amino_acid_dict[amino_acid_key] / tot_protein_of_the_phage
        
        for chemical_element_key in phage_chemical_element_dict.keys():
            dict_of_features['MEAN_CE_' + str(chemical_element_key)] = phage_chemical_element_dict[chemical_element_key] / tot_protein_of_the_phage

        # put weight in dict
        dict_of_features['WEIGHT (Da)'] = protein_weight
        # put aromaticity in dict
        dict_of_features['AROMATICITY'] = protein_aromatiticy
        # put isoelectric_point in dict
        dict_of_features['ISO POINT'] = protein_isoeletric_point

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
