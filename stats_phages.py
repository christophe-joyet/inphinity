import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import os
import network_chart_couples_cl as network
from collections import Counter
from collections import OrderedDict
from scipy.stats import kurtosis

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
        
def getAllAminoAcidForAProtein(protein:ProteinJson):
    """
    get all amino acid of a given protein

    :param protein: protein object

    :type protein: ProteinJson

    :return: all the amino acids of the given protein
    :rtype: dictionnary {name amino acid : occurence}
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

    return protein_amino_acid_dict

def getAllAminoAcidForAPhage(phage:BacteriophageJson, active_percentage:bool=False, get_features:bool=False):
    """
    get all amino acid of a given phage

    :param phage: phage object
    :param active_percentage: provide result in percentage
    :param get_features: add mean and std in the dictionnary
    
    :type phage: BacteriophageJson
    :type active_percentage: bool
    :type get_features: bool


    :return: all the amino acids of the given phage
    :rtype: dictionnary {name amino acid : occurence}
    """
    list_prot = ProteinJson.getByOrganismID(phage.id)
    phage_amino_acid_dict = OrderedDict()

    #source : https://stackoverflow.com/questions/45713887/add-values-from-two-dictionaries
    #get all amino acids of the phage
    for protein in list_prot:
        protein_amino_acids = getAllAminoAcidForAProtein(protein)
        #add values from two dictionnaries in one
        phage_amino_acid_dict = Counter(Counter(phage_amino_acid_dict) + Counter(protein_amino_acids))

    #calcul mean and std
    if get_features == True:
        mean = calculMeanFromDict(phage_amino_acid_dict)
        med = calculMedianFromDict(phage_amino_acid_dict)
        std = calculStdFromDict(phage_amino_acid_dict)
        var = calculVarianceFromDict(phage_amino_acid_dict)
        kurt = calculKurtosisFromDict(phage_amino_acid_dict)
    
    #conversion of the results in percentage
    if active_percentage == True:
        converstionValuesInPercentage(phage_amino_acid_dict)
    
    #add mean and std in dictionnary
    if get_features == True:
        phage_amino_acid_dict['P_AA_MEAN'] = mean
        phage_amino_acid_dict['P_AA_MED'] = med
        phage_amino_acid_dict['P_AA_STD'] = std
        phage_amino_acid_dict['P_AA_VAR'] = var
        phage_amino_acid_dict['P_AA_KURTOSIS'] = kurt

    return phage_amino_acid_dict

def getAllChemicalStructureForAPhage(phage:BacteriophageJson, active_percentage:bool=False, get_features:bool=False):
    """
    get all chemical structure for a phage

    :param phage: phage object
    :param active_percentage: provide result in percentage
    :param get_features: add mean and std in the dictionnary
    
    :type phage: BacteriophageJson
    :type active_percentage: bool
    :type get_features: bool

    :return: all the chemical strucure of the given phage
    :rtype: dictionnary {name amino acid : occurence}
    """
    phage_amino_acid_dict = getAllAminoAcidForAPhage(phage)
    phage_chemical_struct_dict = OrderedDict()

    #source : https://stackoverflow.com/questions/45713887/add-values-from-two-dictionaries
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
        
    #calcul mean and std
    if get_features == True:
        mean = calculMeanFromDict(phage_chemical_struct_dict)
        med = calculMedianFromDict(phage_chemical_struct_dict)
        std = calculStdFromDict(phage_chemical_struct_dict)
        var = calculVarianceFromDict(phage_chemical_struct_dict)
        kurt = calculKurtosisFromDict(phage_chemical_struct_dict)

    #conversion of the results in percentage
    if active_percentage == True:
        converstionValuesInPercentage(phage_chemical_struct_dict)  
    
    #add mean and std in dictionnary
    if get_features == True:
        phage_chemical_struct_dict['P_MOLECULES_MEAN'] = mean
        phage_chemical_struct_dict['P_MOLECULES_MED'] = med
        phage_chemical_struct_dict['P_MOLECULES_STD'] = std
        phage_chemical_struct_dict['P_MOLECULES_VAR'] = var
        phage_chemical_struct_dict['P_MOLECULES_KURTOSIS'] = kurt

    return phage_chemical_struct_dict

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

    
def calculMeanFromDict(dictionnary:dict):
    """
    get mean value from a dictionnary

    :param dictionnary: dictionnary

    :type dictionnary: dict

    :return: mean of all value contained in dictionnary
    :rtype: float
    """
    values = []
    for value in dictionnary.values():
        values.append(value)
    return np.mean(values)

def calculStdFromDict(dictionnary:dict):
    """
    get standart deviation value from a dictionnary

    :param dictionnary: dictionnary

    :type dictionnary: dict

    :return: std of all value contained in dictionnary
    :rtype: float
    """
    values = []
    for value in dictionnary.values():
        values.append(value)
    return np.std(values)

def calculVarianceFromDict(dictionnary:dict):
    """
    get variance value from a dictionnary

    :param dictionnary: dictionnary

    :type dictionnary: dict

    :return: variance of all value contained in dictionnary
    :rtype: float
    """
    values = []
    for value in dictionnary.values():
        values.append(value)
    return np.var(values)

def calculKurtosisFromDict(dictionnary:dict):
    """
    get kurtosis value from a dictionnary

    :param dictionnary: dictionnary

    :type dictionnary: dict

    :return: kurtosis of all value contained in dictionnary
    :rtype: float
    """
    values = []
    for value in dictionnary.values():
        values.append(value)
    return kurtosis(values)

def calculMedianFromDict(dictionnary:dict):
    """
    get median value from a dictionnary

    :param dictionnary: dictionnary

    :type dictionnary: dict

    get median value from a dictionnary
    :return: kurtosis of all value contained in dictionnary
    :rtype: float
    """
    values = []
    for value in dictionnary.values():
        values.append(value)
    return np.median(values)
        
#==============================================================================
#==============================================================================

#list_couple = network.getCouplesInteraction(CoupleJson.getAllAPI(), interaction_type=True)
list_couple = network.getCouplesLysis([5, 8, 9])
list_phages = []
for couple in list_couple:
    if not couple.bacteriophage in list_phages:
        list_phages.append(couple.bacteriophage)

print("nombre de phages différents : " + str(len(list_phages)))
#Supprimer pour inverser
#'''
phages_amino_acids = OrderedDict()
for couple in list_couple:
    if not couple.bacteriophage in phages_amino_acids.keys():
        phage_amino_acid_dict = getAllAminoAcidForAPhage(BacteriophageJson.getByID(couple.bacteriophage), active_percentage=True, get_features=True)
        phages_amino_acids[couple.bacteriophage] = phage_amino_acid_dict

print("amino_acids done")
'''
phages_molecules = OrderedDict()
for couple in list_couple:
    if not couple.bacteriophage in phages_molecules.keys():
        phage_molecule_dict = getAllChemicalStructureForAPhage(BacteriophageJson.getByID(couple.bacteriophage), active_percentage=True, get_features=True)
        phages_molecules[couple.bacteriophage] = phage_molecule_dict

#source : https://stackoverflow.com/questions/31436783/writing-dictionary-of-dictionaries-to-csv-file-in-a-particular-format
'''
df1 = pd.DataFrame.from_dict(data=phages_amino_acids, orient="index")
df1.to_csv(os.path.join('../../statistiques/CSV/', r"fichier_test.csv"))
'''
df2 = pd.DataFrame.from_dict(data=phages_molecules, orient="index")
df2.to_csv(os.path.join('../../statistiques/CSV/', r"fichier_test.csv"))
'''
