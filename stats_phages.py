import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import network_chart_couples_cl as network
from collections import Counter

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.ProteinJ import ProteinJson
from objects_API.WholeDNAJ import WholeDNAJson
from objects_API.ContigJ import ContigJson
from objects_API.FamilyJ import FamilyJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.PersonResponsibleJ import PersonResponsibleJson
from objects_API.StrainJ import StrainJson
from objects_API.SpecieJ import SpecieJson
from objects_API.SourceDataJ import SourceDataJson
from objects_API.OrganismJ import OrganismJson


from objects_API.SourceDataJ import SourceDataJson

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
                                'O' : 2}}  
    
    if amino_acid in molecule_amino_acid.keys():
        return molecule_amino_acid[amino_acid]
        
def getAllAminoAcidForAProtein(protein:ProteinJson):
    """
    get all amino acid of a given protein

    :param protein: protein object

    :type protein: ProteinJson

    :return: all the amino acids of the given protein
    :rtype: dictionnary {name amino acid : occurence}
    """
    protein_AA_dict = {}
    protein_seq_AA = protein.sequence_AA
    for AA in protein_seq_AA:
        if not AA in protein_AA_dict.keys():
            protein_AA_dict[AA] = 1
        else:
            protein_AA_dict[AA] += 1

    return protein_AA_dict

def getAllAminoAcidForAPhage(phage:BacteriophageJson, active_percentage:bool, get_mean_and_std:bool):
    """
    get all amino acid of a given phage

    :param phage: phage object
    :param active_percentage: provide result in percentage
    :param get_mean_and_std: add mean and std in the dictionnary
    


    :type phage: BacteriophageJson
    :type active_percentage: bool
    :type get_mean_and_std: bool


    :return: all the amino acids of the given phage
    :rtype: dictionnary {name amino acid : occurence}
    """
    list_prot = ProteinJson.getByOrganismID(phage.id)
    phage_AA_dict = {}

    #get all amino acids of the phage
    #source : https://stackoverflow.com/questions/45713887/add-values-from-two-dictionaries
    for protein in list_prot:
        protein_AA = getAllAminoAcidForAProtein(protein)
        #add values from two dictionnaries in one
        phage_AA_dict = Counter(Counter(phage_AA_dict) + Counter(protein_AA))
    
    #calcul mean and std
    if get_mean_and_std == True:
        phage_AA_dict['mean'] = calculMeanFromDict(phage_AA_dict)
        phage_AA_dict['std'] = calculStdFromDict(phage_AA_dict)
    
    #conversion of the results in percentage
    if active_percentage == True:
        total_AA = sum(phage_AA_dict.values())
        for key, value in phage_AA_dict.items():
            if key != 'std' and key != 'mean':
                phage_AA_dict[key] = (value * 100) / total_AA
   
    return phage_AA_dict
'''
    def getAllChemicalStructureForAPhage(phage:BacteriophageJson, active_percentage:bool, get_mean_and_std:bool):
        """
        get all chemical structure for a phage

        :param phage: phage object
        :param active_percentage: provide result in percentage
        :param get_mean_and_std: add mean and std in the dictionnary
        
        :type phage: BacteriophageJson
        :type active_percentage: bool
        :type get_mean_and_std: bool

        :return: all the chemical strucure of the given phage
        :rtype: dictionnary {name amino acid : occurence}
        """
        list_prot = ProteinJson.getByOrganismID(phage.id)
        phage_AA_dict = getAllAminoAcidForAPhage(phage, False, False)
        phage_chemical_struct_dict = {}

        #get all amino acids of the phage
        #source : https://stackoverflow.com/questions/45713887/add-values-from-two-dictionaries
        for amino_acid, value in phage_AA_dict.items():
            if not amino_acid in phage_chemical_struct_dict:
                phage_chemical_struct_dict[amino_acid] = 
            phage_chemical_struct_dict = getAllAminoAcidForAProtein(protein)
            #add values from two dictionnaries in one
            phage_chemical_struct_dict = Counter(Counter(phage_chemical_struct_dict) + Counter(phage_chemical_struct_dict))
        
        #calcul mean and std
        if get_mean_and_std == True:
            phage_chemical_struct_dict['mean'] = calculMeanFromDict(phage_chemical_struct_dict)
            phage_chemical_struct_dict['std'] = calculStdFromDict(phage_chemical_struct_dict)
        
        #conversion of the results in percentage
        if active_percentage == True:
            total_atomium = sum(phage_chemical_struct_dict.values())
            for key, value in phage_chemical_struct_dict.items():
                if key != 'std' and key != 'mean':
                    phage_chemical_struct_dict[key] = (value * 100) / total_atomium
    
        return phage_chemical_struct_dict
'''
def calculMeanFromDict(dictionnary:dict):
    """
    get mean value from a dictionnary

    :param dictionnary: dictionnary

    :type dictionnary: dict

    :return: mean of all value contained in dictionnary
    :rtype: float
    """
    #get values from phage_AA_dict to calculate std and mean
    values = []
    for key,value in dictionnary.items():
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
    #get values from phage_AA_dict to calculate std and mean
    values = []
    for key,value in dictionnary.items():
        values.append(value)
    return np.std(values)

#==============================================================================
#==============================================================================

#list_couple = network.getCouplesInteraction(CoupleJson.getAllAPI(), True)
list_couple = network.getCouplesLysis([5,8,9])
list_phages = []
for couple in list_couple:
    if not couple.bacteriophage in list_phages:
        list_phages.append(couple.bacteriophage)

print("nombre de phages différents : " + str(len(list_phages)))

list_AA = {}
for couple in list_couple:
    if not couple.bacteriophage in list_AA.keys():
        phage_AA_dict = getAllAminoAcidForAPhage(BacteriophageJson.getByID(couple.bacteriophage), active_percentage=False, get_mean_and_std=True)
        list_AA[couple.bacteriophage] = phage_AA_dict


#source : https://stackoverflow.com/questions/31436783/writing-dictionary-of-dictionaries-to-csv-file-in-a-particular-format
df = pd.DataFrame.from_dict(data=list_AA, orient="index")
df.to_csv("stats_phages_ALL_CL_number.csv")
