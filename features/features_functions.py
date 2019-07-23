# -*- coding: utf-8 -*-

import os
import sys
sys.path.insert(0, './')
import pandas as pd
import csv

from collections import Counter
from collections import OrderedDict
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import general_functions
from features import constants
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from objects_API.ProteinJ import ProteinJson
from objects_API.CoupleJ import CoupleJson
from objects_API.BacteriophageJ import BacteriophageJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def createMatrixFeatures(csv_file:str, matrix_of_features:list, organisms_designation:list):
    """
    | Read a csv file who contained features and fill matrix_of_features list with them.
    | Get organisms label from a csv and add them into organisms_designation

    :param csv_file: csv file with features
    :param matrix_of_features: matrix of features
    :param organisms_designation: list of organisms' designation

    :type csv_file: str
    :type matrix_of_features: list
    :type organisms_designation: list
    """
    # /!\ If a row description is not in the csv_file, add it /!\ 
    with open(csv_file, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            # Get the indexes of columns of features on the first iteration
            if 'MEAN_AA_M' in row:
                mean_aa_m = row.index('MEAN_AA_M')
            if 'MEAN_AA_E' in row:
                mean_aa_e = row.index('MEAN_AA_E')
            if 'MEAN_AA_N' in row:
                mean_aa_n = row.index('MEAN_AA_N')
            if 'MEAN_AA_Y' in row:
                mean_aa_y = row.index('MEAN_AA_Y')
            if 'MEAN_AA_K' in row:
                mean_aa_k = row.index('MEAN_AA_K')
            if 'MEAN_AA_F' in row:
                mean_aa_f = row.index('MEAN_AA_F')
            if 'MEAN_AA_I' in row:
                mean_aa_i = row.index('MEAN_AA_I')
            if 'MEAN_AA_A' in row:
                mean_aa_a = row.index('MEAN_AA_A')
            if 'MEAN_AA_H' in row:
                mean_aa_h = row.index('MEAN_AA_H')
            if 'MEAN_AA_L' in row:
                mean_aa_l = row.index('MEAN_AA_L')
            if 'MEAN_AA_V' in row:
                mean_aa_v = row.index('MEAN_AA_V')
            if 'MEAN_AA_Q' in row:
                mean_aa_q = row.index('MEAN_AA_Q')
            if 'MEAN_AA_R' in row:
                mean_aa_r = row.index('MEAN_AA_R')
            if 'MEAN_AA_T' in row:
                mean_aa_t = row.index('MEAN_AA_T')
            if 'MEAN_AA_D' in row:
                mean_aa_d = row.index('MEAN_AA_D')
            if 'MEAN_AA_S' in row:
                mean_aa_s = row.index('MEAN_AA_S')
            if 'MEAN_AA_W' in row:
                mean_aa_w = row.index('MEAN_AA_W')
            if 'MEAN_AA_C' in row:
                mean_aa_c = row.index('MEAN_AA_C')
            if 'MEAN_AA_G' in row:
                mean_aa_g = row.index('MEAN_AA_G')
            if 'MEAN_AA_P' in row:
                mean_aa_p = row.index('MEAN_AA_P')
            if 'MEAN_AA_X' in row:
                mean_aa_x = row.index('MEAN_AA_X')
            if 'MEAN_CE_H' in row:
                mean_ce_h = row.index('MEAN_CE_H')
            if 'MEAN_CE_N' in row:
                mean_ce_n = row.index('MEAN_CE_N')
            if 'MEAN_CE_O' in row:
                mean_ce_o = row.index('MEAN_CE_O')
            if 'MEAN_CE_C' in row:
                mean_ce_c = row.index('MEAN_CE_C')
            if 'MEAN_CE_S' in row:
                mean_ce_s = row.index('MEAN_CE_S')
            if 'WEIGHT (Da)' in row:
                weight = row.index('WEIGHT (Da)')
            if 'ISO POINT' in row:
                iso_point = row.index('ISO POINT')
            if 'AROMATICITY' in row:
                aromaticity = row.index('AROMATICITY')
                continue # On the first iteration we are on the first row and there is no data only labels

            # Handle empty case for X
            if row[mean_aa_x] == "":
                row[mean_aa_x] = 0.0

            #the number of features is the dimension of our final matrice (29 features)
            features = [float(row[mean_aa_a]), float(row[mean_aa_c]), float(row[mean_aa_d]), float(row[mean_aa_e]),
                        float(row[mean_aa_f]), float(row[mean_aa_g]), float(row[mean_aa_y]), float(row[mean_aa_h]),
                        float(row[mean_aa_i]), float(row[mean_aa_k]), float(row[mean_aa_l]), float(row[mean_aa_m]),
                        float(row[mean_aa_n]), float(row[mean_aa_p]), float(row[mean_aa_q]), float(row[mean_aa_r]),
                        float(row[mean_aa_s]), float(row[mean_aa_t]), float(row[mean_aa_v]), float(row[mean_aa_w]),
                        float(row[mean_ce_h]), float(row[mean_ce_n]), float(row[mean_ce_o]), float(row[mean_ce_c]),
                        float(row[mean_aa_x]), float(row[mean_ce_s]), float(row[weight]), float(row[iso_point]),
                        float(row[aromaticity])]

            matrix_of_features.append(features) # Add the features in the matrix
            organisms_designation.append(row[0]) # Add organisms designation
    f.close()

def createFeaturesFile(file_name:str, path:str, organism_id:int=None, is_phage:bool=False):
    """
    | Create a csv file with organisms features.
    | The organisms to compare are choosen according to the organism_id param.
    | Takes all couples with organism_id and select only the couples Clear Lysis.
    | Then extracts features from the remaining organisms.
    | If no organism_id referenced take the organisms in list_couple

    The features are :
        - mean of all amino acids
        - mean of all chemical elements
        - weight of the organism
        - iso eletrical point
        - aromaticity

    :param file_name: name of csv file
    :param path: path to save the file
    :param organism_id: id of organism
    :param is_phage: if true, compare bacteria. If false, compare bacteriophages

    :type file_name: str
    :type path: str
    :type organism_id: int
    :type is_phage: boolean
    """
    # Choose what type of lysis we want
    lysis_type = constants.ALL_CLEAR_LYSIS

    file_name = file_name
    path = path

    # Get the phages of all the couples with an interaction type
    # list_couple = network.getCouplesInteraction(CoupleJson.getAllAPI(), interaction_type=True)

    # Get the phages of specific couples 
    liste_couple = general_functions.getCouplesLysis(lysis_type)

    # Get couple from specific bacterie
    organism_dict = {}
    # Research bact or phage by ID
    if is_phage == False and organism_id != None:
            organism_dict['bacterium'] = organism_id
            liste_couple = (CoupleJson.getCouplesByFilterParameter(organism_dict))
    elif organism_id != None:
            organism_dict['bacteriophage'] = organism_id
            liste_couple = (CoupleJson.getCouplesByFilterParameter(organism_dict))

    if organism_id != None:
            # Select couple in function of the lysis
            liste_couple_final = []
            for couple in liste_couple:
                    if couple.lysis in lysis_type:
                            liste_couple_final.append(couple)
    else:
            liste_couple_final = liste_couple

    # phages to extract features will be here
    list_phages = []

    for couple in liste_couple_final:
            if not couple.bacteriophage in list_phages:
                    list_phages.append(couple.bacteriophage)

    phages_amino_acids = OrderedDict()

    for couple in liste_couple_final:
            if not couple.bacteriophage in phages_amino_acids.keys():
                    phage_amino_acid_dict = getFeaturesForAPhage(BacteriophageJson.getByID(couple.bacteriophage), active_percentage=False, get_features=True)
                    phages_amino_acids[BacteriophageJson.getByID(couple.bacteriophage).designation] = phage_amino_acid_dict


    # Source : 
    # ADAM SMITH [Pseudonyme], 2015. This is beyond simple if you can use pandas. 
    # Stackoverflow [en ligne]. 16 Juillet 2015 à 16h38. 
    # [Consulté le 3 Avril 2019]. Disponible à l'adresse : https://stackoverflow.com/questions/31436783/writing-dictionary-of-dictionaries-to-csv-file-in-a-particular-format
    df1 = pd.DataFrame.from_dict(data=phages_amino_acids, orient="index")
    df1.to_csv(os.path.join(path, file_name))
    df1 = df1.rename_axis('Phages designation', axis='columns')
    df1.to_csv(os.path.join(path, file_name), index_label= 'Phages designation')
    ending_message = "file " + file_name + " saved in " + path
    print(ending_message)

def getAllChemicalStructureOfAminoAcids(amino_acid:str):
    """
    | Get chemical strucutre of a given amino acid.

    :param amino_acid: amino acid

    :type amino_acid: str

    :return: chemical strucutre of the given amino acid
    :rtype: dictionnary {chemical element : occurence}
    """
    # Source : 
    # K.REDDY, Michael, 1998. "Amino acid CHEMICAL COMPOUND". 
    # ENCYCLOPAEDIA BRITANNICA [en ligne]. 20 Juillet 1998. 
    # [Consulté le 27 février 2019]. Disponible à l'adresse : https://www.britannica.com/science/amino-acid
    
    # Dict contains 20 proteins with their molecules
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
    | Get the frequence of each amino acid of a given protein.

    :param protein: protein object

    :type protein: ProteinJson

    :return: all the frequences of each amino acids and chemicals elements of the given protein
    :rtype: dictionnary {name amino acid : frequence}
    """
    # Dictionnary {amino acid : occurence} + {chimical element : occurence}
    protein_amino_acid_dict = OrderedDict()
    protein_chemical_element_dict = OrderedDict()

    # Get the proteique sequence of current prot
    sequence_of_amino_acid = protein.sequence_AA

    # Count all the amino acid in the sequence
    for amino_acid in sequence_of_amino_acid:
        if not amino_acid in protein_amino_acid_dict.keys():
            protein_amino_acid_dict[amino_acid] = 1
        else:
            protein_amino_acid_dict[amino_acid] += 1
    
        # Get chemical element of amino acid
        protein_chemical_element_dict = Counter(Counter(getAllChemicalStructureOfAminoAcids(amino_acid)) + Counter(protein_chemical_element_dict))

    # Transform the occurence of each amino acids in frequence
    for amino_acid_key, amino_acid_value in protein_amino_acid_dict.items():
        #divide the occurence by the length of the sequence
        protein_amino_acid_dict[amino_acid_key] = amino_acid_value / len(protein.sequence_AA)
    
    # Transform the occurence of each chemical element in frequence
    total_chemical_elements = protein_chemical_element_dict['H'] + protein_chemical_element_dict['O'] + protein_chemical_element_dict['N'] + protein_chemical_element_dict['C'] + protein_chemical_element_dict['S']

    for chemical_element_key, chemical_element_value in protein_chemical_element_dict.items():
        # Divide the occurence by the length of the sequence
        protein_chemical_element_dict[chemical_element_key] = chemical_element_value / total_chemical_elements 

    return protein_amino_acid_dict, protein_chemical_element_dict

def getFeaturesForAPhage(phage:BacteriophageJson, active_percentage:bool=False, get_features:bool=False):
    """
    | Get features for a given phage.

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

    # Get all amino acids of the phage
    for protein in list_prot:
        protein_amino_acids, protein_chemical_element = getFrequenceAllAminoAcidForAProtein(protein)
         
        # Add values from two dictionnaries in one
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

    # Calcul mean
    tot_protein_of_the_phage = len(ProteinJson.getByOrganismID(phage.id))
    
    # Dict contained mean std and weight for the phage
    dict_of_features = {}
    if get_features == True:
        # Calcul mean of each amino acid and chemical element
        for amino_acid_key in phage_amino_acid_dict.keys():
            dict_of_features['MEAN_AA_' + str(amino_acid_key)] = phage_amino_acid_dict[amino_acid_key] / tot_protein_of_the_phage
        
        for chemical_element_key in phage_chemical_element_dict.keys():
            dict_of_features['MEAN_CE_' + str(chemical_element_key)] = phage_chemical_element_dict[chemical_element_key] / tot_protein_of_the_phage

        # Put weight in dict
        dict_of_features['WEIGHT (Da)'] = protein_weight
        # Put aromaticity in dict
        dict_of_features['AROMATICITY'] = protein_aromatiticy
        # Put isoelectric_point in dict
        dict_of_features['ISO POINT'] = protein_isoeletric_point

        if 'MEAN_AA_X' not in dict_of_features.keys():
            dict_of_features['MEAN_AA_X'] = 0
        if 'MEAN_CE_Unknow' not in dict_of_features.keys():
            dict_of_features['MEAN_CE_Unknow'] = 0

    return dict_of_features

def converstionValuesInPercentage(dictionnary:dict):
    """
    | Convert all value of a dictionnary in percent.

    :param dictionnary: dictionnary

    :type dictionnary: dict

    """
    total_value = sum(dictionnary.values())
    for key, value in dictionnary.items():
        dictionnary[key] = (value * 100) / total_value
