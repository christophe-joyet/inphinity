#!/usr/bin/env python3
# File : stats_phages_script.py
# Author : Christophe Joyet
# Date : Juin 2019

import csv
import os
import pandas as pd
import network_chart_couples_cl as network
import get_features_from_organism as sp


from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from collections import OrderedDict
from objects_API.CoupleJ import CoupleJson
from objects_API.ProteinJ import ProteinJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.OrganismJ import OrganismJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def getFeaturesFromOrganismScript(file_name:str, path:str, organism_id:int=None, is_phage:bool=False):
        #=============================================================================================
        #Constantes declarations
        #=============================================================================================

        #Lysis Type
        CLEAR_LYSIS = 5
        SEMI_CLEAR_LYSIS = 6
        OPAQUE_LYSIS = 7
        #dilution > 1e7
        CLEAR_LYSIS_1E7PLUS = 8
        SEMI_CLEAR_LYSIS_1E7PLUS = 10
        #dilution < 1e7
        CLEAR_LYSIS_1E7MINUS = 9
        SEMI_CLEAR_LYSIS_1E7MINUS = 11

        ALL_CLEAR_LYSIS = [CLEAR_LYSIS, CLEAR_LYSIS_1E7PLUS, CLEAR_LYSIS_1E7MINUS]
        ALL_SEMI_CLEAR_LYSIS = [SEMI_CLEAR_LYSIS, SEMI_CLEAR_LYSIS_1E7PLUS, SEMI_CLEAR_LYSIS_1E7MINUS]

        # =============================================================================================
        # =============================================================================================

        # choose what type of lysis we want
        lysis_type = ALL_CLEAR_LYSIS

        # =============================================================================================
        # =============================================================================================

        file_name = file_name
        path = path

        # =============================================================================================
        # =============================================================================================

        # get the phages of all the couples with an interaction type
        # list_couple = network.getCouplesInteraction(CoupleJson.getAllAPI(), interaction_type=True)

        # get the phages of specific couples 
        liste_couple = network.getCouplesLysis(lysis_type)

        # get couple from specific bacterie
        organism_dict = {}
        # research bact or phage by ID
        if is_phage == False and organism_id != None:
                organism_dict['bacterium'] = organism_id
                liste_couple = (CoupleJson.getCouplesByFilterParameter(organism_dict))
        elif organism_id != None:
                organism_dict['bacteriophage'] = organism_id
                liste_couple = (CoupleJson.getCouplesByFilterParameter(organism_dict))


        if organism_id != None:
                # select couple in function of the lysis
                liste_couple_final = []
                for couple in liste_couple:
                        if couple.lysis in lysis_type:
                                liste_couple_final.append(couple)
        else:
                liste_couple_final = liste_couple

        # =====================================================================================================
        # =====================================================================================================

        # phages to extract features will be here
        list_phages = []

        for couple in liste_couple_final:
                if not couple.bacteriophage in list_phages:
                        list_phages.append(couple.bacteriophage)

        print("nombre de phages différents : " + str(len(list_phages)))

        phages_amino_acids = OrderedDict()

        for couple in liste_couple_final:
                if not couple.bacteriophage in phages_amino_acids.keys():
                        phage_amino_acid_dict = sp.getFeaturesForAPhage(BacteriophageJson.getByID(couple.bacteriophage), active_percentage=False, get_features=True)
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

        # =====================================================================================================
        # =====================================================================================================

#getFeaturesFromOrganismScript(file_name="features_70_phages_clear_lysis.csv", path="../../extraction_features/AA_CE_WEIGHT")