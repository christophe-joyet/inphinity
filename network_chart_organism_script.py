# libraries
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import network_chart_couples_cl as network


from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.StrainJ import StrainJson
from objects_API.SpecieJ import SpecieJson
from objects_API.OrganismJ import OrganismJson
from objects_API.BacteriophageJ import BacteriophageJson


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def networkChartOrganismScript(organism_id:int, is_phage=False):
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

        #=============================================================================================
        #=============================================================================================

        # choose what type of lysis we want
        lysis_type = ALL_CLEAR_LYSIS

        #=============================================================================================
        #=============================================================================================
        
        # get couple from specific bacterie
        organism_dict = {}
        # research bact or phage by ID
        if is_phage == False:
                organism_dict['bacterium'] = organism_id
        else:
                organism_dict['bacteriophage'] = organism_id

        liste_couple = (CoupleJson.getCouplesByFilterParameter(organism_dict))

        # select couple in function of the lysis
        liste_couple_final = []
        for couple in liste_couple:
                if couple.lysis in lysis_type:
                        liste_couple_final.append(couple)
        
        #=============================================================================================
        #=============================================================================================

        # defining two correlation tables between phages and bacteriums
        phages = []
        bacterium = []

        for couple in liste_couple_final:
                # get designation and phage id
                phages.append(BacteriophageJson.getByID(couple.bacteriophage).designation)
                # get the name of bacterium (strain designation + species designation) and his id
                strain_id = BacteriumJson.getByID(couple.bacterium).strain
                strain_designation = StrainJson.getByID(strain_id).designation
                specie_designation = SpecieJson.getByID(StrainJson.getByID(strain_id).specie).designation
                bacterium.append(specie_designation + '\n' +  strain_designation + '\n' + str(couple.bacterium))

        # network graph
        network.draw_graph(phages, bacterium, liste_couple_final, graph_name='essai', is_png=False)
