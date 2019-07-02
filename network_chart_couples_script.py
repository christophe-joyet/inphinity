# libraries
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import network_chart_couples_cl as network
import time


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

def networkChartCouplesScript(parameter_type:str):
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
    if parameter_type == "CLEAR_LYSIS":
        lysis_type = CLEAR_LYSIS
    elif parameter_type == "SEMI_CLEAR_LYSIS":
         lysis_type = SEMI_CLEAR_LYSIS       
    elif parameter_type == "OPAQUE_LYSIS":
         lysis_type = OPAQUE_LYSIS
    elif parameter_type == "CLEAR_LYSIS_1E7PLUS":
        lysis_type = CLEAR_LYSIS_1E7PLUS
    elif parameter_type == "SEMI_CLEAR_LYSIS_1E7PLUS":
        lysis_type = SEMI_CLEAR_LYSIS_1E7PLUS
    elif parameter_type == "CLEAR_LYSIS_1E7MINUS":
        lysis_type = CLEAR_LYSIS_1E7MINUS
    elif parameter_type == "SEMI_CLEAR_LYSIS_1E7MINUS":
        lysis_type = SEMI_CLEAR_LYSIS_1E7MINUS
    elif parameter_type == "ALL_CLEAR_LYSIS":
        lysis_type = ALL_CLEAR_LYSIS
    elif parameter_type == "ALL_SEMI_CLEAR_LYSIS":
        lysis_type = ALL_SEMI_CLEAR_LYSIS

    # couples to analyse will be in list_couples_lysis_type
    list_couples_lysis_type = []
    list_couples_lysis_type = network.getCouplesLysis(lysis_type)

    print("Nombre de couples : " + str(len(list_couples_lysis_type)))

    # defining two correlation tables between phages and bacteriums
    phages = []
    bacterium = []

    for couple in list_couples_lysis_type:
        phages.append(BacteriophageJson.getByID(couple.bacteriophage).designation)
        #get the name of bacterium (strain designation + species designation)
        strain_id = BacteriumJson.getByID(couple.bacterium).strain
        strain_designation = StrainJson.getByID(strain_id).designation
        specie_designation = SpecieJson.getByID(StrainJson.getByID(strain_id).specie).designation
        bacterium.append(specie_designation + '\n' +  strain_designation + '\n' + str(couple.bacterium))

    # network graph
    network.draw_graph(phages, bacterium, list_couples_lysis_type, graph_name='all_clear_lysis', is_png=False)
