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
from objects_API.OrganismJ import OrganismJson
from objects_API.BacteriophageJ import BacteriophageJson


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

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

#choose what type of lysis we want
lysis_type = CLEAR_LYSIS

list_couples_lysis_type = []
list_couples_lysis_type = network.getCouplesLysis(lysis_type)

print("Nombre de couples : " + str(len(list_couples_lysis_type)))

#defining two correlation tables between phages and bacteriums
phages = []
bacterium = []

for couple in list_couples_lysis_type:
    phages.append(BacteriophageJson.getByID(couple.bacteriophage).designation)
    bacterium.append(BacteriumJson.getByID(couple.bacterium).strain)

# network graph
network.draw_graph(phages, bacterium, list_couples_lysis_type, graph_name='network_cl')