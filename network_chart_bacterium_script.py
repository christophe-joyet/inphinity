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

#=============================================================================================
#=============================================================================================

#defining two correlation tables between phages and bacteriums
phages = []
bacterium = []

list_couple = []
bacterium_dict = {}
bacterium_dict['bacterium'] = 70

liste_couple = CoupleJson.getCouplesByFilterParameter(bacterium_dict)
print(liste_couple)

for couple in liste_couple:
    phages.append(BacteriophageJson.getByID(couple.bacteriophage).designation)
    #get the name of bacterium (strain designation + species designation)
    strain_id = BacteriumJson.getByID(couple.bacterium).strain
    strain_designation = StrainJson.getByID(strain_id).designation
    specie_designation = SpecieJson.getByID(StrainJson.getByID(strain_id).specie).designation
    bacterium.append(specie_designation + '-' +  strain_designation)

# network graph
network.draw_graph(phages, bacterium, liste_couple, graph_name='bacterium_try_1')
