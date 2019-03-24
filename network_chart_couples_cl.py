# libraries
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.OrganismJ import OrganismJson
from objects_API.BacteriophageJ import BacteriophageJson

from Bio import Entrez #to import data from GenBank
Entrez.email = "christophe.joyet@heig-vd.ch"

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

#Constantes declarations

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

#==================================IMPLEMENTATION OF FUNCTIONS===============================

#Draw network graph to show relations between phages and bacterium
def draw_graph(phages:list, bacterium:list, list_couples_lysis_type:list,
               node_size=1600, node_alpha=0.5,
               node_text_size=12,
               edge_alpha=0.6, edge_tickness=0.5,
               edge_text_pos=1.0,
               text_font='sans-serif'):

    fig, ax = plt.subplots(figsize=(200, 100))

    ax.set_title('Network between phages and bacteries', fontsize=16)
    
    #get all different phages
    nodes_phages = []
    #get all different bacterium
    nodes_bacterium = []
    
    for couple in list_couples_lysis_type:
        couple_bacterium = BacteriumJson.getByID(couple.bacterium).acc_number
        if not couple_bacterium in nodes_bacterium:
            nodes_bacterium.append(couple_bacterium)
        
        #get phages' designation
        couple_bacteriophage = BacteriophageJson.getByID(couple.bacteriophage).designation
        if not couple_bacteriophage in nodes_phages:
            nodes_phages.append(couple_bacteriophage)

    print("Nombre de phages différents : " + str(len(nodes_phages)))
    print("Nombre de bactéries différentes : " + str(len(nodes_bacterium)))

    # all the nodes in our graph
    nodes = set(nodes_phages + nodes_bacterium)

    # create networkx graph
    G=nx.Graph()

    # add nodes
    for node in nodes:
        G.add_node(node)

    # add edges
    i = 0
    while(i < len(phages)):
        G.add_edge(phages[i], bacterium[i])
        i += 1

    # draw graph
    # layout
    graph_pos=nx.spring_layout(G)

    # draw graph
    #defining nodes features
    nx.draw_networkx_nodes(G,graph_pos,nodelist=nodes_phages,node_size=node_size, 
                           alpha=node_alpha, node_color='r')
    
    nx.draw_networkx_nodes(G,graph_pos,nodelist=nodes_bacterium,node_size=node_size, 
                           alpha=node_alpha, node_color='g')
        
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
                           alpha=edge_alpha,edge_color='b')
   
    #display ID of bacterium and phages
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,
                            font_family=text_font)

    #show graph
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel('Green = Bacterium\n Red = Phages')
    #plt.show()
    #save graph in png
    plt.savefig('network_CL_E7P_E7M.png')

#get all couples of the DB inphinity
def getAllOfCouples():
    list_couples = CoupleJson.getAllAPI()
    return list_couples

#get couples according to their lysis classification
def getCouplesLysis(lysis_type):
    
    #verify if lysis is an int 
    if isinstance(lysis_type, int):
        tmp = lysis_type
        #change lysis_type in list
        lysis_type = []
        lysis_type.append(tmp)

    list_couples = getAllOfCouples()
    list_couples_lysis_type = []

    #Get all couple with lysis type
    for couple in list_couples:
        if couple.lysis in lysis_type:
            list_couples_lysis_type.append(couple)

    return list_couples_lysis_type

#get couples dictionnary according to their taxonomie
def getCouplesTaxonomie(couples_list:list):
    #Level interaction
    STRAIN = 1
    SPECIES = 2
    GENUS = 3
    FAMILY = 4

    taxonomie_dictionnary = {STRAIN:0, SPECIES:0, GENUS:0, FAMILY:0}
    for couples in couples_list:
        if couples.level in taxonomie_dictionnary:
            taxonomie_dictionnary[couples.level] += 1

    return taxonomie_dictionnary

#get accession number
def getAccessionNumber(organism_list:list):
    accessions_number_list = []
    for organism in organism_list:
        accessions_number_list.append(organism.acc_number)
    
    return accessions_number_list

#get definition of an organism using Genbank
'''def getOrganismDefinitionFromGenBank(organism_list:list):
    organism_definition_list = []
    i = 0

    while(i < len(organism_list)):
        query = str(organism_list[i])
        handle = Entrez.esearch(db='nucleotide', term=query)
        record = Entrez.read(handle)
        gi = record['IdList']
        protein_info = Entrez.efetch(db="nucleotide",id=gi,rettype="gb",retmode="xml")
        protein_info = Entrez.read(protein_info)

        organism_definition_list.append()

        i += 1'''
#==================================IMPLEMENTATION OF FUNCTIONS END===============================

#choose what type of lysis we want
lysis_type = CLEAR_LYSIS

list_couples_lysis_type = []
list_couples_lysis_type = getCouplesLysis(lysis_type)

print("Nombre de couples clear lysis : " + str(len(list_couples_lysis_type)))

#defining two correlation tables between phages and bacteriums
phages = []
bacterium = []

for couple in list_couples_lysis_type:
    phages.append(BacteriophageJson.getByID(couple.bacteriophage).designation)
    bacterium.append(BacteriumJson.getByID(couple.bacterium).acc_number)

'''for couple in list_couples_lysis_type:
    phages.append(couple.bacteriophage)
    bacterium.append(couple.bacterium)'''

print("phage" + str(len(phages)))
print("bacterium" + str(len(bacterium)))

couples_tax = getCouplesTaxonomie(list_couples_lysis_type)
print(couples_tax)
print(len(list_couples_lysis_type))
# network graph
draw_graph(phages, bacterium, list_couples_lysis_type)