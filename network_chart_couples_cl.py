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
from objects_API.StrainJ import StrainJson
from objects_API.SpecieJ import SpecieJson

from Bio import Entrez #to import data from GenBank
Entrez.email = "christophe.joyet@heig-vd.ch"

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def draw_graph(phages:list, bacterium:list, list_couples_lysis_type:list,
               is_png=True,
               node_size=1600, node_alpha=0.5,
               node_text_size=8,
               edge_alpha=0.5, edge_tickness=0.5,
               edge_text_pos=1.0,
               text_font='sans-serif',
               graph_name='network_graphic'):
    """
    draw a network graphics to bind phages with bacterium according to their lysis attribute

    :param phages: list of phages
    :param bacterium: list of bacterium
    :param list_couple_lysis_type: list of couples
    :param is_png: if true, save the graph as a png image
    :param node_size: size of node
    :param node_alpha: clearness of node_alpha
    :param node_text_size: size of text in node
    :param edge_alpha: clearness of node_alpha
    :param edge_tickness: tichness of edge,
    :param edge_text_pos: position of text
    :param text_font: text font
    :param graph_name: name of graphic

    :type phages: list
    :type bacterium: list
    :type list_couple_lysis_type: list
    :type is_png: boolean
    :type node_size: int
    :type node_size: int
    :type node_alpha: int
    :type node_text_size: int
    :type edge_alpha: int
    :type edge_tickness: int
    :type edge_text_pos: int
    :type text_font: string
    :type graph_name: string

    """
    fig, ax = plt.subplots(figsize=(200, 100))

    ax.set_title('Network between phages and bacteries', fontsize=16)
    
    #get all different phages
    nodes_phages = []
    #get all different bacterium
    nodes_bacterium = []
    
    #get the name of each bacterium (strain + species)
    for couple in list_couples_lysis_type:
        strain_id = BacteriumJson.getByID(couple.bacterium).strain
        strain_designation = StrainJson.getByID(strain_id).designation
        specie_designation = SpecieJson.getByID(StrainJson.getByID(strain_id).specie).designation
        bacterium_designation = specie_designation + '-' +  strain_designation
        if not bacterium_designation in nodes_bacterium:
            nodes_bacterium.append(bacterium_designation)
        
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
    
    #save graph in png or display it
    if is_png:
        plt.savefig('../../statistiques/NETWORK/' + graph_name + '.png')
    else:
        plt.show()

def getAllOfCouples():
    """
    get all couples of DB Inphinity

    :return: list of couples
    :rtype: list
    """
    list_couples = CoupleJson.getAllAPI()
    return list_couples


def getCouplesLysis(lysis_type):
    """
    get all couples of DB Inphinity according with the lysis type

    :return: list of couples
    :rtype: list
    """

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
    """
    get couples dictionnary according to their taxonomie

    :return: list of couples
    :rtype: list
    """
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

def getCouplesInteraction(couples:list, interaction_type:bool):
    """
    get couples according to their interaction 
    :param couple: list of couples
    :param interaction: type of the interaction (True or False)

    :type couple: CoupleJson
    :type interaction_type: bool

    :return: list of couples
    :rtype: list
    """
    list_couple = []
    for couple in couples:
        if couple.interaction_type == interaction_type:
            list_couple.append(couple)
    
    return list_couple

def getAccessionNumber(organism_list:list):
    """
    get accession number of an organism list

    :return: list of organism
    :rtype: list
    """
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
