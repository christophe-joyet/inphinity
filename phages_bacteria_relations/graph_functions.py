# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, './')

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter

import general_functions as functions
from phages_bacteria_relations import constants
from objects_API.CoupleJ import CoupleJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.OrganismJ import OrganismJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.StrainJ import StrainJson
from objects_API.SpecieJ import SpecieJson
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def relationGraphOrganism(organism_id:int, is_phage=False):
    """
    | Draw a relation graph between phages and bacteria according
    | to the organism_id and clear lysis attribute. 

    :Remark: Change the list_couple_final for other comparisons

    :param organism_id: lysis
    :param is_phage: if the organism is a phage

    :type organism_id: int
    :type is_phage: boolean
    """
    # Choose what type of lysis we want
    lysis_type = constants.ALL_CLEAR_LYSIS
 
    # Get couple from specific bacterie
    organism_dict = {}
    # Research bact or phage by ID
    if is_phage == False:
            organism_dict['bacterium'] = organism_id
    else:
            organism_dict['bacteriophage'] = organism_id

    liste_couple = (CoupleJson.getCouplesByFilterParameter(organism_dict))

    # Select couple in function of the lysis
    liste_couple_final = []
    for couple in liste_couple:
            if couple.lysis in lysis_type:
                    liste_couple_final.append(couple)
    
    # Defining two correlation tables between phages and bacteriums
    phages = []
    bacterium = []

    for couple in liste_couple_final:
            phages.append(BacteriophageJson.getByID(couple.bacteriophage).designation)  # Get designation and phage id
            # Get the name of bacterium (strain designation + species designation) and his id
            strain_id = BacteriumJson.getByID(couple.bacterium).strain
            strain_designation = StrainJson.getByID(strain_id).designation
            specie_designation = SpecieJson.getByID(StrainJson.getByID(strain_id).specie).designation
            bacterium.append(specie_designation + '-' +  strain_designation + '\n' + str(couple.bacterium))

    # Draw network graph
    draw_graph(phages, bacterium, liste_couple_final, graph_name='graph', is_png=False)

def relationGraphCouples(lysis_attribut:str):
    """
    | Draw a relation graph between phages and bacteria according
    | to the lysis attribute.

    :param lysis_attribut: lysis

    :type lysis_attribut: str
    """

    # Choose what type of lysis we want
    if lysis_attribut == 'CLEAR_LYSIS':
        lysis_type = constants.CLEAR_LYSIS
    elif lysis_attribut == 'SEMI_CLEAR_LYSIS':
         lysis_type = constants.SEMI_CLEAR_LYSIS       
    elif lysis_attribut == 'OPAQUE_LYSIS':
         lysis_type = constants.OPAQUE_LYSIS
    elif lysis_attribut == 'CLEAR_LYSIS_1E7PLUS':
        lysis_type = constants.CLEAR_LYSIS_1E7PLUS
    elif lysis_attribut == 'SEMI_CLEAR_LYSIS_1E7PLUS':
        lysis_type = constants.SEMI_CLEAR_LYSIS_1E7PLUS
    elif lysis_attribut == 'CLEAR_LYSIS_1E7MINUS':
        lysis_type = constants.CLEAR_LYSIS_1E7MINUS
    elif lysis_attribut == 'SEMI_CLEAR_LYSIS_1E7MINUS':
        lysis_type = constants.SEMI_CLEAR_LYSIS_1E7MINUS
    elif lysis_attribut == 'ALL_CLEAR_LYSIS':
        lysis_type = constants.ALL_CLEAR_LYSIS
    elif lysis_attribut == 'ALL_SEMI_CLEAR_LYSIS':
        lysis_type = constants.ALL_SEMI_CLEAR_LYSIS

    # Couples to analyse will be in list_couples_lysis_type
    list_couples_lysis_type = []
    list_couples_lysis_type = functions.getCouplesLysis(lysis_type)
    
    # Defining two correlation tables between phages and bacteria
    phages = []
    bacterium = []

    for couple in list_couples_lysis_type:
        phages.append(BacteriophageJson.getByID(couple.bacteriophage).designation)
        # Get the name of bacteria (strain designation + species designation)
        strain_id = BacteriumJson.getByID(couple.bacterium).strain
        strain_designation = StrainJson.getByID(strain_id).designation
        specie_designation = SpecieJson.getByID(StrainJson.getByID(strain_id).specie).designation
        bacterium.append(specie_designation + '-' +  strain_designation + '\n' + str(couple.bacterium))

    # Draw network graph
    draw_graph(phages, bacterium, list_couples_lysis_type, graph_name='graph', is_png=False)

def draw_graph(phages:list, bacteria:list, list_couples_lysis_type:list,
               is_png=False,
               node_size=300, node_alpha=0.5,
               node_text_size=8,
               edge_alpha=0.5, edge_tickness=0.5,
               edge_text_pos=1.0,
               text_font='sans-serif',
               graph_name='network_graphic'):
    """
    | Make a graph to bind phages with bacteria according to their lysis attribute

    :param phages: list of phages
    :param bacteria: list of bacteria
    :param list_couple_lysis_type: list of couples
    :param is_png: if true, save the graph as a png image
    :param node_size: size of node
    :param node_alpha: clearness of node_alpha
    :param node_text_size: size of text in node
    :param edge_alpha: clearness of node_alpha
    :param edge_tickness: tichness of edge,
    :param edge_text_pos: position of text
    :param text_font: text font
    :param graph_name: name of graph

    :type phages: list
    :type bacteria: list
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
    # Plot declaration
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.set_title('Network between phages and bacteria', fontsize=16)
    margin=0.1
    fig.subplots_adjust(margin, margin, 1.-margin, 1.-margin)
    ax.axis('equal')
    
    nodes_phages = []  # All different phages
    nodes_bacteria = []  # All different bacteria
    # Different couples in function of their taxonomy
    nodes_couples_strain_level  = []  
    nodes_couples_species_level = []

    # All species of the current research
    all_species = []

    # Get the name of each bacteria (strain + species)
    for couple in list_couples_lysis_type:
        strain_id = BacteriumJson.getByID(couple.bacterium).strain
        strain_designation = StrainJson.getByID(strain_id).designation
        specie_designation = SpecieJson.getByID(StrainJson.getByID(strain_id).specie).designation
        bacteria_designation = specie_designation + '-' +  strain_designation + '\n' + str(couple.bacterium)
       
        # Get bacteria designation
        if not bacteria_designation in nodes_bacteria:
            nodes_bacteria.append(bacteria_designation)
            
        # Get phages' designation
        phages_designation = BacteriophageJson.getByID(couple.bacteriophage).designation
        if not phages_designation in nodes_phages:
            nodes_phages.append(phages_designation)

        if couple.level == constants.STRAIN_ID:
            if not phages_designation in nodes_couples_strain_level:
                nodes_couples_strain_level.append(phages_designation)
        elif couple.level == constants.SPECIES_ID:
            if not phages_designation in nodes_couples_species_level:
                nodes_couples_species_level.append(phages_designation)

        all_species.append(specie_designation)
   
    designation_of_species, number_of_species = np.unique(all_species, return_counts=True)
    list_of_list = [[] for i in range(len(number_of_species))]

    i = 0
    while(i < len(number_of_species)):
        for bact in nodes_bacteria:
            if bact.split('-')[0] == designation_of_species[i]:
                list_of_list[i].append(bact)
        i += 1
     
    nodes = set(nodes_phages + nodes_bacteria)  # All the nodes in our graph
    G=nx.Graph()  # Create networkx graph

    # Add nodes
    for node in nodes:
        G.add_node(node)

    # Add edges
    i = 0
    while(i < len(phages)):
        G.add_edge(phages[i], bacteria[i])
        i += 1

    graph_pos=nx.spring_layout(G)  # Draw graph
    # Defining nodes features for couples level strain
    nx.draw_networkx_nodes(G,graph_pos,nodelist=nodes_couples_strain_level,node_size=node_size, 
                           alpha=node_alpha, node_color='g')
    
    # Defining nodes features for couples level sepcies
    nx.draw_networkx_nodes(G,graph_pos,nodelist=nodes_couples_species_level,node_size=node_size, 
                            alpha=node_alpha, node_color='black')
    
    # Different colors for different strains
    color = ['red', 'purple', 'blue', 'orange', 'grey']
    i = 0
    for el in list_of_list:
        # Defining nodes features for bacteria
        nx.draw_networkx_nodes(G,graph_pos,nodelist=el,node_size=node_size, 
                            alpha=node_alpha, node_color=color[i])
        i = (i + 1) % 5
        
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
                           alpha=edge_alpha,edge_color='b')
   
    #display ID of bacteria and phages
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,
                            font_family=text_font)

    #show graph
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel('Rouge/Violet/Bleu/Orange/Gris = Bactéries' + 
                    ' ----- Vert = Phages - Couple niveau souche' + 
                    ' ----- Noir = Phages - Couple niveau espèce' +
                    '\nNombre de phages différents : ' + str(len(nodes_phages)) +
                    ' ----- Nombre de bactéries différentes : ' + str(len(nodes_bacteria)) +
                    '\nNombre d\'espèces différentes : ' + str(len(number_of_species)) + 
                    '\n'+ str(designation_of_species))
    
    # Save graph in png or display it
    if is_png:
        plt.savefig('./' + graph_name + '.png')
    else:
        plt.show()

def getCouplesTaxonomie(couples_list:list):
    """
    | Get couples dictionnary according to their taxonomy.

    :param couples_list: list of couples to sort

    :type couples_list: list

    :return: list of couples
    :rtype: list
    """

    taxonomie_dictionnary = {constants.STRAIN_ID: 0,
                             constants.SPECIES_ID: 0, 
                             constants.GENUS_ID: 0, 
                             constants.FAMILY_ID: 0}
    for couples in couples_list:
        if couples.level in taxonomie_dictionnary:
            taxonomie_dictionnary[couples.level] += 1

    return taxonomie_dictionnary

def getCouplesInteraction(couples:list, interaction_type:bool):
    """
    | Get couples according to their interaction.

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
    | Get accession number of all organisms in a list.

    :param organism_list: list of organisms

    :type organism_list: list

    :return: list of organism
    :rtype: list
    """
    accessions_number_list = []
    for organism in organism_list:
        accessions_number_list.append(organism.acc_number)
    
    return accessions_number_list