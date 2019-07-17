
import pandas as pd
import numpy as np
import os
import constants

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.OrganismJ import OrganismJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.StrainJ import StrainJson
from objects_API.SpecieJ import SpecieJson
from objects_API.ProteinJ import ProteinJson
from objects_API.WholeDNAJ import WholeDNAJson
from objects_API.ContigJ import ContigJson
import network_chart_couples_cl as network
import time
import matrix_similarity as ms
from Bio import Entrez
from Bio import Align 
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat.MatrixInfo import blosum62

from skbio.alignment import StripedSmithWaterman 
from skbio.alignment import AlignmentStructure

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()
Entrez.email = "christophe.joyet@heig-vd.ch"

def getbestprot(organism_id_list:list, similarity_min=0.70, similarity_max=1.0):
    """
    
    """
    protein_list = []
    matrix = []

    #récupérer les protéines des phages
    for organism in organism_id_list:
        # get all the protein of the organisms
        protein_list.append(ProteinJson.getByOrganismID(organism))
    
    first_list_of_list_of_couple = [[] for i in range (len(organism_id_list) - 1)]

    # récupérer les couples 
    for protein_phage_1 in protein_list[0]:
        for i in range (1, len(organism_id_list)):
            for protein_phage_n in protein_list[i]:
                #si les on a une correspondance de plus d'un certain pourcentage 
                similarity_score = ms.getSimilarityScoreTwoProteinLocalAlign(protein_phage_1, protein_phage_n)
                if similarity_score >= similarity_min and similarity_score <= similarity_max:
                    # creation d'un couple
                    couple = [protein_phage_1.sequence_AA, protein_phage_n.sequence_AA]
                    first_list_of_list_of_couple[i-1].append(couple)
    
    matrix.append(first_list_of_list_of_couple)

    print(len(matrix[0]))
    for i in range (len(organism_id_list) - 2):
        print("BOU")
        matrix.append(foo(matrix[i], 4, similarity_min=0.70, similarity_max=1.0))


    '''
    first_list_of_list_of_couple = [[] for i in range (len(organism_id_list) - 1)]

    # Comparaison de chaque protéines du premier phage avec les autres
    for protein_phage_1 in protein_list[0]:
        for i in range (1, len(organism_id_list)):
            for protein_phage_n in protein_list[i]:
                #si les on a une correspondance de plus d'un certain pourcentage 
                similarity_score = ms.getSimilarityScoreTwoProteinLocalAlign(protein_phage_1, protein_phage_n)
                if similarity_score >= similarity_min and similarity_score <= similarity_max:
                    # creation d'un couple
                    couple = [protein_phage_1.sequence_AA, protein_phage_n.sequence_AA]
                    first_list_of_list_of_couple[i-1].append(couple)

    # [[(phage1_prot1, phage2_protX), (phage1_prot2, phage2_protV),...], [(phage1_prot1, phage3_protY), (phage1_prot2, phage3_protK), ...]]

    
    second_list_of_list_of_couple = [[] for i in range (len(organism_id_list) - 2)]

    for couple in first_list_of_list_of_couple[0]:
        for i in range (1, len(first_list_of_list_of_couple)):
            for j in range (len(first_list_of_list_of_couple[i])):
                if couple[0] == first_list_of_list_of_couple[i][j][0]: # trouver les séquences du phages_1 identiques dans les deux couples des listes
                    similarity_score = ms.getSimilarityScoreTwoProteinLocalAlignText(couple[1], first_list_of_list_of_couple[i][j][1])
                    if similarity_score >= similarity_min and similarity_score <= similarity_max:
                        # creation d'un couple
                        couple = [couple[0], couple[1], first_list_of_list_of_couple[i][j][1]]
                        second_list_of_list_of_couple[i-1].append(couple)
    '''
    
    print("tot prot phage 1 : " + str(len(protein_list[0])))
    print("tot prot phage 2 : " + str(len(protein_list[1])))
    print("tot prot phage 3 : " + str(len(protein_list[2])))
    '''print("tot couple phage 1 - 2 : " + str(len(first_list_of_list_of_couple[0])))
    print("tot couple phage 1 - 3 : " + str(len(first_list_of_list_of_couple[1])))
    print("tot couple phage 1,2,3 : " + str(len(second_list_of_list_of_couple[0])))'''

def foo(first_list_of_list_of_couple:list, nombre_organismes:int, similarity_min=0.70, similarity_max=1.0):

    second_list_of_list_of_couple = [[] for i in range (len(first_list_of_list_of_couple) - 1)]

    for couple in first_list_of_list_of_couple[0]:
        for i in range (1, len(first_list_of_list_of_couple)):
            for j in range (len(first_list_of_list_of_couple[i])):
                    if couple[nombre_organismes - len(first_list_of_list_of_couple)] == first_list_of_list_of_couple[i][j][nombre_organismes - len(first_list_of_list_of_couple)]: # trouver les séquences du phages_1 identiques dans les deux couples des listes
                        similarity_score = ms.getSimilarityScoreTwoProteinLocalAlignText(couple[nombre_organismes - len(first_list_of_list_of_couple) + 1], first_list_of_list_of_couple[i][j][nombre_organismes - len(first_list_of_list_of_couple) + 1])
                        if similarity_score >= similarity_min and similarity_score <= similarity_max:
                            # creation d'un couple
                            couple = [couple, first_list_of_list_of_couple[i][j][nombre_organismes - len(first_list_of_list_of_couple) + 1]]
                            second_list_of_list_of_couple[i-1].append(couple)
    
    return second_list_of_list_of_couple

getbestprot([5038,5030,5023,5045])