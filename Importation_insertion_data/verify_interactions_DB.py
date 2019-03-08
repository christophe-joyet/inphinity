# -*- coding: utf-8 -*-
"""
Created on Wen Jan 31 14:44:17 2018

@author: Diogo
"""

from Bio import SeqIO
from Bio import Entrez
import re

import os

from objects_new.Organisms_new import *
from files_treatment_new.fasta_protein_NCBI import *
from objects_new.Gene_new import *
from os.path import commonprefix
from objects_new.Protein_dom_new import *
from objects_new.PPI_preview_new import *
from objects_new.Gene_new import *
from objects_new.COGS_prot_new import *
from objects_new.F_score_blast_P_new import *
from objects_new.COGS_preview_new import *


def write_id_fasta_error(id_organism):
    file_name_db_verification ="fast_error_during_verification.csv"
    file = open(file_name_db_verification,'a') 
    file.write(id_organism)
    file.close() 

def write_id_corrected_file(id_organism, qty_prots_falsta, qty_prots_organism):
    file_name_db_verification =""
    content = ""
    if qty_prots_falsta == qty_prots_organism:
        file_name_db_verification ="id_corrects_verified.csv"
        content = str(id_organism) + " " +
    file = open(file_name_db_verification,'a') 
    file.write(id_organism)
    file.close() 


def get_all_organisms():
    list_organism_NCBI = Organism.get_ids_all_organisms_by_params(type_organism = 1, source_data = 1)
    list_organism_RAST = Organism.get_ids_all_organisms_by_params(type_organism = 1, source_data = 2)

    list_organisms = list(set(list_organism_NCBI + list_organism_RAST))
    return list_organisms

def get_prots_ids_by_organism_acc_id_seq(acc_organism):
    """
    Get fasta file and list of all proteins in the organism by its acc

    :param acc_organism: id of the cog score - -1 if unknown

    :type acc_organism: int - required

    :return: list of proteins
    :rtype: List(Proteins)
    """

    #print(os.getpid())
    file_name = 'C:\\Users\\ci4cb\\AppData\\Local\\Temp\\temp_orga_{0}.fasta'.format(os.getpid())

    list_prots = []

    Entrez.email = "diogo1790@hotmail.com"

    fasta_proteins_aa = Entrez.efetch(db="nuccore", id=acc_organism, rettype="fasta_cds_aa", retmode="text")
    fasta_content = fasta_proteins_aa.read()

    temp = open(file_name, 'w')
    temp.write(fasta_content)


    fasta_obj = Fasta_protein_NCBI(file_name)
    list_prots = fasta_obj.get_list_of_proteins()
    temp.close()


    return list_prots

def verify_proteins(list_proteins, organism_id):
    """
    remove duplicates proteins in the database

    :param list_of_proteins: list of proteins object - -1 if unknown
    :param organism_id: id of the organis - -1 if unknown

    :type list_of_proteins: List(Protein) - required
    :type organism_id: int - required

    """
    for protein in list_of_proteins:
        list_ids = Protein.get_protein_id_by_sequence_location(protein.sequence_prot, protein.start_point, protein.end_point, organism_id)
        if len(list_ids) > 1:
            print('noooooooooooooo')


def confirm_orga_in_db(organism_id):
    list_ids_duplicates_prots = Protein.get_duplicates_sequences_ids_by_organism_id(organism_id)
    print(list_ids_duplicates_prots)
    return list_ids_duplicates_prots

def confirm_organism(list_of_proteins, organism_id):
    """
    Confirm if the number of proteins in the fasta are the same that those in the database

    :param list_of_proteins: list of proteins object - -1 if unknown
    :param organism_id: id of the organis - -1 if unknown

    :type list_of_proteins: List(Protein) - required
    :type organism_id: int - required

    """
    quantity_prots = len(list_of_proteins)
    gene_number = Gene.get_number_of_Genes_by_organism_id(organism_id)
    if gene_number != quantity_prots and quantity_prots != 0:
    #if quantity_prots >= (30 + gene_number) or (gene_number - 30) >= quantity_prots:
        verify_proteins(list_of_proteins, organism_id)
        return False

    elif gene_number == quantity_prots:



def get_part_string_save(list_ids_duplicates, organism_id):
    list_ids_dup = []
    list_string_description = []
    for element in list_ids_duplicates:
        list_ids_dup.append(element[0])
        list_ids_dup.append(element[1])

    list_prots = Protein.get_all_Proteins_by_organism_id(organism_id)
    for protein in list_prots:
        if protein.id_protein not in list_ids_dup:
            list_string_description.append(protein.designation)
        else:
            print('noooo')
    common_part = commonprefix(list_string_description)
    return common_part

def delete_prots(list_duplicate, id_organism):
    for element in list_duplicate:
        id_prot = element[1]
        print(id_prot)
        PPI_preview.remove_PPI_preview_by_protein_id(id_prot)
        ProteinDom.remove_prot_dom_by_protein_id(id_prot)
        F_score_blast_P.delete_FK_score_bast_P_by_protein_id(id_prot)
        list_prot_cog_ids = COGSprot.get_all_COGSprot_by_protein_id(id_prot)
        for element in list_prot_cog_ids:
            COGS_preview.delete_COGS_preview_by_cogProt_id(element.id_COG_prot)
        COGSprot.remove_COGProt_by_protein_id(id_prot)
        Gene.delete_gene_from_id_protein(id_prot)
        Protein.remove_protein_by_its_id(id_prot)

        



#list_organisms = get_all_organisms()
#list_organisms = [4648,4651,4653,4663,4674,4683,4689,4696,4742,4793,4828,4848,4946,5021,4677,5027,5034,5043,5054,5067,5073,5110,5124,5142,5150,5170,5205,5225,5294,5314,5346,5467,5468,5518,5523,5580,5601,5614,5621,5622,5640,5652,5660,5670,5675,5718,5727,5861,5866,5883,5899,5960,5968,5983,5984,6013,6036,6091,6114,6155,6156,6186,6199,6215,6228,6295,6334,6336,6341,5904,6414,6425,6431,6472,6528,6537,6574,6593,6755,8191,8219,8256,8959,8958,8960,8957,8964,8965,9360,9359,9361,8622,8812,10093,8773,4791,9278,8666,5025,8001,8783,9638,9490,7814,7095,6219,9538,9942,8785,9657,6354,5811,8361,10079,6060,8787,8833,9748,9630,10143,9739,10113,9994,7765,10130,10036,10142,8874,8784,8744,10136,10131,10567,7832,9659,8626,10121,8085,8664,8774,9891,4718,8782,4752,7613,7720,8284,8802,8293]

list_ids_organisms_verified =[4648,4651,4653,4663,4674,4683,4689,4696,4742,4793,4828,4848,4946,5021,4677,5027,5034,5043,5054,5067,5073,5110,5124,5142,5150,5170,5205,5225,5294,5314,5346,5467,5468,5518,5523,5580,5601,5614,5621,5622,5640,5652,5660,5670,5675,5718,5727,5861,5866,5883,5899,5960,5968,5983,5984,6013,6036,6091,6114,6155,6156,6186,6199,6215,6228,6295,6334,6336,6341,5904,6414,6425,6431,6472,6528,6537,6574,6593,6755,8191,8219,8256,8959,8958,8960,8957,8964,8965,9360,9359,9361,8622,8812,10093,8773,4791,9278,8666,5025,8001,8783,9638,9490,7814,7095,6219,9538,9942,8785,9657,6354,5811,8361,10079,6060,8787,8833,9748,9630,10143,9739,10113,9994,7765,10130,10036,10142,8874,8784,8744,10136,10131,10567,7832,9659,8626,10121,8085,8664,8774,9891,4718,8782,4752,7613,7720,8284,8802,8293]

list_organisms = Organism.get_all_Organisms()


for organism in list_organisms:
    organism_id = organism.id_organism
    if organism_id not in list_ids_organisms_verified and organism_id > 4786:
        print("Id organism: {0}".format(organism_id))
        list_duplicate_ids = confirm_orga_in_db(organism_id)
        #part_description = get_part_string_save(list_duplicate_ids, organism_id)
        #list_proteins_ids, list_prots_ids_error = get_prots_ids_by_organism_acc(organims_obj.acc_num)
        organims_obj = Organism.get_organism_by_id(organism_id)
        delete_prots(list_duplicate_ids, organism_id)
        #list_of_proteins = get_prots_ids_by_organism_acc_id_seq(organims_obj.acc_num)
        #confirm_organism(list_of_proteins, organism_id)


