# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 09:56:10 2017

@author: Stage
"""


from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
import yaml

import time

from objects_new.Families_new import *
from objects_new.Genus_new import *
from objects_new.Species_new import *
from objects_new.Strains_new import *

from objects_new.WholeDNA_new import *
from objects_new.Organisms_new import *
from objects_new.Contigs_new import *

from objects_new.Proteins_new import *
from objects_new.Gene_new import *

from files_treatment.fasta_parsing import *

import os

import pymysql
pymysql.install_as_MySQLdb()


list_files_protein = os.listdir('./sequence_files_phage_interact/Fasta_AA')
list_files_nc = os.listdir('./sequence_files_phage_interact/Fasta')
list_files_gnk = os.listdir('./sequence_files_phage_interact/GenBank')

cwd = os.getcwd()
print(cwd)
cwd_fasta_protein = cwd + "/sequence_files_phage_interact/Fasta_AA/"
cwd_fasta = cwd + "/sequence_files_phage_interact/Fasta/"
cwd_genBank = cwd + "/sequence_files_phage_interact/GenBank/"
print("---------------")
print(len(list_files_protein))

aux_index_files = 0


insert_type = False
insert_organism = False

informations_json = None



def get_aa_sequence(dict_values_proteins, dict_values_nucleotids, value):
#    print(dict_values_proteins.keys())
#    print(dict_values_nucleotids['fig|470.1122.peg.1'])
    key_value = (list(dict_values_proteins.keys())[list(dict_values_proteins.values()).index(value)]) # Prints george
    nucleotid_sequence = dict_values_nucleotids[key_value]
    return nucleotid_sequence


def parse_organism(file_path_name, dict_prot, dict_nn):
#%%Insert taxonomy Bacterim
    #informations_json = yaml.load(repr(record.annotations))
    #taxonomy = informations_json['taxonomy']
    #print(len(taxonomy))
    #organisme_name = informations_json['organism']
    
    ########### For bacterium ###########
    #print(file_path_name)
    #strain = file_path_name.lower().split("streptococcus_tigurinus_")[1]
    #strain = strain[:-4]
    
    
    ########### For phage ###########
    strain = file_path_name.lower().split("Staphylococcus_phage_")[1]
    strain = strain[:-4]
    gi_acc_design = strain.replace("_", " ")
    strain = "Staphylococcus phage" + strain
    
    #assert len(taxonomy) == 7
    family_obj = Family(designation='Phage no family')
    genus_obj = Genus(designation='Phage no genuse')
    specie_obj = Specie(designation='Phage no Specie')
    strain_obj = Strain(designation=strain)
    print("Family : " + str(family_obj.designation))
    print("Genus : " + str(genus_obj.designation))
    print("Specie: " + str(specie_obj.designation))
    print("Strain: " + str(strain_obj.designation))
    

    #Taxonomy insertion
    id_family = family_obj.create_family(family_obj.designation)
    id_genus = genus_obj.create_genus(genus_obj.designation, id_family)
    id_specie = specie_obj.create_specie(specie_obj.designation, id_genus)
    id_strain = strain_obj.create_strain(strain_obj.designation, id_specie)   
    
    
    #%% Whole_Genome
    whole_genome_obj = WholeDNA(head = "NA", head_id = "NA", sequence="NA")
    id_whole_genome = whole_genome_obj.create_whole_dna_no_verification("NA","NA", "NA")
    
    #%% Organism     
    print(id_strain)
    print(id_whole_genome)
    time.sleep(5)
    print('***-*-*-*-*-*-*-*-*****')
    organism_obj = Organism(-1, 'Greg' + gi_acc_design, 'Greg' + gi_acc_design, -1, True, -1, 2, id_strain, 2, id_whole_genome, 3)
    id_organism = organism_obj.create_organism()
    print(id_organism)
    #(id_organism, gi, acc_num, qty_proteins, assembled, qty_contig, fk_source = -1, fk_strain = -1, fk_type = -1, fk_whole_genome = -1
    #Insert taxonomy

    #inserer famille, genus, strain.
    list_of_proteins = []
#%% Protein treatment
    aux = 0
    for protein_key in dict_prot.keys():
        prot_seq = dict_prot[protein_key]
        prot_design = protein_key
        
        nucleotid_sequence = get_aa_sequence(dict_prot, dict_nn, prot_seq)
        protein_obj = Protein(-1, -1, prot_design, str(prot_seq), str(nucleotid_sequence), -1, -1, -1, -1, 1)
        list_of_proteins.append(protein_obj)
    print("Proteins: " + str(len(list_of_proteins)))
    for protein in list_of_proteins:
        id_protein = protein.create_protein()
        gene_obj = Gene(FK_id_organism=id_organism, FK_id_protein=id_protein)
        value_id_gene = gene_obj.create_gene()
        
        aux += 1
        if aux % 1500 == 0:
            print("sleep")
            time.sleep(3)
            print("End sleep")
    
#%% main function

list_files_name = []
bact_count = 0
for file_gen_name in list_files_protein:
    #For bacterium
    #if "phage" not in file_gen_name.lower() and "phi" not in file_gen_name.lower() and "Acinetobacter_bacteriophage_" in file_gen_name.lower():
    #For Phage
    if "Staphylococcus_phage_" in file_gen_name.lower():
        print("Hello")

        
        insert_type = False

        insert_organism = False

        
        informations_json = None

        id_strain = -1
        id_organism = -1
        id_whole_genome = -1        
        file_proteins = list_files_protein[aux_index_files]
        file_nc = file_proteins[0:-4] + '.fna'
        file_genBank = file_proteins[0:-4] + '.gbk'
        
        print("Files")
        print(file_proteins)
        print(file_nc)
        print(file_genBank)
        
        fasta_file_protein = Fasta_parsing(cwd_fasta_protein + file_proteins)
        dict_values_proteins = fasta_file_protein.parse_fasta()
        
        fasta_file = Fasta_parsing(cwd_fasta + file_nc)
        dict_values = fasta_file.parse_fasta()
        
        fasta_file_path = cwd_fasta_protein + file_proteins
        
        list_files_name.append(fasta_file_path)
        bact_count += 1
        
        
        print(file_proteins)
        print(fasta_file)
        print(file_genBank)
        
        print(len(dict_values_proteins))
        print(len(dict_values))
        
        print(fasta_file_protein)
        parse_organism(fasta_file_path, dict_values_proteins, dict_values)
        #print("Hello")
        
        #print("NON")
    aux_index_files += 1

print(bact_count)