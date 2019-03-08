# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 12:00:13 2017

@author: Stage
"""

# DONT FORGET TO UNCOMMENT EVERYTHING...
# DONT FORGET TO UNCOMMENT EVERYTHING...
# DONT FORGET TO UNCOMMENT EVERYTHING...
# DONT FORGET TO UNCOMMENT EVERYTHING...
# DONT FORGET TO UNCOMMENT EVERYTHING...

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

from objects_new.Contigs_new import *

from files_treatment.fasta_parsing import *

import os

import pymysql
pymysql.install_as_MySQLdb()


list_files_gnk = os.listdir('./sequence_files/GenBank')

cwd = os.getcwd()
print(cwd)

cwd_genBank = cwd + "/sequence_files/GenBank/"
print("---------------")

aux_index_files = 0


insert_type = False
insert_organism = False

informations_json = None



#%%Function parse file
def parse_organism(gnb_path_name_file, genus, specie, strain):
#%% checke if the organism exists
    organism_obj = Organism(id_organism = -1, gi = "", acc_num = "", qty_proteins = 150, assembled = False, qty_contig = 5)    
    id_organism = organism_obj.get_id_organism_by_taxonomy(genus, specie, strain)
    id_fk_whole_dna = organism_obj.get_id_fk_whole_dna_by_id(id_organism)
    print(id_organism)
    print(id_fk_whole_dna)
    
    
    if id_organism != -1:
        with open(gnb_path_name_file, "rU") as input_handle:
            for record in SeqIO.parse(input_handle, "genbank"):
                print("------------")
                print(record.seq[0:10])
                #print(record.features[0])
                print(len(record.features))
                
                
                
                #Insert contig
                contig_obj = Contig(sequence = record.seq, fk_id_whole_genome=id_fk_whole_dna)
                id_contig = contig_obj.create_contig()
                
                
                if len(record.features) > 0:
                    for feature in record.features:
                        if 'db_xref' in feature.qualifiers:
                            designation_prot = feature.qualifiers['db_xref'][0].replace('SEED:','')
                            prot_obj = Protein(designation = designation_prot)
                            id_prot = prot_obj.get_id_prot_by_designation()
                            if id_prot != -1:
                                print(id_prot)
                                print(feature.qualifiers['db_xref'][0].replace('SEED:',''))
                                print(feature.location.start)
                                print(feature.location.end)
                                
                                
                                prot_obj.id_protein = id_prot
                                prot_obj.start_point_cnt = int(feature.location.start)
                                prot_obj.end_point_cnt = int(feature.location.end)
                                prot_obj.fk_id_contig = id_contig
                                prot_obj.update_protein_contig()
                    #print(record)
                    print("------------")
#    global insert_type
#    global insert_organism
#    global informations_json
#    for record in SeqIO.parse(gnb_path_name_file, "genbank"):
#        print(gnb_path_name_file)
#        try:
#            json_features = record.features
#            informations_json = yaml.load(repr(record.annotations))
#            #%% Create the organism        
#            if insert_organism == False:
#                #informations_json = yaml.load(repr(record.annotations))
#                #taxonomy = informations_json['taxonomy']
#                #print(len(taxonomy))
#                #organisme_name = informations_json['organism']
#                strain = gnb_path_name_file.split("Acinetobacter_baumannii_")[1]
#                strain = strain[:-4]
#                #assert len(taxonomy) == 7
#                family_obj = Family(designation='Enterobacteriaceae')
#                genus_obj = Genus(designation='	Escherichia')
#                specie_obj = Specie(designation='E. coli')
#                strain_obj = Strain(designation=strain)
#                print("Family : " + str(family_obj.designation))
#                print("Genus : " + str(genus_obj.designation))
#                print("Specie: " + str(specie_obj.designation))
#                print("Strain: " + str(strain_obj.designation))
#                
#                
#                response = input("The information seems correct (y/n) ")
#                if response == 'n':
#                    break
#                #%% Taxonomy        
#    
#                #Taxonomy insertion
#                id_family = family_obj.create_family(family_obj.designation)
#                id_genus = genus_obj.create_genus(genus_obj.designation, id_family)
#                id_specie = specie_obj.create_specie(specie_obj.designation, id_genus)
#                id_strain = strain_obj.create_strain(strain_obj.designation, id_specie)
#                
#                
#                
#    
#                
#                #%% Whole_Genome
#                whole_genome_obj = WholeDNA(head = "NA", head_id = "NA", sequence="NA")
#                id_whole_genome = whole_genome_obj.create_whole_dna("NA","NA", "NA")
#                
#                
#                #%% Organism     
#                print(id_strain)
#                print(id_whole_genome)
#                print('***-*-*-*-*-*-*-*-*****')
#                organism_obj = Organism(-1, 'Greg' + strain_obj.designation, -1, -1, True, -1, 2, id_strain, 1, id_whole_genome)
#                id_organism = organism_obj.create_organism()
#                #(id_organism, gi, acc_num, qty_proteins, assembled, qty_contig, fk_source = -1, fk_strain = -1, fk_type = -1, fk_whole_genome = -1
#                #Insert taxonomy
#            
#                #inserer famille, genus, strain.
#                insert_organism = True
#                
#                
#            num_proteins = 0
#            sequence_contig = record.seq
#            contig_obj = Contig(fk_id_whole_genome = id_whole_genome, sequence = str(sequence_contig))
#            id_contig = contig_obj.create_contig()
#            
#            #print(strain)
#            print('----------------------')
#            
#            
#            
#    #%% proteins
#            list_of_proteins = []
#            #print(len(record.features))
#            #print(record.features[2].qualifiers['translation'])
#            for protein in record.features :
#                if 'translation' in protein.qualifiers.keys():
#                    
#                    #self, id_protein = -1, id_prot_DB_online = -1, designation = "", sequence_prot = "", sequence_dna = "", start_point = -1, end_point = -1, start_point_cnt = -1, end_point_cnt = -1, fk_id_contig= -1
#                    
#                    num_proteins += 1
#            
#                    #print('+++++-----+++')
#                    #print(protein.qualifiers)
#                    #print('++++++++++++')
#                    #print(protein.qualifiers['translation'])
#                    #print('++++++++++++')
#                    prot_id = protein.qualifiers['db_xref'][0]
#                    prot_seq = protein.qualifiers['translation'][0]
#                    prot_cnt_start = protein.location.start
#                    prot_cnt_end = protein.location.end
#                    prot_designation = protein.qualifiers['product'][0]
#                    
#                    
#                    print(len(dict_prot))
#                    print(len(dict_nn))
#                    
#                    nucleotid_sequence = get_aa_sequence(dict_prot, dict_nn, prot_seq)
#                    print("---Creating proteins---")
#                    protein_obj = Protein(-1, -1, prot_designation, str(prot_seq), str(nucleotid_sequence), -1, -1, prot_cnt_start, prot_cnt_end, id_contig)
#                    list_of_proteins.append(protein_obj)
#                    #print(prot_seq)
#            print(num_proteins)
#            print(len(list_of_proteins))
#            
#            aux = 0
#            print("---Insertion of genes---")
#            for protein in list_of_proteins:
#                id_protein = protein.create_protein()
#                print(id_protein)
#                gene_obj = Gene(FK_id_organism=id_organism, FK_id_protein=id_protein)
#                value_id_gene = gene_obj.create_gene()
#                
#                aux += 1
#                if aux % 1500 == 0:
#                    print("sleep")
#                    time.sleep(3)
#    
#    
#            
#    
#            
#    
#    
#    #%% Exception zone        
#        except AssertionError:
#            print('La taxonomie est plus detaillee que ce qui est pretendu')
#            
#    
#            #print(repr(record.annotations))
#            #aux = repr(record.annotations)
#            #print(aux)
#            
#            
#            
#            #Info de l'organisme
#            #aux = repr(record.annotations)
#            #auxB = d = yaml.load(aux)
#            #print(type(auxB))
#            #Obtenir la taxonomie de l'organism
#            #print(auxB['taxonomy'])
#            #oubtenir la source
#            #print(auxB['source'])
#            #obtenir le nom de l'organism
#            #print(auxB['organism'])
#        #except KeyError:
#         #   print('Taxonomy doesn\'t exist')
#        
     


#%% for each organism (Bacterium)
bact_count = 0
list_files_name = []
list_files_errors = []
for file_gen_name in list_files_gnk:
    
#%% bacteria part    
    if "phage" not in file_gen_name.lower() and "phi" not in file_gen_name.lower():
        try:       
            genbank_file_path = cwd_genBank + file_gen_name
            
            list_files_name.append(genbank_file_path)
            
            
            print(genbank_file_path)
            
            taxonomy_filne_name = file_gen_name[0:-4].split("_")
            
            assert len(taxonomy_filne_name) == 3
            
               
            genus = taxonomy_filne_name[0]
            specie = taxonomy_filne_name[1]
            strain = taxonomy_filne_name[2]
            
            print("Taxonomy: ")
            print("genus: " + genus)
            print("Specie: " + specie)
            print("Strain: " + strain)
            
            parse_organism(genbank_file_path, genus, specie, strain)

            bact_count += 1
        except AssertionError:
            print("It is not possible to obtained the taxonomy based on it name")
            print(file_gen_name)
            list_files_errors.append(file_gen_name)

        
        

        #print("Hello")
        
        #print("NON")
    aux_index_files += 1

print(bact_count)
#
#insert_type = False
#insert_organism = False
#
#informations_json = None
#
#id_strain = -1
#id_organism = -1
#id_whole_genome = -1
#
#path_protein_fasta_file = "83333.257.faa"
#path_nucleotid_fasta_file = "83333.257.fna"





#fasta_obj_proteins = Fasta_parsing(path_protein_fasta_file)
#dict_values_proteins = fasta_obj_proteins.parse_fasta()
#
#fasta_obj_nucleotids = Fasta_parsing(path_nucleotid_fasta_file)
#dict_values_nucleotids = fasta_obj_nucleotids.parse_fasta()



   
