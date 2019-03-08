# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 15:04:50 2017

@author: Stage
"""

from files_treatment.fasta_proteins import *
from files_treatment.fasta_parsing import *
from files_treatment.genBank_file import *

from objects_new.Families_new import *
from objects_new.Genus_new import *
from objects_new.Species_new import *
from objects_new.Strains_new import *

from objects_new.WholeDNA_new import *

from objects_new.Organisms_new import *
from objects_new.Gene_new import *


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

#Here i parse fasta files for phages

def get_aa_sequence(dict_values_proteins, dict_values_nucleotids, value):
#    print(dict_values_proteins.keys())
#    print(dict_values_nucleotids['fig|470.1122.peg.1'])
    key_value = (list(dict_values_proteins.keys())[list(dict_values_proteins.values()).index(value)]) # Prints george
    nucleotid_sequence = dict_values_nucleotids[key_value]
    return nucleotid_sequence

fasta_parsin_obj = Fasta_parsing("")


for file_prot, file_nuc, file_genbank in zip(list_files_protein, list_files_nc, list_files_gnk):
    if "staphylococcus_phage_" in file_prot.lower():
        print(file_prot)

        fasta_parsin_obj.path_file = (cwd_fasta_protein + file_prot)
        fasta_prot = fasta_parsin_obj.parse_fasta_all_informations()
        
        fasta_parsin_obj.path_file = (cwd_fasta + file_nuc)
        fasta_nuc = fasta_parsin_obj.parse_fasta_all_informations()
        
        print(type(fasta_prot))
        
        
        fasta_obj = proteins_fasta(fasta_prot, fasta_nuc)
        list_prot = fasta_obj.parse_fasta_format()
        
        
        
        
        print(list_prot[0].id_protein)
        print(list_prot[0].id_prot_DB_online)
        print(list_prot[0].designation)
        print(str(list_prot[0].sequence_dna))
        print(str(list_prot[0].sequence_prot))
        print(list_prot[0].start_point)
        print(list_prot[0].end_point)



        
        

        
        id_strain = -1
        id_organism = -1
        id_whole_genome = -1        
        file_proteins = file_prot
        file_nc = file_proteins[0:-4] + '.fna'
        file_genBank = file_proteins[0:-4] + '.gbk'
        
        strain = file_prot.split("_")[2][0:-4]
        
#%% taxonomy zone
        family_obj = Family(designation='Phage no family')
        genus_obj = Genus(designation='Phage no genuse')
        specie_obj = Specie(designation='Phage no Specie')
        strain_obj = Strain(designation="Staphylococcus Phage " + strain)
        print("Family : " + str(family_obj.designation))
        print("Genus : " + str(genus_obj.designation))
        print("Specie: " + str(specie_obj.designation))
        print("Strain: " + str(strain_obj.designation))
        
        id_family = family_obj.create_family(family_obj.designation)
        id_genus = genus_obj.create_genus(genus_obj.designation, id_family)
        id_specie = specie_obj.create_specie(specie_obj.designation, id_genus)
        id_strain = strain_obj.create_strain(strain_obj.designation, id_specie)   
        
#%%Whole genome zone
        whole_genome_obj = WholeDNA(head = "NA", head_id = "NA", sequence="NA")
        id_whole_genome = whole_genome_obj.create_whole_dna_no_verification("NA","NA", "NA")
   
     
#%% Organism zone
        genBank_obj = gen_bank_parsing_file(cwd_genBank + file_genbank)
        genBank_obj.parse_genebank_file()
    
        print(id_strain)
        print(id_whole_genome)
        print('***-*-*-*-*-*-*-*-***')
        #self, id_organism, gi, acc_num, qty_proteins, assembled, qty_contig, fk_source = -1, fk_strain = -1, fk_type = -1, fk_whole_genome = -1, fk_source_data = "NULL"
        organism_obj = Organism(-1, 'Greg' + genBank_obj.GI, genBank_obj.Acc, -1, True, -1, 2, id_strain, 2, id_whole_genome, 1)
        id_organism = organism_obj.create_organism()
        print(id_organism)
        
        print("Organism GI: %s \nOrganism Acc: %s" % (genBank_obj.GI, genBank_obj.Acc))
#%% Genes and proteins insertion
        aux = 0
        for protein in list_prot:
            id_protein = protein.create_protein()
            gene_obj = Gene(FK_id_organism=id_organism, FK_id_protein=id_protein)
            value_id_gene = gene_obj.create_gene()
            
            aux += 1
            if aux % 1500 == 0:
                print("sleep")
                time.sleep(3)
                print("End sleep")
        


    

