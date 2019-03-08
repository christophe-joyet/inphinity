# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 13:52:10 2018

@author: Diogo Leite
"""



from files_treatment_new.fasta_nucleotid_NCBI import *
from files_treatment_new.fasta_protein_NCBI import *
from files_treatment_new.fasta_whole_genome_NCBI import *

from objects_new.WholeDNA_new import * 

from objects_new.Gene_new import * 
from objects_new.Proteins_new import * 


id_organism = 8434 
id_whole_genome = 8434  

path_file_prots_aa = 'C:/Users/Stage/Desktop/xavier' + '/NC_026600.1_PROTEIN_CS.fasta'
path_file_prots_nuc =  'C:/Users/Stage/Desktop/xavier' + '/NC_005880.2_DNA_CS.fasta'
path_file_prots_WG  =  'C:/Users/Stage/Desktop/xavier' + '/NC_005880.2_COMPLETE.fasta'


#fasta_nuc = Fasta_nucleotid_NCBI(path_file_prots_nuc)
fasta_prot = Fasta_protein_NCBI(path_file_prots_aa)
#fasta_whole_dna = Fasta_whole_genome_NCBI(path_file_prots_WG)

#whole_genome_obj = fasta_whole_dna.get_whole_genome()
#whole_genome_obj.id_wholeDNA = id_whole_genome
#print(whole_genome_obj)
#a=0


#protein treatment
aux = 0
list_prots = fasta_prot.get_list_of_proteins()
#for protein in list_prots:
#    sequence_aa = fasta_nuc.get_nucleotid_sequence_prot_id(protein.id_accession)
#    list_prots[aux].sequence_dna = sequence_aa
#    a +=1

#update whole_genome
#whole_genome_obj.update_whole_dna_by_id(id_whole_genome)

#creation of proteins and genes
for protein in list_prots:
    gene_obj = Gene(FK_id_organism = id_organism)
    id_prot = protein.create_protein()
    gene_obj.FK_id_protein = id_prot
    gene_obj.create_gene()


