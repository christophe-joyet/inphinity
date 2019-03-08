# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 08:57:47 2018

@author: Diogo
"""

# NOTE: this script does not take in consideration the data for the datasets


from objects_new.Couples_new import *
from objects_new.Proteins_new import *
from objects_new.Gene_new import *
from objects_new.Protein_dom_new import *
from objects_new.Contigs_new import *
from objects_new.Organisms_new import *
from objects_new.WholeDNA_new import *
from objects_new.Strains_new import *

organism_id = 10608

#Organism
organism_obj = Organism.get_organism_by_id(organism_id)
print(organism_obj)

#couples
list_couples = Couple.get_all_couples_by_bacterium(organism_id)
print(len(list_couples))

#list genes
list_genes = Gene.get_all_Genes_by_organism_id(organism_id)
print(len(list_genes))


#list prots
list_proteins = Protein.get_all_Proteins_by_organism_id(organism_id)
print(len(list_proteins))

#list_prot_doms
try:
    list_prot_doms = ProteinDom.get_all_protein_domain_by_protein_id(list_proteins[0].id_protein)
    print(len(list_prot_doms))
except:
    print("no proteins")

#list contigs
try:
    list_contig = Contig.get_all_Contigs_bi_whole_DNA_id(organism_obj.fk_whole_genome)
    print(len(list_contig))
except:
    print("no whole contigs")

#whole dna
try:
    whole_dna_obj = WholeDNA.get_whole_dna_by_id(organism_obj.fk_whole_genome)
    print(whole_dna_obj)
except:
    print("No whole dna")


########### delete part ############
Couple.remove_couple_by_id_bacterium(organism_obj.id_organism)
Gene.delete_gene_from_id_organism(organism_obj.id_organism)
for protein in list_proteins:
    ProteinDom.remove_prot_dom_by_protein_id(protein.id_protein)
    Protein.remove_protein_by_its_id(protein.id_protein)

Contig.remove_contig_by_FK_whole_dna(whole_dna_obj.id_wholeDNA)
Organism.remove_organism_by_id(organism_id)
Strain.remove_strain_by_id(organism_obj.fk_strain)
WholeDNA.remove_whole_dna_by_id(whole_dna_obj.id_wholeDNA)





