# -*- coding: utf-8 -*-
"""
Created on Wen Fev 27 15:01:05 2018

@author: Diogo
"""

from objects_new.Organisms_new import *
from objects_new.Strains_new import *
from objects_new.Species_new import *
from objects_new.Genus_new import *
from objects_new.Families_new import *

from objects_new.WholeDNA_new import *
from objects_new.Contigs_new import *
from objects_new.Proteins_new import *
from objects_new.Gene_new import *

list_ids_organisms = Organism.get_ids_all_organisms()
print(list_ids_organisms)

print(len(list_ids_organisms))

for id_organism in list_ids_organisms[10:]:
    print("------------------")
    #Get the organism
    organims_obj = Organism.get_organism_by_id(id_organism)
    print(organims_obj)

    #TAXONOMY
    strain_obj = Strain.get_strain_by_id(organims_obj.fk_strain)
    print(strain_obj)

    specie_obj = Specie.get_specie_by_id(strain_obj.fk_specie)
    print(specie_obj)

    genus_obj = Genus.get_genus_by_id(specie_obj.fk_genus)
    print(genus_obj)

    family_obj = Family.get_family_by_id(genus_obj.fk_family)
    print(family_obj)

    #get Whole_DNA
    whole_dna_obj = WholeDNA.get_whole_dna_by_id(organims_obj.fk_whole_genome)
    print(whole_dna_obj)

    #All contigs
    list_contigs_objs = Contig.get_all_Contigs_bi_whole_DNA_id(whole_dna_obj.id_wholeDNA)
    print(len(list_contigs_objs))

    #Get all proteins
    for contig_obj in list_contigs_objs:
        list_proteins = Protein.get_all_Proteins_by_fk_contig(contig_obj.id_contig)
        qty_prots = len(list_proteins)
        if qty_prots == 0:
            print(len(list_proteins))


    ##Insertion dans nouvelle base de données - CHANGER DANS FACTORY
    ##taxonomy
    id_family = family_obj.create_family()
    genus_obj.fk_family = id_family
    id_genus = genus_obj.create_genus()
    specie_obj.fk_genus = id_genus
    id_specie = specie_obj.create_specie()
    strain_obj.fk_specie = id_specie
    id_strain = strain_obj.create_strain()

    organims_obj.fk_strain = id_strain

    #whole_genome

    id_whole_genome = whole_dna_obj.create_whole_dna_no_verification()
    organims_obj.fk_whole_genome = id_whole_genome

    id_organism = organims_obj.create_organism()

    #CHANGER DB POUR RéCUPéRER les PROTEINS
    for contig_obj in list_contigs_objs:
        contig_obj.fk_id_whole_genome = id_whole_genome

        list_proteins = Protein.get_all_Proteins_by_fk_contig(contig_obj.id_contig)

        #CHANGER DB
        id_contig = contig_obj.create_contig_no_verification()
        for protein in list_proteins:
            protein.fk_id_contig = id_contig
            id_protein = protein.create_protein_all_details()

            gene_obj = Gene(gene_number = -1, dna_head = "No head", 
                        dna_sequence = "No sequence", start_position = -1, 
                        end_position = -1, FK_id_organism = id_organism, FK_id_protein = id_protein)

            gene_obj.create_gene()
   


