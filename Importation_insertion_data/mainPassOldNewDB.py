#import DAL





# %% Test Zone

##Test old bacteria
import time
from objects_old.Bacteria_old import *

from files_treatment.bacteria_files import *

from objects_new.Families_new import * 
from objects_new.Genus_new import *
from objects_new.Species_new import *
from objects_new.Strains_new import *
from objects_new.Organisms_new import *
from objects_new.Genus_new import *
from objects_new.Gene_new import *

from SQL_obj_new.Organism_sql_new import *

bact_obj = Bacteria_old()
organism_sql_obj = Organisms_sql_new()


list_bacts = bact_obj.get_all_Bacteria()


print(len(list_bacts))
print(list_bacts[0].GI)


print("qweqweqweqwe")
list_bacts[0].complete_bacteria_from_old_DB()


print("111111111")


## %% final production list_phages
#from objects_old.Phages_old import *
#
#phage_old_obj = Phage_old()
#
#list_phages = phage_old_obj.get_all_phages()
#
#
#list_phages[0].complete_phage_from_old_DB()
#
#print(list_phages[0].phage_id)
#print(list_phages[0].name)
#print(list_phages[0].GI)
#print(list_phages[0].nb_proteins)
#print(list_phages[0].whole_genome)
#print(list_phages[0].dna_code_sequence)
#print(list_phages[0].prot_sequence)
#
#aux = 0
#
#for phage in list_phages:
#    phage.complete_phage_from_old_DB()
#    print('--- Start name instead of strain---')
#    strain = Strain()
#    id_strain_inserted = strain.create_strain(designation=phage.name, fk_specie= 1)
#    print("---End name instead of strain---")    
#    
#    
#    print("--Start Whole Genome---")
#    print(phage.phage_id)
#    dna_whole_genome = WholeDNA()
#    id_whole_genome_inserted = dna_whole_genome.create_whole_dna(phage.whole_genom_obj.head, phage.whole_genom_obj.head_id, phage.whole_genom_obj.sequence)
#    print(id_whole_genome_inserted)
#    print("--End Whole Genome---")    
#
#    print("---Start Organism---")
#    organism = Organism(-1, phage.GI, 'NA', -1, 0, -1, 1, 1, 2, id_whole_genome_inserted)
#    id_organism_inserted = organism.create_organism()
#    print(id_organism_inserted)
#    print("---End Organism---")
#
#
#    print("---Start Proteins and Genes---")
#
#    print("Il y a: " + str(len(phage.proteins_list)))
#    aux = 0
#    for protein_values in phage.proteins_list:
#        value_id_protein = protein_values.create_protein()
#        gene_obj = Gene(FK_id_organism = id_organism_inserted, FK_id_protein = value_id_protein)
#        value_id_gene = gene_obj.create_gene()
#        aux += 1
#        if aux % 1500 == 0:
#            print("sleep")
#            time.sleep(3)
#        #print(value_id_protein)
#        #print(value_id_gene)
#    print("Il fut insere: " + str(aux))
#    print("il y avait nb prots nucleo counted: " + str(phage.number_prots_counted_nucleo))
#    print("il y avait nb prots aa counted: " + str(phage.number_prots_counted_aa))
#    time.sleep(1)


# %% final production list_bacts
for bacterium in list_bacts:
    bacterium.complete_bacteria_from_old_DB()
    print("-----------------")
    print(bacterium.GI)
    print("verify if organisme exists")
    organism_exists = organism_sql_obj.get_organisme_by_GI(bacterium.GI)
    if organism_exists == True:
        print("Organisme exists so we pass(")
        continue
    
    bacteria_file = BacteriaFile(bacterium)
    bacteria_file.get_obj_id()
    print(bacteria_file.bacterium)
    print(bacteria_file.file_sequece_info)
    print(bacteria_file.file_taxonomy)
    print(bacteria_file.objectId)
    print(bacteria_file.genus)
    print(bacteria_file.specie)
    print(bacteria_file.family)
    print(bacteria_file.strain)
    print(bacteria_file.accession_num)
    print("-----------------")
    print("---Family---")
    fam = Family()
    id_fam_inserted = fam.create_family(designation = bacteria_file.family)
    print(id_fam_inserted)
    print("---End family---")
    print("---Start Genus---")
    gen = Genus()
    id_genus_inserted = gen.create_genus(designation = bacteria_file.genus, fk_family = id_fam_inserted)
    print("ID: " + str(id_genus_inserted))
    print("---End Genus---")
    print("---Start Specie---")
    specie = Specie()
    id_specie_inserted = specie.create_specie(designation=bacteria_file.specie, fk_genus= id_genus_inserted)
    print("---End Specie---")
    print("---Start Strain---")
    strain = Strain()
    id_strain_inserted = strain.create_strain(designation=bacteria_file.strain, fk_specie= id_specie_inserted)
    print("---End Strain---")
   
    print("--Start Whole Genome---")
    print(bacterium.GI)
    dna_whole_genome = WholeDNA()
    id_whole_genome_inserted = dna_whole_genome.create_whole_dna(bacterium.whole_genom_obj.head, bacterium.whole_genom_obj.head_id, bacterium.whole_genom_obj.sequence)
    print(id_whole_genome_inserted)
    print("--End Whole Genome---")
    print("---Start Organism---")
    organism = Organism(-1, bacterium.GI, bacteria_file.accession_num, -1, 0, -1, 1, id_strain_inserted, 1, id_whole_genome_inserted)
    id_organism_inserted = organism.create_organism()
    print(id_organism_inserted)
    print("---End Organism---")
    print("---Start Proteins and Genes---")

    print("Il y a: " + str(len(bacterium.proteins_list)))
    aux = 0
    for protein_values in bacterium.proteins_list:
        value_id_protein = protein_values.create_protein()
        gene_obj = Gene(FK_id_organism = id_organism_inserted, FK_id_protein = value_id_protein)
        value_id_gene = gene_obj.create_gene()
        aux += 1
        if aux % 1500 == 0:
            print("sleep")
            time.sleep(3)
        #print(value_id_protein)
        #print(value_id_gene)
    print("Il fut insere: " + str(aux))
    print("il y avait nb prots nucleo counted: " + str(bacterium.number_prots_counted_nucleo))
    print("il y avait nb prots aa counted: " + str(bacterium.number_prots_counted_aa))
    
    
    time.sleep(3)
    
    
    
    print("---End Proteins and Genes---")
