from files_treatment_new.xls_gen_bank_rast import Xls_gen_bank
from files_treatment_new.fasta_contigs_RAST import Fasta_contigs_RAST

from files_treatment_new.genbank_file_RAST import Genbank_proteic_RAST

from files_treatment_new.fasta_contigs_RAST import Fasta_contigs_RAST

from objects_new.WholeDNA_new import WholeDNA

from objects_new.Proteins_new import *

import csv

#taxonomy
from objects_new.Families_new import *
from objects_new.Genus_new import *
from objects_new.Species_new import *
from objects_new.Strains_new import *
from objects_new.Organisms_new import *
from objects_new.Gene_new import *
from objects_new.Proteins_new import *


import os
from os import listdir
from os.path import isfile, join

import time


def get_list_ids_files_in_path(path):
    """
    This method list all files in a given path and return a list with these names

    :param path: path where it necessary to list the files

    :type path: string - required

    :return list with the files paths
    :rtype list(str)

    :note when the start point is smaller than end point (int the contig), it is because the "Strand field int excel file is negative
    """

    current_path = os.getcwd() + path
    list_files = os.listdir(current_path)

    return list_files

def check_file_exits(file_path):
    """
    This method just verify if a given file exists (it is necessary to give the complete path)

    :param file_path: complete path of the file

    :type file_path: string - required

    :return True or False according the existance
    :rtype boolean

    """
    file_exists = os.path.exists(file_path)
    if file_exists is True:
        return True
    else:
        return False


#Def write list into file
def write_list_into_file(list, file_name):
    csvfile = file_name
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in list:
            writer.writerow([val])  


def get_files_from_error_csv(path_csv):
    list_files = None
    list_element = []
    with open(path_csv, 'r') as f:
        reader = csv.reader(f)
        list_files = list(reader)
    for element in list_files:
        list_element.append(element[0])
    return list_element
#organism_obj_a = Organism(id_organism = -1, gi = 'qweqwe', acc_num = 'NC_004679', qty_proteins = 4556, assembled = True, qty_contig = 456, fk_source = -1, fk_strain = -1, fk_type = 2)



#print("--------------------")
#asdasd = organism_obj_a.get_id_organism_by_acc()
#print(asdasd)

#eeeee = Organism.get_id_organism_by_designation('Mycobacterium phage Nyxis')
#print(eeeee)

#qqqq = Organism.get_id_organism_by_acc('NC_004679')
#print(qqqq)


#qweqweqwe = Organism.get_organism_id_by_acc_or_designation('Mycobacterium phage Nyxis', 'na')

#print(eeeee)
#print("--------------------")

cwd = os.getcwd()



list_of_proteins = []

#this line takes the files from csv_error_files
list_files = get_files_from_error_csv('error_files_list.csv')


#Uncomment this line for the first insertion
#list_files = get_list_ids_files_in_path('/RAST/SPREADSHEET/')
list_files = reversed(list_files)
list_files_error = []
list_files_done = []
for file_name in list_files:
    print(file_name)
    #Load taxonomy
    #try:
    path_file_genBank = cwd + '/RAST/GEN_BANK/' + file_name[:-3] + 'gbk'
    gen_bank_obj = Genbank_proteic_RAST(path_file_genBank)
    taxo = gen_bank_obj.get_taxonomy_array()
    print(gen_bank_obj.get_family())
    if len(taxo) == 7 or 'bacteriophage' in file_name.lower():

        family_obj = None
        genus_obj = None
        specie_obj = None
        strain_obj = None
        if ('_phi' not in file_name.lower() or 'phage' not in file_name.lower()) and len(taxo) == 7:
            family_obj = Family(designation = gen_bank_obj.get_family())
            genus_obj = Genus(designation = gen_bank_obj.get_genus())
            specie_obj = Specie(designation = gen_bank_obj.get_specie())
            strain_obj = Strain(designation = gen_bank_obj.get_strain())
        else:
            name_phage = file_name.split('-')[1]
            name_phage = name_phage.split('.')[0]
            name_phage = name_phage.replace('_', ' ')
            print(name_phage)
            family_obj = Family(designation = 'Phage no family')
            genus_obj = Genus(designation = 'Phage no genuse')
            specie_obj = Specie(designation = 'Phage no Specie')
            strain_obj = Strain(designation = name_phage)

        family_obj.create_family()
        genus_obj.fk_family = family_obj.id_family
        genus_obj.create_genus()
        specie_obj.fk_genus = genus_obj.id_genus
        specie_obj.create_specie()
        strain_obj.fk_specie = specie_obj.id_specie
        strain_obj.create_strain()





        #Test proteins from file
        path_file_xls = cwd + '/RAST/SPREADSHEET/' + file_name
        path_file_contig = cwd + '/RAST/CONTIGS/' + file_name[:-3] + 'contigs.fa'

        xls_obj = Xls_gen_bank(path_file_xls)

        value = check_file_exits(path_file_xls)

        contig_file_exist = check_file_exits(path_file_contig)

        #Test contigs from file
        if contig_file_exist is True:
            fasta_contig_obj = Fasta_contigs_RAST(path_file_contig)
            qtity_cntg = fasta_contig_obj.get_qty_of_contigs()
            list_cnt = fasta_contig_obj.create_contigs_from_file()
            print("Hello")

        list_of_proteins = xls_obj.create_proteins_from_file()
        qty_proteins_loaded = len(list_of_proteins)
        qty_contigs_loaded = xls_obj.get_number_different_contigs()
        print(qty_contigs_loaded)
        print(qty_proteins_loaded)
   
        #Test empty WHole Genome
        #, id_wholeDNA = -1, head = "", head_id = "", sequence = ""
        whole_dna_obj = WholeDNA(head = "Unknown", head_id = "Unknown", sequence = "Unknown")
        whole_dna_obj.create_whole_dna_no_verification()


        #Create an organism:
        #id_organism, gi, acc_num, qty_proteins, assembled, qty_contig, fk_source = -1, fk_strain = -1, fk_type = -1, fk_whole_genome = -1, fk_source_data = "NULL"
        gi_name = "Greg_" + strain_obj.designation
        acc_value = "Greg_" + strain_obj.designation
        id_strain = strain_obj.id_strain

        print("insert the Whole DNA")
        id_whole_genom = whole_dna_obj.id_wholeDNA

        fk_type_organism = -1
        if "_phi" not in file_name or "phage" not in file_name:
            fk_type_organism = 1
        else:
            fk_type_organism = 2
        organism_obj = Organism(id_organism = -1, gi = gi_name, acc_num = acc_value, qty_proteins = qty_proteins_loaded, assembled = True, fk_source = 2, fk_strain = id_strain, fk_type = fk_type_organism, fk_whole_genome = id_whole_genom, fk_source_data = 3, qty_contig = qty_contigs_loaded)
        print(organism_obj)

        organism_obj.create_organism()




        print("Start the insertion")
        print("Insert the organism")
        #id_organism_inserted = organism_obj.create_organism()

        print("Insert contigs")
        for contig in list_cnt:
            contig.fk_id_whole_genome = id_whole_genom
            contig.create_contig_no_verification()
            print(contig)
            list_of_proteins = xls_obj.get_proteins_objects_by_contig_id(contig.head)
            print(len(list_of_proteins))
            for protein in list_of_proteins:
                protein.fk_id_contig = contig.id_contig
                id_protein = protein.create_protein()
                print(protein.start_point_cnt)
                print(protein.end_point_cnt)
                print(protein.fk_id_contig)

                print(type(protein.start_point_cnt))
                print(type(protein.end_point_cnt))
                print(type(protein.fk_id_contig))
                protein.update_protein_contig()

                gene_obj = Gene(gene_number = -1, dna_head = "No head", 
                 dna_sequence = "No sequence", start_position = -1, 
                 end_position = -1, FK_id_organism = organism_obj.id_organism, FK_id_protein = id_protein)

                gene_obj.create_gene()


        list_files_done.append(file_name)
        write_list_into_file(list_files_done, "success_files_list.csv")

    else:
        list_files_error.append(file_name)
    #except ValueError:
    #    list_files_error.append(file_name)
        write_list_into_file(list_files_error, "error_files_list.csv")


print(list_files_error)

