from files_treatment_new.fasta_nucleotid_NCBI import Fasta_nucleotid_NCBI
from files_treatment_new.fasta_protein_NCBI import Fasta_protein_NCBI
from files_treatment_new.genbank_file_NCBI import Genbank_proteic_NCBI

from objects_new.Proteins_new import Protein

#taxonomy
from objects_new.Families_new import *
from objects_new.Genus_new import *
from objects_new.Species_new import *
from objects_new.Strains_new import *
from objects_new.Organisms_new import *
from objects_new.Gene_new import *
from objects_new.Contigs_new import *
from objects_new.WholeDNA_new import WholeDNA

import os
from os import listdir
from os.path import isfile, join




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

cwd = os.getcwd()

list_files = get_list_ids_files_in_path('/Phages_NCBI/PROTEINS/')


for file_name in list_files:
    path_file_protein = cwd + '/Phages_NCBI/PROTEINS/' + file_name
    path_file_nucleotid = cwd + '/Phages_NCBI/NUCLEOTIDS/' + file_name[:-4] + '.fna'
    path_file_genbank = cwd + '/Phages_NCBI/WHOLE_GENOME_GENBANK/Whole_genome_' + file_name[:-4] + '.gb'

    ncbi_file_data_prot = Fasta_protein_NCBI(path_file_protein)
    ncbi_file_data_nucleotid = Fasta_nucleotid_NCBI(path_file_nucleotid)
    ncbi_file_genbank = Genbank_proteic_NCBI(path_file_genbank)

    accession_phage = ncbi_file_genbank.get_accession_number()
    contig_sequence_phage = ncbi_file_genbank.get_contig_of_the_phage()

    list_ids_prot = ncbi_file_data_prot.get_list_protein_name()

    list_of_proteins = []
    for id_prot in list_ids_prot:
        sequence_prot = ncbi_file_data_prot.get_protein_sequence_prot_id(id_prot)
        sequence_nuc = ncbi_file_data_nucleotid.get_nucleotid_sequence_prot_id(id_prot)
        vec_location = ncbi_file_data_nucleotid.get_location_by_prot_id(id_prot)

        protein_obj = Protein(id_accession = id_prot, sequence_prot=sequence_prot, sequence_dna= sequence_nuc, start_point=vec_location[0], end_point=vec_location[1] )
        list_of_proteins.append(protein_obj)

    family_obj = None
    genus_obj = None
    specie_obj = None
    strain_obj = None

    family_obj = Family(designation = 'Phage no family')
    genus_obj = Genus(designation = 'Phage no genuse')
    specie_obj = Specie(designation = 'Phage no Specie')
    strain_obj = Strain(designation = file_name[:-4])
    print(strain_obj)

    family_obj.create_family()
    genus_obj.fk_family = family_obj.id_family
    genus_obj.create_genus()
    specie_obj.fk_genus = genus_obj.id_genus
    specie_obj.create_specie()
    strain_obj.fk_specie = specie_obj.id_specie
    strain_obj.create_strain()
    id_strain = strain_obj.id_strain

    qty_prots_loaded = len(list_of_proteins)

    #Whole genome
    whole_dna_obj = WholeDNA(head = "Unknown " + file_name[:-4], head_id = "Unknown", sequence = "Unknown")
    whole_dna_obj.create_whole_dna_no_verification()
    id_whole_genom = whole_dna_obj.id_wholeDNA

    #Contig
    contig_obj = Contig(head = 'no head ' + file_name[:-4], sequence = contig_sequence_phage)
    contig_obj.fk_id_whole_genome = id_whole_genom
    contig_obj.create_contig_no_verification()
    id_contig = contig_obj.id_contig

    organism_obj = Organism(id_organism = -1, gi = 'no GI', acc_num = accession_phage, qty_proteins = qty_prots_loaded, assembled = True, fk_source = 2, fk_strain = id_strain, fk_type = 2, fk_whole_genome = id_whole_genom, fk_source_data = 1, qty_contig = 1)

    organism_obj.create_organism()

    print(len(list_of_proteins))

    for protein in list_of_proteins:
        protein.fk_id_contig = contig_obj.id_contig
        id_protein = protein.create_protein()
        print(protein.start_point_cnt)
        print(protein.end_point_cnt)
        print(protein.fk_id_contig)

        #id_protein = protein.create_protein()

        gene_obj = Gene(gene_number = -1, dna_head = "No head", 
                        dna_sequence = "No sequence", start_position = -1, 
                        end_position = -1, FK_id_organism = organism_obj.id_organism, FK_id_protein = id_protein)

        gene_obj.create_gene()
