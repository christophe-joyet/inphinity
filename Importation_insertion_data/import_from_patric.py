# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:01:54 2018

@author: Diogo
"""

class importation_patric(object):
    """
    This class import the data from the files obtaines in PATRIC platform

    """

    def __init__(self, xls_file, sheet_xls_name, contig_fasta_file):

        """
        Constructor of the importation from patric class object. All the parameters have a default value

        :param xls_file: name of the xls file (that contain information about the proteins in the contig)
        :param sheet_xls_name: name of the sheet that contain the data
        :param contig_fasta_file:name of the contig fasta file (that contain the contig sequences)

        :type xls_file: string - required
        :type sheet_xls_name: string - required
        :type contig_fasta_file: string - required 
        """
        self.xls_file = xls_file
        self.sheet_xls_name = sheet_xls_name
        self.contig_fasta_file = contig_fasta_file

        self.xls_data = None
        self.contig_fasta = None
        self.qty_proteins = 0
        self.qty_contig = 0

    def load_data(self):
        """
        Load data from the exls and fasta files (obtained from Patric)
        """
        try:
            self.xls_object = Xls_gen_bank(path_file = self.xls_file.encode('string-escape'), sheet_name = self.sheet_xls_name.encode('string-escape'))
            self.contig_file = Fasta_contigs_RAST(path_file = self.contig_fasta_file.encode('string-escape'))
        except FileNotFoundError:
            print("Files not found...")

        data_validation_rslt = self.data_validation()
        assert data_validation == True

    def data_validation(self):
        """
        This method is used to validate the correct proteins importation according the Contigs. True is return in case of all the data are correct and False in the contrary 

        :return: True or False according the data validation
        :rtype boolean
        """
        data_validated = False
        qty_prots_contigs = 0
        qty_proteins_organism = self.xls_object.get_number_of_proteins()
        list_contigs_obj = self.contig_file.create_contigs_from_file()


        self.qty_contig = len(list_contigs_obj)

        for contig in list_contigs_obj:
            list_proteins = xls_object.get_proteins_objects_by_contig_id(contig.head)
            qty_prots_contigs += len(list_proteins)

        self.qty_proteins = qty_proteins_organism

        if qty_proteins_organism == qty_prots_contigs and qty_proteins_organism > 0 and qty_contig > 0:
            data_validated = True
        else:
            data_validated = False
        return data_validated

    def insert_data_db(self, strain_name, fk_specie, type_of_organism):
        #verifié si strain exist
        strain_obj = Strain(designation = designation_strain, fk_specie = fk_specie)

        id_strain = strain_obj.create_strain()

        #Creation of the Whole DNA
        whole_dna_obj = WholeDNA(id_wholeDNA = -1, head = "Unknown_test", head_id = "Unknown", sequence = "Unknown")
        id_whole_dna = whole_dna_obj.create_whole_dna_no_verification()

        #Creation of the Organism
        organism_obj = Organism(gi = "", acc_num = "", qty_proteins = self.qty_proteins, assembled = True, qty_contig = self.qty_contig, fk_source = 2, fk_strain = id_strain, fk_type = type_of_organism, fk_whole_genome = id_whole_dna, fk_source_data = 5)
        id_organism_Created = organism_obj.create_organism()

        #Creation of the proteins
        for contig in list_contigs_obj:
            list_proteins = xls_object.get_proteins_objects_by_contig_id(contig.head)
            contig.fk_id_whole_genome = id_whole_dna
            contig_id = contig.create_contig_no_verification()
            for protein_obj in list_proteins:
                protein_obj.fk_id_contig = contig_id
                id_protein = protein_obj.create_protein()
                gene_obj = Gene(FK_id_organism = id_organism_Created, FK_id_protein = id_protein)
                gene_obj.create_gene()
