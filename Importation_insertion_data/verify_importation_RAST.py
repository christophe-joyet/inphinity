# -*- coding: utf-8 -*-
"""
Created on Wen Jan 31 14:44:17 2018

@author: Diogo
"""

import csv
import os

from files_treatment_new.fasta_protein_RAST import Fasta_protein_RAST
from files_treatment_new.fasta_nucleotid_RAST import Fasta_nucleotid_RAST

from objects_new.Proteins_new import *

list_files = []

with open('success_files_list.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        list_files.append(row[0][:-4])

print(list_files)

cwd = os.getcwd()


file_name = list_files[0]

#confirm proteins:
path_file_proteic_fasta = cwd + '/RAST/PROTEIC_FAA/' + file_name + '.faa'
protein_fasta_obj = Fasta_protein_RAST(path_file_proteic_fasta)
list_prots = protein_fasta_obj.create_list_proteins()
#complete proteins with DNA sequence
path_file_proteic_fasta = cwd + '/RAST/NUCLEIC_FNA/' + file_name + '.fna'
protein_fasta_nuc_obj = Fasta_nucleotid_RAST(path_file_proteic_fasta)
list_prots = protein_fasta_nuc_obj.complete_list_proteins(list_prots)

id_prot = list_prots[0].get_id_by_acc_seqProt_seqNuc()

a = 1
