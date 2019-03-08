# -*- coding: utf-8 -*-
"""
Created on Wed Mac 18 14:30:47 2017

@author: Diogo
"""

#/HMMER/hmmer-3.1b2/db$

from objects_new.Proteins_new import *
from objects_new.Protein_dom_new import *
from objects_new.Domains_new import *
from objects_new.temp_prot import *
import tempfile
import os
import re


from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed

import multiprocessing

#Regular expretion to match the domains:
expression = r"PF[0-9]{5}"


temp_prot = id_prot_temp()
list_temp_ids_prot = temp_prot.get_all_ids_prot()





prot_obj = Protein()
print("Start get all proteins")
listOfProteins = prot_obj.get_all_Proteins_limited(1000000, 999999)
print("End get all proteins")
print(len(listOfProteins))




#Search all proteins:
aux_index = 0
def serach_domains(protein):
    global aux_index
    if protein.id_protein not in list_temp_ids_prot:
        print("Protein number %d with id: %d" % (aux_index, protein.id_protein))
        #Creation of temporary files
        fp = tempfile.NamedTemporaryFile()
        fp_out = tempfile.NamedTemporaryFile()
        fp_out_content = tempfile.NamedTemporaryFile()

        fp.write(('>' + protein.designation + '\n').encode())
        fp.write(str(protein.sequence_prot).encode())

        fp.seek(0)
        aux = fp.read()
        temp_name_fasta = fp.name
        temp_name_out = fp_out.name
        temp_name_output = fp_out_content.name

        #Ubuntu command used for the domaine discovery
        temp_command = r"hmmscan --tblout " + temp_name_out + r" --acc --noali -E 1e-10 --qformat fasta /home/diogo.leite/HMMER/hmmer-3.1b2/db/Pfam-A.hmm " + temp_name_fasta

        results_cmd = os.system(temp_command + " > " + temp_name_output)

        temporary_content = fp_out.read().decode()


        #Find domains
        matches = re.findall(expression, temporary_content)
        prot_dom_obj = ProteinDom()
        domain_obj = Domain()

        if len(matches) > 0:
            prot_dom_obj.Fk_id_protein = protein.id_protein
            for domain in matches:
                print(domain)
                domain_obj.designation = domain
                id_domain = domain_obj.create_domain()
                prot_dom_obj.Fk_id_domain = id_domain
                prot_dom_obj.Fk_id_protein = protein.id_protein
                prot_dom_obj.create_protDom_if_not_exist()


        fp.close()
        fp_out.close()


        #insert id_protein
        temp_prot = id_prot_temp()
        temp_prot.id_prot_temp = protein.id_protein

        temp_prot.create_id_prot()
        aux_index += 1
    else:
        print("Protein number %d with id: %d already done" % (aux_index, protein.id_protein))



# %% Start multiprocess pool
max_proce = multiprocessing.cpu_count()
print("------------------------------")
print("you have these processores: " + str(max_proce))
print("------------------------------")
max_proce = max_proce - 2


with ThreadPoolExecutor(max_workers=max_proce) as executor:
    futures = [executor.submit(serach_domains, protein_obj) for protein_obj in listOfProteins]
    for future in as_completed(futures):
        print(future.result())
