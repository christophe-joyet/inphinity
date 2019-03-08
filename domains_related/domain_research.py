# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 14:30:47 2017

@author: Diogo
"""


#http://pydoc.net/bcbio-nextgen/1.0.4/bcbio.hmmer.search/

import time

from objects_new.Proteins_new import *
from objects_new.Protein_dom_new import *
from objects_new.Domains_new import *
from objects_new.temp_prot import *

from objects_new.COGS_new import *

from pfam_detection import *

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed

import multiprocessing


temp_prot = id_prot_temp()
temp_prot.id_prot_temp = -4567

temp_prot.create_id_prot()

#list_temp_ids_prot = temp_prot.get_all_ids_prot()

list_temp_ids_prot = [0]
print(list_temp_ids_prot)


#list_ids_organisms = [4648,4651,4653,4663]

list_ids_organisms = [4663,4674,4683,4689,4696,4742,4793,4828,4848,4946,5021,4677,5027,5034,5043,5054,5067,5073,5110,5124,5142,5150,5170,5205,5225,5294,5314,5346,5467,5468,5518,5523,5580,5601,5614,5621,5622,5640,5652,5660,5670,5675,5718,5727,5861,5866,5883,5899,5960,5968,5983,5984,6013,6036,6091,6114,6155,6156,6186,6199,6215,6228,6295,6334,6336,6341,5904,6414,6425,6431,6472,6528,6537,6574,6593,6755,8191,8219,8256,8959,8958,8960,8957,8964,8965,9360,9359,9361,8622,8812,10093,8773,4791,9278,8666,5025,8001,8783,9638,9490,7814,7095,6219,9538,9942,8785,9657,6354,5811,8361,10079,6060,8787,8833,9748,9630,10143,9739,10113,9994,7765,10130,10036,10142,8874,8784,8744,10136,10131,10567,7832,9659,8626,10121,8085,8664,8774,9891,4718,8782,4752,7613,7720,8284,8802,8293]

list_ids_organisms = [8802, 8293]

prot_obj = Protein()
#print("Start get all proteins")
#listOfProteins = Protein.get_all_Proteins_by_organism_id(10504)
#print("End get all proteins")
#print(len(listOfProteins))
listOfProteins = []
print("Start get all proteins")
for id_organism in list_ids_organisms:
    print(id_organism)
    listOfProteins_aux = Protein.get_all_Proteins_by_organism_id(id_organism)
    listOfProteins = listOfProteins + listOfProteins_aux
print(len(listOfProteins))
print("End get all proteins")
#print(listOfProteins[5].sequence_prot)
#print(listOfProteins[5].designation)
#print(listOfProteins[5].id_accession)

#print(listOfProteins[5].id_protein)


#pfamDetectObj = PFAM_detection(listOfProteins[10].sequence_prot)

#list_dom = pfamDetectObj.detecterPFAM()

#print(list_dom)


domain_obj = Domain()
list_dom = Domain.get_all_Domains()
#domain_obj.designation = list_dom[0]
#id_domain = domain_obj.create_domain()

#prot_dom_obj = ProteinDom()
#prot_dom_obj.Fk_id_protein = listOfProteins[10].id_protein
#prot_dom_obj.Fk_id_domain = id_domain
#prot_dom_obj.create_protDom_if_not_exist()

print("----------")

#********************** thread manipulation
# %% Search for all proteins (if find same domains they arn't add)
#def serach_domains(prot_obj):
#    print("Protein id: " + str(prot_obj.id_protein))
#    
#    pfamDetectObj = PFAM_detection(prot_obj.sequence_prot)
#    prot_dom_obj = ProteinDom()
#    
#    qty_protein = prot_dom_obj.protein_verified(prot_obj.id_protein)
#    if qty_protein > 0:
#        print("Tested")
#    else:
#        print("Not Tested-----------------------------------")
#        list_dom = pfamDetectObj.detecterPFAM()
#        id_prot = prot_obj.id_protein
#        for domain in list_dom:
#            id_domain = domain_obj.create_domain(domain)
#            prot_dom_obj.create_protDom_if_not_exist(id_prot, id_domain)
#        print("domains created: " + str(len(list_dom)))
        
# %% Search for new proteins
def serach_domains(prot_obj):
    print("Protein id: " + str(prot_obj.id_protein))
    if prot_obj.id_protein not in list_temp_ids_prot:
        pfamDetectObj = PFAM_detection(prot_obj.sequence_prot, prot_obj.id_protein)
        prot_dom_obj = ProteinDom()
        prot_dom_obj.Fk_id_protein = prot_obj.id_protein
        print(prot_obj.id_protein)
        qty_protein = prot_dom_obj.protein_number_domains()
        if qty_protein > 0:
            print("Tested")
        else:
            pfamDetectObj.id_prot = prot_obj.id_protein
            list_dom = pfamDetectObj.detecterPFAM()
            id_prot = prot_obj.id_protein
            for domain in list_dom:
                
                domain = domain.split('.')
                

                domain_obj.designation = domain[0]


                id_domain = domain_obj.create_domain()

                prot_dom_obj.Fk_id_domain = id_domain
                prot_dom_obj.Fk_id_protein = id_prot
                prot_dom_obj.create_protDom_if_not_exist()
            print("domains created: " + str(len(list_dom)))
    else:
        print("protein %s already tested".format(str(prot_obj.id_protein)))
    
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

print("FIni")
        
        
        
        


















