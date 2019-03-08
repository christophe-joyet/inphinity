# -*- coding: utf-8 -*-
"""
Created on Fri May 4 08:27:54 2018

@author: Diogo
"""

from objects_new.Couples_new import *
from objects_new.Proteins_new import *
from objects_new.Protein_dom_new import *
from objects_new.DDI_interactions_DB_new import *
from objects_new.PPI_preview_new import *
import time



listasd = DDI_interaction_DB.get_DDI_interact_db_id_by_id_domains(10165, 10165)

#qty_val = PPI_preview.get_number_ppi_score_by_bact_phage_prots(1977507, 27085575)



def calcule_score_DDI(list_scores_ddi):
    #calcule the DDI score iPFAM, 3DID and ME have the maximum score ->9
    score_ddi = 0
    type_score = 0
    if 1 in list_scores_ddi or 2 in list_scores_ddi:
        score_ddi = 9
        type_score = 1 
    elif 3 in list_scores_ddi:
        score_ddi = 9
        type_score = 2
    else:
        score_ddi = len(list_scores_ddi)
        type_score = 3
    return score_ddi, type_score


def get_couples_NCBI_phagesDB():
    #Get all couples from NCBI and Phages DB
    list_couples = Couple.get_all_couples_by_type_level_source(1, 4, 1)
    list_couples = list_couples + Couple.get_all_couples_by_type_level_source(1, 4, 2)
    list_couples = list_couples + Couple.get_all_couples_by_type_level_source(0, 4, 4)
    return list_couples

def get_couples_Greg():
    #Get all couples from Greg
    list_couples = Couple.get_all_couples_by_type_level_source(1, 4, 3)
    list_couples = list_couples + Couple.get_all_couples_by_type_level_source(0, 4, 3)
    return list_couples


list_couples = get_couples_NCBI_phagesDB()
print(len(list_couples))
print("qwqwqwqwqwqw")





#id_couple = list_couples[0].id_couple
auxaa = 0
for couple in list_couples:
    if couple.id_couple > 34943:

        print(couple)
        #for each couple
        list_prot_bact = Protein.get_all_Proteins_with_domains_by_organism_id(couple.fk_bacteria)
        list_prot_phage = Protein.get_all_Proteins_with_domains_by_organism_id(couple.fk_phage)
        print("It has {0:d} proteins for the bacterium".format(len(list_prot_bact)))
        print("It has {0:d} proteins for the Phages".format(len(list_prot_phage)))

        qty_ppi = PPI_preview.get_ppi_preview_scores_grouped_by_couple_id(couple.id_couple)

        aux_value = 0
        print("Quantity of PPI {0}".format(len(qty_ppi)))
        if len(qty_ppi) == 0:
            #for all prots in bact combine with those of the phage
            for prot_bact in list_prot_bact:
                #get all domains by protein id
                list_domain_id_prot_Bact = ProteinDom.get_all_protein_domain_by_protein_id(prot_bact.id_protein)
                #for all proteins in phage
                for prot_phage in list_prot_phage:
                    #get all domains by protein id
                    list_domain_id_prot_phage = ProteinDom.get_all_protein_domain_by_protein_id(prot_phage.id_protein)
                    #combiner tous les DDI et faire les insérer dans le tableau PPI preview
                    for dom_bact in list_domain_id_prot_Bact:
                        for dom_phage in list_domain_id_prot_phage:
                            liste_scores_ddi = DDI_interaction_DB.get_DDI_interact_db_id_by_id_domains(dom_bact.Fk_id_domain, dom_phage.Fk_id_domain)
                            if len(liste_scores_ddi) > 0:
                                #print(couple)
                                #print(aux_value)
                                score_ddi, type_score = calcule_score_DDI(liste_scores_ddi)
                                ppi_prev_obj = PPI_preview(score_ppi_prev = score_ddi, type_ppi_prev = type_score, fk_couple = couple.id_couple, fk_prot_bact = prot_bact.id_protein, fk_prot_phage = prot_phage.id_protein)
                                ppi_prev_obj.create_ppi_preview()
                                auxaa += 1
                                if auxaa %100 == 0:
                                    time.sleep(5)
                                print(ppi_prev_obj)
                    aux_value += 1
        time.sleep(3)

print(aux_value)