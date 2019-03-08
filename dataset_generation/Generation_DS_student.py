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
from objects_new.Ddi_interactions_new import *
import numpy as np


def get_domains_by_list_prots(list_proteins):
    dict_prot_domains = {}
    list_domains_id = []
    for protein in list_proteins:
        list_domain_id_prot_phage = ProteinDom.get_all_protein_domain_by_protein_id(protein.id_protein)
        list_domains_id = []
        for domain in list_domain_id_prot_phage:
            list_domains_id.append(domain.Fk_id_domain)
        if len(list_domains_id ) > 0:
            dict_prot_domains[protein.id_protein] = list_domains_id
    return dict_prot_domains

def get_interact_ddi_id(key_dom_bact, key_dom_phage):
    global dict_ddi_interact
    id_interaction = -1
    pair_key = (key_dom_bact, key_dom_phage)
    pair_key_invert = (key_dom_phage, key_dom_bact)
    if pair_key in dict_ddi_interact:
        id_interaction = dict_ddi_interact[pair_key]
    elif pair_key_invert in dict_ddi_interact:
        id_interaction = dict_ddi_interact[pair_key_invert]
    return id_interaction

def get_dict_key_interact():
    dict_DDI = {}
    list_ddi_interact = DDI_interaction.get_all_DDI_interaction()
    for interaction in list_ddi_interact:
        key = (interaction.FK_id_domain_A, interaction.FK_id_domain_B)
        value = interaction.id_ddi_interaction
        dict_DDI[key] = value
    return dict_DDI

def get_DDI_by_DB_in_dic():
    dict_id_inter_db = {}
    list_ddi_inter_db = DDI_interaction_DB.get_all_DDI_interactionDB()
    for ddin_ter_db in list_ddi_inter_db:
        if ddin_ter_db.FK_DDI_interaction in dict_id_inter_db and ddin_ter_db.FK_DB_source != 16:
            dict_id_inter_db[ddin_ter_db.FK_DDI_interaction].append(ddin_ter_db.FK_DB_source)
        elif ddin_ter_db.FK_DB_source != 16:
            dict_id_inter_db[ddin_ter_db.FK_DDI_interaction] = [ddin_ter_db.FK_DB_source]
    return dict_id_inter_db



def get_couples_NCBI_phagesDB():
    #Get all couples from NCBI and Phages DB
    list_couples = Couple.get_all_couples_by_type_level_source(1, 4, 1)
    list_couples = list_couples + Couple.get_all_couples_by_type_level_source(1, 4, 2)
    list_couples = list_couples + Couple.get_all_couples_by_type_level_source(0, 4, 4)
    size_couples = len(list_couples)
    new_size = len(list(set(list_couples)))
    assert size_couples == new_size


    return list_couples


def get_couples_Greg():
    #Get all couples from NCBI and Phages DB
    list_couples = Couple.get_all_couples_by_type_level_source(1, 4, 3)
    list_couples = list_couples + Couple.get_all_couples_by_type_level_source(0, 4, 3)
    return list_couples

def calculat_in_rules_ddi(id_interaction):
    global dict_interactions_DB
    score_ddi = 0
    type_score = 0
    vec_values_ddi_scores = dict_interactions_DB[id_interaction]
    if 1 in vec_values_ddi_scores or 2 in vec_values_ddi_scores:
        score_ddi = 9
        type_score = 1
    elif 3 in vec_values_ddi_scores:
        score_ddi = 9
        type_score = 2
    else:
        score_ddi = len(vec_values_ddi_scores)
        type_score = 3
    return score_ddi, type_score


def search_id_interactions(vec_dom_prot_bac, vec_dom_prot_phage, id_couple, id_prot_bact, id_prot_phage):
    global dict_ddi_interact
    list_ppi_preview = []
    interaction_ddi_id = -1
    for dom_bact in vec_dom_prot_bac:
        for dom_phage in vec_dom_prot_phage:
            if (dom_bact, dom_phage) in dict_ddi_interact:
                interaction_ddi_id = dict_ddi_interact[(dom_bact, dom_phage)]
                if interaction_ddi_id not in dict_interactions_DB:
                    interaction_ddi_id = -1
            elif (dom_phage, dom_bact) in dict_ddi_interact:
                interaction_ddi_id = dict_ddi_interact[(dom_phage, dom_bact)]
                if interaction_ddi_id not in dict_interactions_DB:
                    interaction_ddi_id = -1
            if interaction_ddi_id > -1:
                score_ddi, type_score = calculat_in_rules_ddi(interaction_ddi_id)
                interaction_ddi_id = -1
                ppi_prev_obj = PPI_preview(score_ppi_prev = score_ddi, type_ppi_prev = type_score, fk_couple = couple.id_couple, fk_prot_bact = id_prot_bact, fk_prot_phage = id_prot_phage)
                list_ppi_preview.append(ppi_prev_obj)
    return list_ppi_preview

#list_couples = get_couples_NCBI_phagesDB()

list_couples = get_couples_Greg()
print(len(list_couples))
print("qwqwqwqwqwqw")

dict_ddi_interact = get_dict_key_interact()
id_key = get_interact_ddi_id(434,14319)
id_key_b = get_interact_ddi_id(14319, 434)
id_key_c = get_interact_ddi_id(14319, 999)

dict_interactions_DB = get_DDI_by_DB_in_dic()

print(len(dict_interactions_DB))
print(len(dict_ddi_interact))
id_bact = -1
id_phage = -1

list_prot_bact = []
dict_dom_prot_bact = {}
list_prot_phage = []
dict_dom_prot_phage = {}
list_all_ppi_preview = []

new_list = []
not_new_list = []

print(len(list_couples))
list_id_fk = PPI_preview.get_all_ppi_preview_couple()
for couple in list_couples:
    print(couple)
    list_all_ppi_preview = []
    if couple.id_couple not in list_id_fk:
        new_list.append(couple.id_couple)
        #for each couple
        if couple.fk_bacteria != id_bact:
            list_prot_bact = Protein.get_all_Proteins_with_domains_by_organism_id(couple.fk_bacteria)
            dict_dom_prot_bact = get_domains_by_list_prots(list_prot_bact)
            id_bact = couple.fk_bacteria

        if couple.fk_phage != id_phage:
            list_prot_phage = Protein.get_all_Proteins_with_domains_by_organism_id(couple.fk_phage)
            dict_dom_prot_phage = get_domains_by_list_prots(list_prot_phage)
            id_phage = couple.fk_phage

        print("It has {0:d} proteins for the bacterium".format(len(list_prot_bact)))
        print("It has {0:d} proteins for the Phages".format(len(list_prot_phage)))


        #parcourir domains by prot de la bact
        for key_id_prot_bact in dict_dom_prot_bact.keys():
            vec_dom_prot_bact = dict_dom_prot_bact[key_id_prot_bact]
            for key_id_prot_phage in dict_dom_prot_phage.keys():
                vec_dom_prot_phage = dict_dom_prot_phage[key_id_prot_phage]
                list_ppi_preview = search_id_interactions(vec_dom_prot_bact, vec_dom_prot_phage, couple.id_couple, key_id_prot_bact, key_id_prot_phage)
                list_all_ppi_preview = list_ppi_preview + list_all_ppi_preview
        print(len(list_all_ppi_preview))
        for ppi_obj in list_all_ppi_preview:
            ppi_obj.create_ppi_preview()

        aux_value = 0
    else:
        print("Couple treated {0:d}".format(couple.id_couple))
        not_new_list.append(couple.id_couple)

np.savetxt('no_interactions.csv', new_list, delimiter=',')
