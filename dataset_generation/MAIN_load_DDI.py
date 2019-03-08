# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 10:39:47 2018

@author: Diogo
"""

from objects_Pfam.PfamA import *
from objects_Pfam.Pfam_interactions import *

from objects_3did.Domain import *
from objects_3did.DDI_interaction_view import *

from objects_DOMINE.PFAM import *
from objects_DOMINE.INTERACTION import *

from objects_new.DDI_interactions_DB_new import *

from objects_new.Domains_new import *
from objects_new.Ddi_interactions_new import *
from objects_new.DDI_interactions_DB_new import *



#### DOMINE ####
#dom_obj = Domain(database_name = 'INPH_proj_out')
#list_dom_proj = dom_obj.get_all_Domains()
#listDomains = PFAM_ddi.get_all_domains()
#for dom in listDomains:
#    dom_obj = Domain(designation = dom.domain, database_name = 'INPH_proj_out')
#    dom_obj.create_domain()

#dom_obj = Domain(designation = listDomains[0])
#print(dom_obj)

#From 3 did
list_domains_3did = Domain_3did.get_all_domains()
#for dom in list_domains_3did:
#    print(dom)
#    dom_obj = Domain(designation = dom.domain, database_name = 'INPH_proj_out')
#    dom_obj.create_domain()

#From Pfam
list_domains_pfam = Pfam_pfamDB.get_all_pfam_only_Acc()
#for dom in list_domains_pfam:
#    print(dom)
#    dom_obj = Domain(designation = dom.pfam_designation, database_name = 'INPH_proj_out')
#    dom_obj.create_domain()

#Insert iPfam DDI
#list_pfam_interactions = Pfam_interaction.get_all_pfam_interactions()
#for interact in list_pfam_interactions:
#    #print(interact)
#    dom_obj_a = Domain(designation = interact.pfam_a)
#    dom_obj_b = Domain(designation = interact.pfam_b)
#    idomA = dom_obj_a.get_id_domain_by_description()
#    idomB = dom_obj_b.get_id_domain_by_description()
#    interaction_obj = DDI_interaction(FK_id_domain_A = idomA, FK_id_domain_B = idomB)
#    #print(dom_obj_a)
#    #print(dom_obj_b)
#    id_interaction = interaction_obj.create_DDI_interaction_if_not_exists()
#    ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 1)
#    id_interaction_b = ddi_interact_db_obj.create_DDI_interactions_DB()
#    print(id_interaction_b)

#Insert 3DID DDI
list_3did_interactions = DDI_interaction_view.get_all_pfam_interactions()
#for interact in list_3did_interactions:
#    print(interact)
#    dom_obj_a = Domain(designation = interact.domain_A)
#    dom_obj_b = Domain(designation = interact.domain_B)
#    idomA = dom_obj_a.get_id_domain_by_description()
#    idomB = dom_obj_b.get_id_domain_by_description()
#    interaction_obj = DDI_interaction(FK_id_domain_A = idomA, FK_id_domain_B = idomB)
#    print(dom_obj_a)
#    print(dom_obj_b)
#    id_interaction = interaction_obj.create_DDI_interaction_if_not_exists()
#    ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 2)
#    ddi_interact_db_obj.create_DDI_interactions_DB()


#insert for DOMINE
list_DOMINE_interactions = interaction_ddi.get_all_pfam_interactions()
for interact in list_DOMINE_interactions:
    print(interact)
    dom_obj_a = Domain(designation = interact.domain_A)
    dom_obj_b = Domain(designation = interact.domain_B)
    idomA = dom_obj_a.get_id_domain_by_description()
    idomB = dom_obj_b.get_id_domain_by_description()
    interaction_obj = DDI_interaction(FK_id_domain_A = idomA, FK_id_domain_B = idomB)
    print(dom_obj_a)
    print(dom_obj_b)
    id_interaction = interaction_obj.create_DDI_interaction_if_not_exists()
    print(interaction_obj)
    #ipfam
    if interact.iPfam == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 1)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #3did
    if interact.did3 == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 2)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #ME
    if interact.ME == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 3)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #RCDP
    if interact.RCDP == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 4)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #Pvalue
    if interact.Pvalue == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 5)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #Fusion
    if interact.Fusion == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 6)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #DPEA
    if interact.DPEA == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 7)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #PE
    if interact.PE == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 8)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #GPE
    if interact.GPE == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 9)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #DIPD
    if interact.DIPD == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 10)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #RDFF
    if interact.RDFF == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 11)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #KGIDDI
    if interact.KGIDDI == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 12)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #INSITE
    if interact.INSITE == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 13)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #DomainGA
    if interact.DomainGA == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 14)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #PP
    if interact.PP == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 15)
        ddi_interact_db_obj.create_DDI_interactions_DB()
    #SameGO
    if interact.SameGO == 1:
        ddi_interact_db_obj = DDI_interaction_DB(FK_DDI_interaction = id_interaction, FK_DB_source = 16)
        ddi_interact_db_obj.create_DDI_interactions_DB()
