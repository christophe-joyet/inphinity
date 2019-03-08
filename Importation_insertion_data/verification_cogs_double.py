from objects_new.COGS_interactions_new import *
from objects_new.COGS_source_interact_new import *
from objects_new.COGS_preview_new import *


import time



print("Start to collect the COGS")
list_of_cog_interact = COGS_interactions.get_all_interaction_cogs_limit(245680, 10000)
print("Finish to collect the COGS")
waiting_time = 0
for couple_cog_inter in list_of_cog_interact:
    waiting_time += 1
    print("COG interaction id: {0}".format(couple_cog_inter.id_cog))
    couple_inverted = COGS_interactions(FK_cog_group_A = couple_cog_inter.FK_cog_group_B, FK_cog_group_B = couple_cog_inter.FK_cog_group_A)
    id_cog_inverted = COGS_interactions.get_cog_interaction_id_by_fks_cogs(couple_inverted.FK_cog_group_A,  couple_inverted.FK_cog_group_B)
    if id_cog_inverted != -1:
        #obtention of the id of the normal COG (not the inverted)
        id_cog_interaction_normal = couple_cog_inter.id_cog
        ###UPDATE HERE
        list_cogs_source_interact = Cog_source_interact.get_all_source_interact_by_fk_if_cog_couple(id_cog_inverted)
        if len(list_cogs_source_interact) > 0:
            for cog_source_interact in list_cogs_source_interact:
                list_cog_preview_invert = COGS_preview.get_all_COGS_preview_couple_give_fk_cog_interaction(cog_source_interact.id_source_interact)
                for elements_cog in list_cog_preview_invert:
                    cog_preview_obj_confirm = COGS_preview(elements_cog.FK_id_prot_bact, elements_cog.FK_id_prot_phage, couple_cog_inter.id_cog, elements_cog.FK_id_couple)
                    qty_exitst = cog_preview_obj_confirm.verify_COG_preview_exist()

                    if qty_exitst == 0:
                        #If any exists it is updated
                        COGS_preview.update_id_cog_interaction_by_cog_preview_id(elements_cog.id_cog_preview, id_cog_interaction_normal)
                    elif qty_exitst == 1:
                        COGS_preview.remove_COG_preview_by_its_id(elements_cog.id_cog_preview)
                        #remove the invert
                    else:
                        print("It exists an error with: {0}".format(elements_cog.id_cog_preview))
                        #error
        
        
        Cog_source_interact.delete_cog_score_interaction_score_give_it_id(id_cog_inverted)

        id_cog = COGS_interactions.remove_cog_by_fk_a_fk_b(couple_cog_inter.FK_cog_group_B, couple_cog_inter.FK_cog_group_A)
        print("id cog original {0}: COA_a: {1} - COG_b: {2}".format(couple_cog_inter.id_cog, couple_cog_inter.FK_cog_group_A, couple_cog_inter.FK_cog_group_B))
    if (waiting_time % 1000 == 0):
        time.sleep(5)
        print("SLEEP")
