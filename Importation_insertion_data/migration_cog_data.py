# -*- coding: utf-8 -*-
"""
Created on Tue May 22 11:02:46 2018

@author: Diogo Leite
"""

import threading
import time
import concurrent.futures
import multiprocessing

from objects_new.COGS_new import *
from objects_new.Cog_interactions_new import *
from objects_new.COGS_source_interact_new import *
from objects_new.COGS_interactions_new import *


class myThread(threading.Thread):
    def __init__(self, cog_obj):
        threading.Thread.__init__(self)
        self.cog_obj = cog_obj
    def run(self):
      print ("Start cod ID {0}".format(self.cog_obj.id_cog_score))

      insert_cog(self.cog_obj)

      print("End cod ID {0}".format(self.cog_obj.id_cog_score))


def insert_cog(cog_obj):
    print ("Start cod ID {0}".format(cog_obj.id_cog_score))
    FK_source_data = 1
    score_interaction = cog_obj.score_value
    design_cog_a = cog_obj.grp_1
    design_cog_b = cog_obj.grp_2
    id_cog_1 = COG.get_id_COG_by_description(design_cog_a)
    id_cog_2 = COG.get_id_COG_by_description(design_cog_b)

    if id_cog_1 == -1:
        print("create new COG")
        cog_obj = COG(designation = design_cog_a)
        id_cog_1 = cog_obj.create_COG()
    if id_cog_2 == -1:
        print("create new COG")
        cog_obj = COG(designation = design_cog_a)
        id_cog_2 = cog_obj.create_COG()

    if (id_cog_1 != -1 and id_cog_2 != -1):
        existance_cogs = COGS_interactions.get_cog_interaction_id_by_fks_cogs(id_cog_1, id_cog_2)
        if existance_cogs == -1:
            cog_interaction_obj = COGS_interactions(FK_cog_group_A = id_cog_1, FK_cog_group_B = id_cog_2)
            id_interaction_cog = cog_interaction_obj.create_cog_interaction()
            cog_source_score = Cog_source_interact(score = score_interaction, FK_source = FK_source_data, FK_interact = id_interaction_cog)
            cog_source_score.create_cog_score_interaction_score()
    print("End cod ID {0}".format(cog_obj.id_cog_score))
#list_cog_interactions = Cog_Interaction.get_all_score_cogs_designation_by_group(1)
#print(len(list_cog_interactions))
#list_cog_interactions = list_cog_interactions + Cog_Interaction.get_all_score_cogs_designation_by_group(2)
#print(len(list_cog_interactions))
#list_cog_interactions = list(set(list_cog_interactions))
#print(len(list_cog_interactions))

#for element in list_cog_interactions:
#    cog_obj = COG(designation=element)
#    cog_obj.create_COG()


##### WHEN FINISH UNCOMENT AND STAR THIS PART
print("Load from the database")
list_cog_interactions_all = Cog_Interaction.get_all_score_cogs()

aux = 0

print("Start the insertion")

nb_processor = multiprocessing.cpu_count()
nb_processor -=2


with concurrent.futures.ThreadPoolExecutor(max_workers=nb_processor) as executor:
    future_to_url = {executor.submit(insert_cog, cog_obj): cog_obj for cog_obj in list_cog_interactions_all}