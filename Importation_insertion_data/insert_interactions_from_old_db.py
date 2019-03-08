# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 11:08:09 2017

@author: Stage
"""

import csv
import numpy as np
import time

from SQL_obj_new.Organism_sql_new import *

from objects_new.Couples_new import *



path_positive_interactions = "./tables_csv/positive_interactions.csv"
path_negative_interactions = "./tables_csv/negative_interactions.csv"


values_positive = np.loadtxt(open(path_positive_interactions, "rb"), delimiter=",", skiprows=0)
values_negative = np.loadtxt(open(path_negative_interactions, "rb"), delimiter=",", skiprows=0)


#
#
sql_org_obj = Organisms_sql_new()
#
#print(int(values_positive[0][1]))
#
#id_bact = sql_org_obj.get_bacterium_id_by_gi(int(values_positive[0][1]))
#id_phage = values_positive[0][0]
#
#couple_obj = Couple(-1, id_bact, id_phage, 1, 4, True)
#
#id_couple = couple_obj.create_couple()


#print(id_couple)


for elements in values_negative:
    id_bact = sql_org_obj.get_bacterium_id_by_gi(int(elements[1]))
    id_phage = elements[0]
    
    couple_obj = Couple(-1, id_bact, id_phage, 1, 4, False)
    
    id_couple = couple_obj.create_couple()
    
    print(id_couple)
    time.sleep(0.1)

        
