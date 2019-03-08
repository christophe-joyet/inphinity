# -*- coding: utf-8 -*-
"""
Created on Tue Mai 15 11:17:32 2018

@author: Diogo
"""

from objects_new.Datasets_configurations_types_new import *
from objects_new.Dataset_config_dataset_new import *

list_ds_Types_config = DS_configurations_types.get_all_configurations_ds_types()


for element in list_ds_Types_config:
    print(element)

ds_conf_ds_obj = Dataset_conf_ds(value_configuration = 50, FK_id_configuration_DCT_DCD = 1, FK_id_dataset_DS_DCD = 1)

ds_conf_ds_obj.create_ds_config_ds()


