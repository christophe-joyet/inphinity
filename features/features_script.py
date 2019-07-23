#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, './')
from features import features_functions

features_functions.createFeaturesFile('Bacterium_id_5190.csv', './features/CSV_files', 5190)