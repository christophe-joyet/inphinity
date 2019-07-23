#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, './')

import pandas as pd
import numpy as np
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.ensemble import ExtraTreesClassifier
import matplotlib.pyplot as plt

data = pd.read_csv("../Bacterium_id_5190")
X = data.iloc[1:,1:31]   # Features
y = data.iloc[1:,0]      # Target column i.e price range
# k is the number of features you want to select
model = ExtraTreesClassifier()
model.fit(X,y)

bestfeatures = SelectKBest(score_func=chi2, k=30)
fit = bestfeatures.fit(X,y)
dfscores = pd.DataFrame(fit.scores_)
dfcolumns = pd.DataFrame(X.columns)
# Concat two dataframes for better visualization 
featureScores = pd.concat([dfcolumns,dfscores],axis=1)
featureScores.columns = ['Features','Score']  # Naming the dataframe columns
print(featureScores.nlargest(30,'Score'))  # Print 10 best features

feat_importances = pd.Series(model.feature_importances_, index=X.columns)
feat_importances.nlargest(30).plot(kind='barh')
plt.show()