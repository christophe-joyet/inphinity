#!/usr/bin/python
# -*- coding: utf-8 -*-
import similar_proteins

# Seek similar proteins in organisms 5038,5030,...  with more than 70% of similarity
list_protein = similar_proteins.getSimilarProt([5038,5030,5023,5045,5024], 0.7, 0.99)

# Put the similar sequences in a fasta format file
file = open('essai.fasta', 'a+')
for j in range (len(list_protein[3][0])):
    for k in range (5):
        file.write("\n>proteine " +str(k) + " ensemble " + str(j) + "\n")
        file.write(list_protein[3][0][j][k] + "\n")
            
            