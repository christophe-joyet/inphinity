# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 11:18:16 2018

@author: Diogo
"""

from objects_new.Organisms_new import *
from objects_new.Couples_new import *
from objects_new.Species_new import *

list_organisms_phage = Organism.get_all_organisms_by_type(2)

print(len(list_organisms_phage))

list_couple = Couple.get_all_couples_by_phage_id(list_organisms_phage[0].id_organism)


list_id_organism = []
for couple in list_couple:
    list_id_organism.append(couple.fk_bacteria)

print(list_id_organism)
list_specie = []
for fk_bact in list_id_organism:
    strin_id = Specie.get_specie_id_by_organism_id(fk_bact)
    list_specie.append(strin_id)

print(list_specie)

list_bac_id_cant_interact = []

for specie in list_specie:
    print(specie)
    list_bact_species = Organism.get_all_bacteria_by_id_specie(specie.id_specie)
    for organism in list_bact_species:
        list_bac_id_cant_interact.append(organism.id_organism)

print(list_bac_id_cant_interact)

