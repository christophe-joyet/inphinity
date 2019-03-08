import pandas as pd
import numpy as np
from objects_API.BacteriumJ import BacteriumJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.StrainJ import StrainJson
from objects_API.SpecieJ import SpecieJson
from objects_API.CoupleJ import CoupleJson

from objects_new.Couples_new import Couple
from objects_new.Organisms_new import Organism
from objects_new.Strains_new import Strain
from objects_new.Species_new import Specie

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI


def readCSVToDataFrame(path_file):
    dataframe_interactions = pd.read_csv(filepath_or_buffer = path_file, sep=',')
    return dataframe_interactions

def getCoupleFromInteractionIdOldDB(interaction_id):
    couple_obj = Couple.get_couples_by_list_id([interaction_id])[0]
    return couple_obj

def getPhageDesignationById(id_phage):
    organism_obj = Organism.get_organism_by_id(id_phage)
    id_strain = organism_obj.fk_strain
    strain_obj = Strain.get_strain_by_id(id_strain)
    strain_designation = strain_obj.designation
    return strain_designation


def getBacteriumStrainSpecieDesignationById(id_bacterium):
    organism_obj = Organism.get_organism_by_id(id_bacterium)
    id_strain = organism_obj.fk_strain
    strain_obj = Strain.get_strain_by_id(id_strain)


    id_specie = strain_obj.fk_specie


    specie_obj = Specie.get_specie_by_id(id_specie)
    specie_designation = specie_obj.designation

    strain_designation = strain_obj.designation

    taxonomy_bacterium = 'Specie: ' + specie_designation + ' Strain: ' + strain_designation

    return taxonomy_bacterium


def getBacteriophageDesignationNewDBById(id_new_phage:int):
    bacteriophage_obj = BacteriophageJson.getByID(id_new_phage)
    phage_designation = bacteriophage_obj.designation

    return phage_designation


def getBacteriumDesignationNewByID(id_new_bacterium:int):
    bacterium_obj = BacteriumJson.getByID(id_new_bacterium)
    strain_id = bacterium_obj.strain
    strain_obj = StrainJson.getByID(strain_id)
    strain_designation = strain_obj.designation

    specie_id = strain_obj.specie
    specie_obj = SpecieJson.getByID(specie_id)
    specie_designation = specie_obj.designation

    taxonomy_bacterium = 'Specie: ' + specie_designation + ' Strain: ' + strain_designation

    return taxonomy_bacterium

def addInteractionsNewDB(interaction_type:bool, bacterium_id:int, phage_id:int, level_id:int, lysis_id:int, persone_responsible:int, source_data_id:int, validity_id:bool):


    couple_obj = CoupleJson(interaction_type = interaction_type, bacteriophage = phage_id, bacterium = bacterium_id, level = level_id, lysis = lysis_id, person_responsible = persone_responsible, source_data = source_data_id, validity = validity_id)


    couple_obj_json = couple_obj.setCouple()
    return couple_obj_json


def writeIdsInserted(id_interaction_old, id_interaction_new):
    str_write = str(id_interaction_old) + ',' + str(id_interaction_new)
    with open("ids_couples_inserted.txt", "a") as myfile:
        myfile.write(str_write)

path_file_name = 'data_ids_interaction_new_old_db.csv'

dataframe_data = readCSVToDataFrame(path_file_name)


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

for index, row in dataframe_data.iterrows():
    #Old couple treatment
    interaction_id = row['interaction_id_old_db']
    couple_obj = getCoupleFromInteractionIdOldDB(interaction_id)
    id_phage = couple_obj.fk_phage
    id_bacterium = couple_obj.fk_bacteria

    phage_designation = getPhageDesignationById(id_phage)
    bacterium_designation = getBacteriumStrainSpecieDesignationById(id_bacterium)

    #New couple treatment
    id_new_phage = row['id_phage']
    id_new_bacterium = row['id_bact']

    phage_designation_new = getBacteriophageDesignationNewDBById(id_new_phage)
    bacterium_designation_new = getBacteriumDesignationNewByID(id_new_bacterium)

    type_interaction = couple_obj.interact_pn
    source_data_id = couple_obj.fk_source_data
    person_resposible_id = 3 #Xavier
    level_interact_id = couple_obj.fk_level_interact
    lysis_inter_id = couple_obj.fk_lysis_inter
    validity_id = 4
    #Couple object insertion
    print('------------------------')
    print('Interaction Informations')
    print('Old phage: {0}'.format(phage_designation))
    print('New phage: {0}'.format(phage_designation_new))
    print('Old bacterium: {0}'.format(bacterium_designation))
    print('New bacterium: {0}'.format(bacterium_designation_new))

    if phage_designation == phage_designation_new and bacterium_designation == bacterium_designation_new:
        couple_obj = addInteractionsNewDB(interaction_type = type_interaction, bacterium_id = id_new_bacterium, phage_id = id_new_phage, level_id = level_interact_id, lysis_id = lysis_inter_id, persone_responsible = person_resposible_id, source_data_id = source_data_id, validity_id = validity_id)
        print(couple_obj)
        id_new_couple = couple_obj.id
        id_old_obj = interaction_id
        writeIdsInserted(id_old_obj, id_new_couple)
    else:
        input_value = input("Insert (1 = Yes; other = No) ")
        if input_value == '1':
            couple_obj = addInteractionsNewDB(interaction_type = type_interaction, bacterium_id = id_new_bacterium, phage_id = id_new_phage, level_id = level_interact_id, lysis_id = lysis_inter_id, persone_responsible = person_resposible_id, source_data_id = source_data_id, validity_id = validity_id)
            print(couple_obj)
            id_new_couple = couple_obj.id
            id_old_obj = interaction_id
            writeIdsInserted(id_old_obj, id_new_couple)



    
print(dataframe_data)