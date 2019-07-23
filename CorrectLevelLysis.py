import pandas as pd 

from objects_API.StrainJ import StrainJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.CoupleJ import CoupleJson

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_new.Couples_new import Couple
from objects_new.Organisms_new import Organism

def readCSVToDF(path_csv):
    data = pd.read_csv(path_csv) 
    return data


def getBacteriophageByACCNEWDB(acc_value:str):
    bacteriophage_obj_new = BacteriophageJson.getByAccnumber(acc_value)
    id_bacteriophage_new = bacteriophage_obj_new.id
    return id_bacteriophage_new

def getBacteriumByACCNEWDB(acc_value:str):
    bacterium_obj_new = BacteriumJson.getByAccnumber(acc_value)
    id_bacterium_new = bacterium_obj_new.id
    return id_bacterium_new

def obtainBacteriumIdFromOldDBId(id_bacterium_old_db:int):
    organism_obj = Organism.get_organism_by_id(id_bacterium_old_db)
    fk_strain_old_db = organism_obj.fk_strain
    return fk_strain_old_db

def obtainBacteriumACCnumberFromOldDBId(id_bacterium_old_db:int):
    organism_obj = Organism.get_organism_by_id(id_bacterium_old_db)
    acc_number = organism_obj.acc_num
    return acc_number

def obtainphageACCnumberFromOldDBId(id_bacteriophage_old_db:int):
    organism_obj = Organism.get_organism_by_id(id_bacteriophage_old_db)
    acc_number = organism_obj.acc_num
    return acc_number

def getIdStrainNewDBByStrainBactOldDB(strain_id_old:int, dataframe_strains_id):
    dataframe_line = dataframe_strains_id.loc[dataframe_strains_id['strain_db'] == strain_id_old]
    new_db_id_strain = int(dataframe_line['strain_api'].values[0])
    return new_db_id_strain

def getBacteriumListIdsByStrainId(strain_id:int):
    list_bacterium_ids_treated = []
    strain_obj = StrainJson.getByID(strain_id)
    list_bacterium_ids = strain_obj.bacteria
    for bacterium in list_bacterium_ids:
        bacterium = bacterium.replace('http://trex.lan.iict.ch:8080/api/bacterium/','')[:-1]
        list_bacterium_ids_treated.append(bacterium)
    return list_bacterium_ids_treated

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

path_csv_strains_correspondence = 'correspondenceIDSStrains2.csv'
datafram_csv = readCSVToDF(path_csv_strains_correspondence)

list_couples_old_db = Couple.get_all_couples()
count_error = 0
count_many = 0
count_pos_list = 0
#list_couples_old_db = list_couples_old_db[9602:]
#list_couples_old_db = list_couples_old_db[-2:]
dict_convert_phages_id = {}
dict_convert_phages_id[4656] = 6265

for couple_element in list_couples_old_db:
    if couple_element.fk_source_data == 1 and couple_element.fk_level_interact == 3:
        print('It is the {0} : couple id {1} : pos in list total {2}'.format(count_many, couple_element.id_couple, count_pos_list))
        id_bacterium = couple_element.fk_bacteria
        acc_bacterium = obtainBacteriumACCnumberFromOldDBId(id_bacterium)

        id_new_bacterium_db = -1
        try:
            id_new_bacterium_db = getBacteriumByACCNEWDB(acc_bacterium)
        except:
            id_new_bacterium_db = -1
        #If necessary, check if they have the acc for public data for bacterium
        #strain_id = obtainBacteriumIdFromOldDBId(id_bacterium)
        #strain_id_new_db = getIdStrainNewDBByStrainBactOldDB(strain_id, datafram_csv)
        #list_bacterium_id = getBacteriumListIdsByStrainId(strain_id_new_db)
        #

        id_phage = couple_element.fk_phage
        acc_bacteriophage = obtainphageACCnumberFromOldDBId(id_phage)
        id_new_phage_db = -1
        if id_phage in dict_convert_phages_id:
            id_new_phage_db = dict_convert_phages_id[id_phage]
        else:
            try:
                id_new_phage_db = getBacteriophageByACCNEWDB(acc_bacteriophage)

            except:
                id_new_phage_db = -1
        if id_new_phage_db != -1 and id_new_bacterium_db != -1:
            try:
                couple_obj = CoupleJson.getByBacteriumPhageIds(id_new_bacterium_db, id_new_phage_db)
            except:
                count_error += 1
                interaction_type_cp = couple_element.interact_pn
                id_bacterium_cp = id_new_bacterium_db
                id_phage_cp = id_new_phage_db
                validity_id_cp = 4 #not validate
                level_interaction_cp = 2
                source_data_cp = 1
                person_responsible_cp = 3

                couple_obj_json = CoupleJson(interaction_type = interaction_type_cp,
                                             bacteriophage = id_phage_cp,
                                             bacterium = id_bacterium_cp,
                                             level = level_interaction_cp,
                                             person_responsible = person_responsible_cp,
                                             source_data = source_data_cp,
                                             validity = validity_id_cp)
                couple_obj = couple_obj_json.setCouple()
                print(couple_obj)
                print('INSERTEDDD NEW')
        else:
            count_error += 1
        #print(couple_obj)
        count_many += 1
    count_pos_list += 1

print(len(list_couples_old_db))
print('Hello')