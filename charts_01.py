import matplotlib.pyplot as plt

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.StrainJ import StrainJson
from objects_API.SpecieJ import SpecieJson
from objects_API.GenusJ import GenusJson
from objects_API.FamilyJ import FamilyJson
from objects_API.SourceDataJ import SourceDataJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def getAllOfCouples():
        list_couples = CoupleJson.getAllAPI()
        return list_couples

def getAllOfStrain():
        list_strain = StrainJson.getAllAPI()
        return list_strain

def getAllOfGenus():
        list_genus = GenusJson.getAllAPI()
        return list_genus

def getAllOfFamily():
        list_family = FamilyJson.getAllAPI()
        return list_family

def getAllOfSourceData():
        list_source_data = SourceDataJson.getAllAPI()
        return list_source_data

#fig_1, axs = plt.subplots(1,1)
fig, axs = plt.subplots(1,3, figsize=(40, 10))

#get all couples in a list
list_of_couples = getAllOfCouples()
bacterie_found = 0

#dictionnary {bacterie : nbr of bacteries}
bacteries_dictionnary = {}
for couple in list_of_couples:
        bacterie_found = 0
        for bacterie in bacteries_dictionnary.keys():
                #check if there is already the id of the bacterium in the dictionnary
                if couple.bacterium == bacterie.id:
                        #increase de number of the bacterie whith this id
                        bacteries_dictionnary[bacterie] += 1
                        bacterie_found = 1
                        break

        #if we hadn't found the bacterie we add it in the dictionnary
        if bacterie_found != 1 :
               #add the key-value : bacterie
               bacteries_dictionnary[BacteriumJson.getByID(couple.bacterium)] = 1

#dictionnary {bacteriophage_id : nbr of bacteries}
phages_dictionnary = {}
for couple in list_of_couples:
        #check if there is already the id of the phage in the dictionnary
        if not couple.bacteriophage in phages_dictionnary.keys():
                #add the key-value : pahges_id
                phages_dictionnary[couple.bacteriophage] = 0
        #increase de number of the phage whith this id
        phages_dictionnary[couple.bacteriophage] += 1

print("Number of different bacteries : ")
print(len(bacteries_dictionnary))
print("Number of different phages    : ")
print(len(phages_dictionnary))

print("Strain")

#dictionnary {strain_name : number of bacteries in this strain}
strain_dictionnary = {}
strain_found = 0
for bacterie in bacteries_dictionnary.keys():
        #get the strain designation from the bacterie
        strain_id = bacterie.strain
        for strain in strain_dictionnary.keys():
                strain_found = 0
                #check if there is already the id of the strain in the dictionnary
                if strain_id == strain.id:
                         #increase de number of the bacterie whith this id
                                strain_dictionnary[strain] += 1
                                strain_found = 1
                                break
        #if we hadn't found the bacterie we add it in the dictionnary
        if strain_found != 1 :
                #add the key-value : bacterie
                strain_dictionnary[StrainJson.getByID(bacterie.strain)] = 1
        
print(len(strain_dictionnary))

#dictionnary {strain_name : number of bacteries in this strain}
specie_dictionnary = {}
specie_found = 0
for strain in strain_dictionnary.keys():
        #get the strain designation from the bacterie
        specie_id = strain.specie
        for specie in specie_dictionnary.keys():
                specie_found = 0
                #check if there is already the id of the strain in the dictionnary
                if specie_id == specie.id:
                         #increase de number of the bacterie whith this id
                                specie_dictionnary[specie] += 1
                                specie_found = 1
                                break
        #if we hadn't found the bacterie we add it in the dictionnary
        if specie_found != 1 :
                #add the key-value : bacterie
                specie_dictionnary[SpecieJson.getByID(strain.specie)] = 1

print(len(specie_dictionnary)) 

#species repartition charts
number_of_species_tab = []
species_name = []
for specie, number_of_species in specie_dictionnary.items():
        #take only specie with more than one souche
        if number_of_species > 1:
                number_of_species_tab.append(number_of_species)
                species_name.append(specie.designation)

axs[0].set_title("Species repartition in couples")
axs[0].pie(number_of_species_tab, labels=species_name, autopct='%1.0f%%', startangle=140)

#===================================================================================

print("Genus")

#dictionnary {genus : number of species in this genus}
genus_dictionnary = {}
genus_found = 0

for specie in specie_dictionnary.keys():
        #get the strain designation from the bacterie
        genus_id = specie.genus
        for genus in genus_dictionnary.keys():
                genus_found = 0
                #check if there is already the id of the strain in the dictionnary
                if genus_id == genus.id:
                        #increase de number of the bacterie whith this id
                        genus_dictionnary[genus] += 1
                        genus_found = 1
                        break
        #if we hadn't found the bacterie we add it in the dictionnary
        if genus_found != 1 :
                #add the key-value : bacterie
                genus_dictionnary[GenusJson.getByID(specie.genus)] = 1

print(len(genus_dictionnary)) 

'''
genus_dic = {}
genus_list = getAllOfGenus()
genus_dic_list = []

for specie_id in specie_dic:
        #get the genus id from the strain
        genus_id = SpecieJson.getByID(specie_id).genus
        if not genus_id in genus_dic:
                for genus in genus_list:
                        if genus.id == genus_id:
                                #add the genus id in the dictionnary
                                genus_dic[genus_id] = 0
                                #add the genus object in the list
                                genus_dic_list.append(genus)
        #add the number of the species to the specific genus
        genus_dic[genus_id] += specie_dic[specie_id]

#For the plot
fig_3, axs = plt.subplots(1,1)
#reset
names_for_charts = []
values_for_charts = []
names_for_charts.append('Others')
values_for_charts.append(0)

for genus_id in genus_dic:
        
        if(genus_dic[genus_id] > 50):
                names_for_charts.append(genus_dic_list[list(genus_dic.keys()).index(genus_id)].designation)
                values_for_charts.append(genus_dic[genus_id])
        else:
                #add the value to the section 'Others'
                values_for_charts[0] += genus_dic[genus_id]
  
#display plot
plt.title("Genus")
plt.pie(values_for_charts, labels=names_for_charts, autopct='%1.1f%%', startangle=140)

#plt.show()

print("Family")

family_dic = {}
family_list = getAllOfFamily()
family_dic_list = []

for genus_id in genus_dic:
        #get the family id from the genus
        for genus in genus_dic_list:
                if genus_id == genus.id:
                        family_id = genus.family
                        break
        
        if not family_id in family_dic:
                for family in family_list:
                        if family.id == family_id:
                                #add the family id in the dictionnary
                                family_dic[family_id] = 0
                                #add the family object in the list
                                family_dic_list.append(family)
        #add the number of the species to the specific genus
        family_dic[family_id] += genus_dic[genus_id]

#For the plot
fig_4, axs = plt.subplots(1,1)
#reset
names_for_charts = []
values_for_charts = []
names_for_charts.append('Others')
values_for_charts.append(0)

for family_id in family_dic:
        
        if(family_dic[family_id] > 50):
                names_for_charts.append(family_dic_list[list(family_dic.keys()).index(family_id)].designation)
                values_for_charts.append(family_dic[family_id])
        else:
                #add the value to the section 'Others'
                values_for_charts[0] += family_dic[family_id]
  
#display plot
plt.title("Family")
plt.pie(values_for_charts, labels=names_for_charts, autopct='%1.1f%%', startangle=140)

plt.show()'''
#===================================================================================

source_data_list = getAllOfSourceData()
source_data_dictionnary = {}

for source_data in source_data_list:
        source_data_dictionnary[source_data] = 0
        
for bacterie in bacteries_dictionnary.keys():
        #get the strain designation from the bacterie
        source_data_id = bacterie.source_data
        for source_data in source_data_dictionnary.keys():
                if source_data_id == source_data.id:
                        #increase de number of the bacteries managed by this source_data id
                        source_data_dictionnary[source_data] += 1
                        
#For the charts
number_of_bacterie_tab1 = []
source_data_name1 = []
for source_data, number_of_bacterie in source_data_dictionnary.items():
        #take only specie with more than one bacterie
        if number_of_bacterie > 0:
                number_of_bacterie_tab1.append(number_of_bacterie)
                source_data_name1.append(source_data.designation)

#===================================================================================

for source_data in source_data_list:
        source_data_dictionnary[source_data]= 0

for couple in list_of_couples:
        #get the strain designation from the bacterie
        source_data_id = couple.source_data
        for source_data in source_data_dictionnary.keys():
                source_data_found = 0
                if source_data_id == source_data.id:
                        #increase de number of the couple managed by this source_data id
                                source_data_dictionnary[source_data] += 1  
                            
#For the charts
number_of_bacterie_tab2 = []
source_data_name2 = []
for source_data, number_of_bacterie in source_data_dictionnary.items():
        #take only specie with more than one bacterie
        if number_of_bacterie > 0:
                number_of_bacterie_tab2.append(number_of_bacterie)
                source_data_name2.append(source_data.designation)

#===================================================================================

#dictionnary : public couples (couples where source_data = NCBI or PhagesDB)
public_couples_dictionnary = {}

for couple in list_of_couples:
        if couple 