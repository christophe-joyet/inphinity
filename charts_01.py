import matplotlib.pyplot as plt

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.StrainJ import StrainJson
from objects_API.SpecieJ import SpecieJson
from objects_API.GenusJ import GenusJson
from objects_API.FamilyJ import FamilyJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def getAllOfCouples():
    list_couples = CoupleJson.getAllAPI()
    return list_couples

def getAllOfStrain():
    list_strain = StrainJson.getAllAPI()
    return list_strain

#for the plot
data_name  = []
data_value = []
fig_1, axs = plt.subplots(1,1)

print("bacterium")
list_of_couples = getAllOfCouples()

#dictionnary bacterie_id : number of this bacterie in all couples
bact_dic = {}

for couple in list_of_couples:
        #check if there is already the id of the bacterium in the dictionnary
        if not couple.bacterium in bact_dic.keys():
                #add the key-value : bacterie_id - number of this bacterie
                bact_dic[couple.bacterium] = 0
        #increase de number of the bacterie who has this id
        bact_dic[couple.bacterium] += 1

print("strain")

#dictionnary strain_name : number of bacteries in this strain
strain_dic = {}

for bact_id in bact_dic.keys():
        #get the strain designation from the bacterie
        strain_id = BacteriumJson.getByID(bact_id).strain
        if not strain_id in strain_dic:
                #add the object strain in the dictionnary
                strain_dic[strain_id] = 0
        #add the number of the bacteries to the specific strain
        strain_dic[strain_id] += bact_dic[bact_id]

#For the plot
data_name.append('Others')
data_value.append(0)

for strain in strain_dic:
        #keep only strain with more than 50 bacteries
        if(strain_dic[strain] > 50):
                data_name.append(StrainJson.getByID(strain).designation)
                data_value.append(strain_dic[strain])
        else:
                #add the value to the section 'Others'
                data_value[0] += strain_dic[strain]
        
#display plot
plt.pie(data_value, labels=data_name, autopct='%1.1f%%', startangle=140)
#plt.show()

specie_dic = {}

for strain_id in strain_dic:
        #get the specie designation from the strain
        specie_id = StrainJson.getByID(strain_id).specie
        if not specie_id in specie_dic:
                #add the object specie in the dictionnary
                specie_dic[specie_id] = 0
        #add the number of the bacteries to the specific specie
        specie_dic[specie_id] += strain_dic[strain_id]

#For the plot
fig_2, axs = plt.subplots(1,1)
#reset
data_name = []
data_value = []
data_name.append('Others')
data_value.append(0)

for specie in specie_dic:
        #keep only specie with more than 30 strain
        if(specie_dic[specie] > 30):
                data_name.append(SpecieJson.getByID(specie).designation)
                data_value.append(specie_dic[specie])
        else:
                #add the value to the section 'Others'
                data_value[0] += specie_dic[specie]
  
#display plot
plt.pie(data_value, labels=data_name, autopct='%1.1f%%', startangle=140)

plt.show()

'''
#to get all differents Bacteries' id
list_bacterium_id = []
#to 
list_number_of_bacterium = []
list_strain = getAllOfStrain()
list_species = SpecieJson.getAllAPI()
print(len(list_strain))
print(len(list_species))

for couple in list_of_couples:
        #check if there is already the id of the bacterium in the list
        if not couple.bacterium in list_bacterium_id:
                #add the id of de bacterium
                list_bacterium_id.append(couple.bacterium)
                list_number_of_bacterium.append(0)
        list_number_of_bacterium[list_bacterium_id.index(couple.bacterium)] += 1
         
print(list_bacterium_id)
print(list_number_of_bacterium)
sum = 0
for x in list_number_of_bacterium:
        sum += x
print(sum)
# Strain

print("strain")
list_strain_in_couples = []
list_count_strain = []
list_designation_strain = []
list_strain_id = []

for bacterium_id in list_bacterium_id:
        strain_id = BacteriumJson.getByID(bacterium_id).strain
        #check if there is already the strain in the list
        if not strain_id in list_strain_id:
          #add the strain in the list
                list_strain_id.append(strain_id)
                strain = StrainJson.getByID(strain_id)
                list_strain_in_couples.append(strain)
                list_count_strain.append(0)
        list_count_strain[list_strain_in_couples.index(strain)] += list_number_of_bacterium[list_bacterium_id.index(bacterium_id)]

for strain in list_strain_in_couples:
        list_designation_strain.append(strain.designation)


#plot Strain
axs[0,0].set_title("Strain")
wedges, texts, autotexts = axs[0,0].pie(list_species_count, autopct='%.0f%%')
axs[0,0].legend(wedges, list_designation_strain,
          title="Source data",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
plt.show()

# Species

print("Species")
list_species_in_couples = []
list_species_count = []
list_species_id = []


for strain in list_strain_in_couples:
        if not strain.specie in list_species_id:
                list_species_id.append(strain.specie)
                specie = SpecieJson.getByID(strain.specie)
                list_species_in_couples.append(specie)
                list_species_count.append(0)
        list_species_count[list_species_in_couples.index(specie)] += list_count_strain[list_strain_in_couples.index(strain)]

list_species_designation = []

for species in list_species_in_couples:
        list_species_designation.append(species.designation)

print(list_species_designation)
print(list_species_count)


# Genus
print("Genus")

list_genus_in_couples = []
list_genus_count = []
list_genus = GenusJson.getAllAPI()

for specie in list_species_in_couples:
        if not specie.genus in list_genus_in_couples:
                #add genus id
                list_genus_in_couples.append(specie.genus)
                list_genus_count.append(0)
        list_genus_count[list_genus_in_couples.index(specie.genus)] += list_species_count[list_species_in_couples.index(specie)]

list_genus_designation = []

for species in list_species_in_couples:
        list_species_designation.append(species.designation)

print(len(list_species_in_couples))
print(len(list_species_count))
'''