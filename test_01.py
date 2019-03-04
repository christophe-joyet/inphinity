import matplotlib.pyplot as plt

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.ProteinJ import ProteinJson
from objects_API.WholeDNAJ import WholeDNAJson
from objects_API.ContigJ import ContigJson
from objects_API.FamilyJ import FamilyJson
from objects_API.BacteriophageJ import BacteriophageJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.PersonResponsibleJ import PersonResponsibleJson
from objects_API.StrainJ import StrainJson
from objects_API.SpecieJ import SpecieJson
from objects_API.SourceDataJ import SourceDataJson
from objects_API.OrganismJ import OrganismJson
from objects_API.SourceDataJ import SourceDataJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def getAllCouples():
    list_couple = CoupleJson.getAllAPI()
    return list_couple

def getAllBacteriophages():
    list_bacteriophage = BacteriophageJson.getAllAPI()
    return list_bacteriophage

def getAllBacterium():
    list_bacterium = BacteriumJson.getAllAPI()
    return list_bacterium

def getAllPersonResponsible():
    list_person = PersonResponsibleJson.getAllAPI()
    return list_person

def getAllSourceData():
    list_Source_data = SourceDataJson.getAllAPI()
    return list_Source_data

# Make figure and axes
list_bacterium = getAllBacterium()
list_bacteriophage = getAllBacteriophages()
data = [len(list_bacteriophage), len(list_bacterium)]
data_name = ['Phages', 'Bacteries']

#plot pie
fig1, ax1 = plt.subplots()
ax1.pie(data, labels=data_name, autopct='%1.1f%%', startangle=90)
plt.show()
'''
list_bacterium = []
list_bacteriophage = []

list_couple = getAllCouples()

for couple in list_couple:

  if not couple.bacteriophage in list_bacteriophage:
    list_bacteriophage.append(couple.bacteriophage)

  if not couple.bacterium in list_bacterium:
    list_bacterium.append(couple.bacterium)

data = [len(list_bacteriophage), len(list_bacterium)]
data_name = ['Phages', 'Bacteries']

#plot pie
fig1, ax1 = plt.subplots()
ax1.pie(data, labels=data_name, startangle=90)
plt.show()
'''
'''
fig, axs = plt.subplots(3, 2)

# Create chart to compare number of phages in the database
# and the number of phages who have strain interaction with Bacterium
# bacteriophages = getAllBacteriophages()
#
#=============================================================

#chart to compare types of the interactions between couples
labels = ['type 0', 'type 1', 'type 2', 'type 3', 'type 4']
phages_with_interaction_type = [0] * len(labels)
list_couple = getAllCouples()

for couple in list_couple:
    phages_with_interaction_type[couple.interaction_type] += 1

#delete null value
i = 0
while(i < len(labels)):
  if(phages_with_interaction_type[i] == 0):
    del phages_with_interaction_type[i]
    del labels[i]
    continue
  i += 1

#plot pie
axs[0,0].set_title("Type of the interaction for all the couples")
wedges, texts, autotexts = axs[0,0].pie(phages_with_interaction_type, autopct='%.0f%%')
axs[0,0].legend(wedges, labels,
          title="Type",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))

#=============================================================

#chart to know who has insterted couples
list_person = getAllPersonResponsible()
list_person_name = []

for person in list_person:
  list_person_name.append(person.name)

number_couple_by_person = [0] * len(list_person)

for couple in list_couple:
    for person in list_person:
        if(couple.person_responsible == person.id):
            #increase the number of the couple who was inserted by this person
            number_couple_by_person[person.id - 1] += 1

#delete null value
i = 0
while(i < len(number_couple_by_person)):
  if(number_couple_by_person[i] == 0):
    del number_couple_by_person[i]
    del list_person_name[i]
    continue
  i += 1

axs[0,1].set_title("Person who has insterted couples")
wedges, texts, autotexts = axs[0,1].pie(number_couple_by_person, autopct='%.0f%%')
axs[0,1].legend(wedges, list_person_name,
          title="responsible person",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))

#plt.show()

#=============================================================
#chart source data

list_source_data = getAllSourceData()
list_source_data_id = []
list_source_data_name = []
list_values = []
list_name = []

#list of source id
for source in list_source_data:
  list_source_data_id.append(source.id)
  list_source_data_name.append(source.designation)

#dictionnary {source_id : number of couples}
dic = {key: 0 for key in list_source_data_id}

for couple in list_couple:
  if couple.source_data in dic:
    dic[couple.source_data] += 1

#take all the positive value
for values in dic.values():
  if(values != 0):
    list_values.append(values)

#take the name of the source
for source in dic.keys():
  if(dic.get(source) != 0):
    #append the name of the source
    list_name.append(list_source_data_name[source])

#plot
axs[1,0].set_title("Source Data ")
wedges, texts, autotexts = axs[1,0].pie(list_values, autopct='%.0f%%')
axs[1,0].legend(wedges, list_name,
          title="Source data",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
#plt.show()

#=============================================================
#chart all species there are in couples

list_species = SpecieJson.getAllAPI()
list_species_name = []
list_bact_id = []

for bact in list_couple:
  list_bact_id.append(bact.id)

for specie in list_species:
  list_species_name.append(specie.designation)


#dictionnary {bact_id : number of bact}
dic = {bact_id: 0 for bact_id in list_bact_id}

for couple in list_couple:
  if couple.bacterium in dic:
    dic[couple.bacterium] += 1

list_value = []

for value in dic.values():
  list_value.append(value)

axs[1,1].set_title("Species of bacterium")
axs[1,1].pie(list_value)

#plt.show()
#=============================================================
#chart to know level interaction
taxonomy = ['Family', 'Gender', 'Species', 'Strain']
data = [0] * len(taxonomy)

for couple in list_couple:
    data[couple.level - 1] += 1 

i = 0
while(i < len(data)):
  if(data[i] == 0):
    del data[i]
    del taxonomy[i]
    continue
  i += 1

axs[2,0].set_title("Level of the interaction in terms of taxonomy")
axs[2,0].pie(data, autopct='%.0f%%')
wedges, texts, autotexts = axs[2,0].pie(data, autopct='%.0f%%')
axs[2,0].legend(wedges, taxonomy,
          title="Taxonomy",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
#plt.show()

#=============================================================
#validity of each couple

validity = []
number_bact_validity = []
for couple in list_couple:
  if couple.validity in validity:
    number_bact_validity[validity.index(couple.validity)] += 1
  else:
    validity.append(couple.validity)
    number_bact_validity.append(0)
    number_bact_validity[validity.index(couple.validity)] += 1

axs[2,1].set_title("Validity of Each Couple")
wedges, texts, autotexts = axs[2,1].pie(number_bact_validity, autopct='%.0f%%')
axs[2,1].legend(wedges, validity,
          title="Validity",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))

#plt.show()'''

