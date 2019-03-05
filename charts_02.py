import matplotlib.pyplot as plt

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.OrganismJ import OrganismJson
from objects_API.SourceDataJ import SourceDataJson
from objects_API.ProteinJ import ProteinJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.ContigJ import ContigJson
from objects_API.WholeDNAJ import WholeDNAJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

#for the plots
data_name = []
data_values = []

list_couple = CoupleJson.getAllAPI()
list_oragnism_bacteries = []
list_organism_phages = []
protein_list = []
protein_dic = {}

#get bacteries from the couple
for couple in CoupleJson.getAllAPI():
    if not couple.bacterium in list_oragnism_bacteries:
        list_oragnism_bacteries.append(couple.bacterium)

#get phages from the couple
for couple in CoupleJson.getAllAPI():
    if not couple.bacteriophage in list_organism_phages:
        list_organism_phages.append(couple.bacteriophage)

'''    
print(len((list_oragnism)))
i = 0
for organism_id in list_oragnism:
    print(i)
    i += 1
    protein_list = ProteinJson.getByOrganismID(organism_id)
    for protein in protein_list:
        if not protein.description in protein_dic.keys():
            protein_dic[protein.description] = 0
        protein_dic[protein.description] += 1

#From where the source comes
print("Source Data")

data_name = []
data_values = []

source_data_dic = {}
list_source_data = SourceDataJson.getAllAPI()

for couple in list_couple:
    if not couple.source_data in source_data_dic.keys():
        source_data_dic[couple.source_data] = 0
    source_data_dic[couple.source_data] += 1
    print(BacteriumJson.getByID(couple.bacterium).source_data)

for source_data_id in source_data_dic.keys():
    for source in list_source_data:
        if source.id == source_data_id:
            data_name.append(source.designation)
            data_values.append(source_data_dic[source_data_id])

#plot
fig_1, axs = plt.subplots()
axs.set_title("Source Data")
wedges, texts, autotexts = axs.pie(data_values, autopct='%.0f%%')
axs.legend(wedges, data_name,
          title="Source data",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))'''

#Length of contig
dna_lengths = []
i = 0
while(i < len(list_oragnism_bacteries)):
    print(len(WholeDNAJson.getByOrganismID(list_oragnism_bacteries[i]).sequence_DNA))
    dna_lengths.append(len(WholeDNAJson.getByOrganismID(list_oragnism_bacteries[i]).sequence_DNA))
    i += 1
    
list.sort(dna_lengths)

