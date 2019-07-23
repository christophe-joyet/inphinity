from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.DomainInteractionPairJ import DomainInteractionPairJson
from objects_API.DomainSourceInformationJ import DomainSourceInformationJson
from objects_API.DomainInteractionSourceJ import DomainInteractionSourceJson

from objects_API.FamilyJ import FamilyJson
from objects_API.GenusJ import GenusJson
from objects_API.StrainJ import StrainJson

from objects_API.CoupleJ import CoupleJson





conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

family_obj = FamilyJson.getByID(147)
genus_obj = GenusJson.getByID(98)
strain_obj = StrainJson.getByID(17144)

print(family_obj)
print(genus_obj)
print(strain_obj)

list_couple = CoupleJson.getAllAPI()
list_family = FamilyJson.getAllAPI()
list_gesy = GenusJson.getAllAPI()
list_strain = StrainJson.getAllAPI()

list_couple_ddi = DomainInteractionPairJson.getAllAPI()
list_db_names = DomainSourceInformationJson.getAllAPI()
list_locationsDDI_source = DomainInteractionSourceJson.getAllAPI()
print(len(list_couple_ddi))


#list_couples = getAllCouples()
print('hello')