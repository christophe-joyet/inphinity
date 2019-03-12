from objects_API.FamilyJ import FamilyJson
from objects_API.GenusJ import GenusJson

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()


family_obj = FamilyJson.getByID(147)
print(family_obj)

genus_obj = GenusJson.getByID(200)
print(genus_obj)
print("Gene id: {0}".format(genus_obj.id))
print('hello')