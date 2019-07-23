from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI


from objects_API.ProteinPFAMJ import ProteinPFAMJson
from objects_API.ProteinJ import ProteinJson

from objects_API.CogScoreJ import CogScoreJson

from objects_API.PPICogScoreJ import PPICogScoreJson

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()




list_cog_score = CogScoreJson.getAllAPI()
print(len(list_cog_score))


print(list_cog_score[2])

cog_score_obj = CogScoreJson(score=23.4, cog_method_score=1, cog_interaction_source=1)
cog_score_obj = cog_score_obj.setCogScore()
print(cog_score_obj)


dict_param = {}
dict_param['organism'] = 123
dict_param['organism'] = 123


list_prots = ProteinJson.getProteinByFilterParameter(dict_param)
