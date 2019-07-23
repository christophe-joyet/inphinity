from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.ProteinPFAMJ import ProteinPFAMJson
from objects_API.ProteinJ import ProteinJson
from objects_API.CoupleJ import CoupleJson
from objects_Api.ProteinJ import ProteinJson


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()


def getAllCouplesByParameters(parameter_dic:dict):
    list_couples = CoupleJson.getCouplesByFilterParameter(parameter_dic)
    return list_couples


def getAllProteinsByOrganism(organism_id:int):
    dict_params_protein = {}
    dict_params_protein['organism_id'] = organism_id
    list_proteins = ProteinJson.getProteinByFilterParameter(dict_params_protein)
    return list_proteins



dict_param_couple = {}
dict_param_couple['level_id'] = 1
list_couples = getAllCouplesByParameters(dict_param_couple)

id_bacterium = list_couples[0].bacteriophage
list_prots_bacteriophage = getAllProteinsByOrganism(id_bacterium)
print(len(list_prots_bacteriophage))

print(len(list_couples))
dict_parameter = {}
dict_parameter['protein'] = 11245

list_proteinPfam = ProteinPFAMJson.getProteinPfamByFilterParameter(dict_parameter)
print('Hello')