


import re
import datetime

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.DomainJ import DomainJson
from objects_API.DomainInteractionPairJ import DomainInteractionPairJson
from objects_API.DomainInteractionSourceJ import DomainInteractionSourceJson









#Corect from here


def readFileContent(file_path:str, enconding_value = "ISO-8859-1"):
    """
    read the content of a file and convert it into a long string

    :param file_path: file path name
    :param enconding_value: enconding system that you want to use

    :type file_path: string
    :type enconding_value: string 


    :return: a long string with the content read
    :rtype: str

    """
    f = open(path_file, encoding = enconding_value)
    file_content = f.read()
    return file_content


def getPFAMDomainInLine(line_content:str, regex_compile:re):
    """
    apply the regey in order to found the PFAM domains

    :param content_file: content of the line (str string)
    :param regex_expression: regex used to found the PFam domains

    :type content_file: string
    :type regex_expression: string 


    :return: vector with two PFAMS
    :rtype: vec[str,str]

    """
    results_pfam = re.findall(regex_compile, line_content)
    return results_pfam

def findLinesWithDomains(content_file:str, regex_expression:str, dict_domain:dict):
    """
    read line and search the PFAM domains

    :param content_file: content of the line (str string)
    :param regex_expression: regex used to found the PFam domains

    :type content_file: string
    :type regex_expression: string 


    :return: vector with two PFAMS
    :rtype: vec[str,str]

    """
    for record in content_file.split('\n'):
        #for i in record:
        if record.startswith('#=ID'):
            print(record)
            vec_PFAM_pairs = getPFAMDomainInLine(record, regex_expression)
            assert len(vec_PFAM_pairs) == 2
            domains_a_designation = vec_PFAM_pairs[0]
            domains_b_designation = vec_PFAM_pairs[1]
            id_domain_a = checkDomain(domains_a_designation, dict_domain)
            id_domain_b = checkDomain(domains_b_designation, dict_domain)

            id_ddi_pair = checkAddDDIInteractionPairExists(id_domain_a, id_domain_b)
            id_ddi_pair_source = checkAddDDIPairSource(id_ddi_pair, 2)
            print('Id of the new ddi pair {0}'.format(str(id_ddi_pair_source)))
                


def checkAddDDIPairSource(id_ddi_pair:int, id_source:int):
    """
    check if the ddi insteraction pair already exists for 3did source

    :param id_ddi_pair: id of the ddi_interaction_pair
    :param id_source: source id

    :type id_ddi_pair: int
    :type id_source: int 


    :return: id of the ddi interaction source
    :rtype: int

    """

    id_interaction_source_pair = DomainInteractionSourceJson.verifyDDIpairSourceExistence(id_ddi_pair, 2)
    if id_interaction_source_pair == -1:

        actual_date_time = datetime.datetime.now().date()
        ddi_interaction_source_pair = DomainInteractionSourceJson(actual_date_time, id_ddi_pair, id_source)
        ddi_interaction_source_pair = ddi_interaction_source_pair.setDomainInteractionSource()
        id_interaction_source_pair = ddi_interaction_source_pair.id
    return id_interaction_source_pair

def checkAddDDIInteractionPairExists(id_domain_a:int, id_domain_b:int):
    """
    check if a given DDI pair already exits, if not insert it and return it is ID if yes return the ID

    :param id_domain_a: id of the domain a
    :param id_domain_b: id of the domain b

    :type id_domain_a: int
    :type id_domain_b: int 


    :return: id of the DDI pair
    :rtype: id

    """
    ddi_interact_pair = DomainInteractionPairJson(id_domain_a, id_domain_b)
    ddi_interact_pair_id = ddi_interact_pair.verifyDDIpairExistenceID()
    if ddi_interact_pair_id == -1:
        
        ddi_interact_pair = ddi_interact_pair.setDomainInteractionPair()
        ddi_interact_pair_id = ddi_interact_pair.id
    return ddi_interact_pair_id


#Domains part
def checkDomain(domain_designation:str, dict_domains:dict):
    """
    check if a give domain already exists in the database. If not insert it and add to the dictionnary.
    Return the id of the domain

    :param domain_designation: designation of the domains (PFxxxxx)
    :param dict_domains: dictionary with the domains

    :type domain_designation: string
    :type dict_domains: dict


    :return: id of the domain
    :rtype: int

    """

    id_domain = -1
    if domain_designation in dict_domains.keys():
        id_domain = dict_domains[domain_designation]
    else:
        domain_obj = DomainJson(designation = domain_designation)
        domain_obj = domain_obj.setDomain()
        id_domain = domain_obj.id
        dict_domains[domain_designation] = id_domain
        print(domain_obj)
    return id_domain


def convertListDomainToDict(list_domains:list):
    """
    convert a list of domains objects to a dictionary with the designation in the key and the id in value

    :param list_domains: list of domains json objects

    :type list_domains: list[DomainJ]


    :return: dictionary
    :rtype: dict[designation] = id

    """
    dict_domains = {x.designation: x.id for x in list_domains}
    return dict_domains


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

#18005 domains
#id max = 18632


list_domains_db = DomainJson.getAllAPI()
print(list_domains_db[0:10])
dict_domains = convertListDomainToDict(list_domains_db)

regex_expression = 'PF\d{5}'
path_file = '3did_flat_Mar_10_2019_UTF_only_ids.txt'
#path_file = '/home/diogo/Desktop/3did_flat_Mar_10_2019_UTF.dat.txt'
content_file = readFileContent(path_file)

findLinesWithDomains(content_file ,regex_expression, dict_domains)
