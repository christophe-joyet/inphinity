
import datetime
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI


from objects_API.DomainSourceInformationJ import DomainSourceInformationJson
from objects_API.DomainInteractionPairJ import DomainInteractionPairJson

from objects_API.DomainInteractionSourceJ import DomainInteractionSourceJson

from objects_API.DomainJ import DomainJson

from objects_3did.DDI_interaction_view import DDI_interaction_view
from objects_Pfam.Pfam_interactions import Pfam_interaction
from objects_DOMINE.INTERACTION import interaction_ddi

#===============================================
# Script used the domains into the new database
#===============================================


##Test domains Source information
#list_domains_source = DomainSourceInformationJson.getAllAPI()

#print(list_domains_source[0])
#print(list_domains_source[0])


##Test Domains information Pairs
#dom_int_pai_obj = DomainInteractionPairJson(domain_a = 869,domain_b = 887)
##dom_int_pai_obj.setDomainInteractionPair()

#list_interaction_pairs = DomainInteractionPairJson.getAllAPI()
#print(list_interaction_pairs[0])

##Test Domain Interaction Source
#data_day = datetime.datetime.now().date()

#doma_interact_source = DomainInteractionSourceJson(date_creation = data_day, domain_interaction = 28, information_source = 1)

##doma_interact_source.setDomainInteractionSource()

#list_interaction_source = DomainInteractionSourceJson.getAllAPI()
#print(list_interaction_pairs[0])





def get_list_domains():
    """
    This method get all the pfam domains in the django database and put them into a dictionary

    :return dictionary{pfam_designation: id}
    :rtype dictionary

    """
    dictionary_pfam = {}
    list_domains = DomainJson.getAllAPI()
    dictionary_pfam = { pfam_element.designation: pfam_element.id for pfam_element in list_domains }

    return dictionary_pfam


def get_pfam_interactions_3did():
    """
    This method get all the ddi couples in 3did DB

    :return List[tuple(domain_a, domain_b)
    :rtype List(tuple)

    """
    list_interactions = DDI_interaction_view.get_all_pfam_interactions()

    data_couples = [(ddi_couple.domain_A, ddi_couple.domain_B) for ddi_couple in list_interactions]

    return data_couples

def get_pfam_interactions_iPfam():
    """
    This method get all the ddi couples in iPfam DB

    :return List(tuple(domain_a, domain_b))
    :rtype List(tuple)

    """
    list_interactions = Pfam_interaction.get_all_pfam_interactions()

    data_couples = [(ddi_couple.pfam_a, ddi_couple.pfam_b) for ddi_couple in list_interactions]

    return data_couples

def get_pfam_interactions_DOMINE():
    """
    This method get all the ddi couples in DOMINE DB

    :return List of interactions objects
    :rtype List(INTERACTION)

    """
    list_interactions = interaction_ddi.get_all_pfam_interactions()

    return list_interactions


def get_ddi_by_DOMINE_technic_name(source_ddi, list_interactions_DOMINE):

    """
    get list of tuples given a sorouce data in DOMINE


    :param source_ddi: name of the source in DOMINE database
    :param list_interactions_DOMINE: list of INTERACTION objects

    :type source_ddi: str - required   
    :type list_interactions_DOMINE: List(INTERACTION) - required   

    :return List(tuple(domain_a, domain_b))
    :rtype List(tuple)
    """
    print(source_ddi)
    assert source_ddi in ['iPfam', '3did', 'ME', 'RCDP', 'Pvalue', 'Fusion', 'DPEA',
     'PE', 'GPE', 'DIPD', 'RDFF', 'KGIDDI', 'INSITE', 'DomainGA', 'PP',
      'PredictionConfidence', 'SameGO']

    if source_ddi == '3did':
        source_ddi = 'did3'
    list_tuples_interactions = []
    for interactions_value in list_interactions_DOMINE:
        value = getattr(interactions_value, source_ddi)
        if value == 1:
            list_tuples_interactions.append((interactions_value.domain_A, interactions_value.domain_B))
    return list_tuples_interactions


def getListSourceDomainsDict():
    """
    get list of domains sources


    :return dictionary with all the sources and id
    :rtype Dictionary{source_name: id}
    """
    dict_sources = {}
    list_domains_source = DomainSourceInformationJson.getAllAPI()

    data_couples = {source_ddi.designation : source_ddi.id for source_ddi in list_domains_source}

    return data_couples


def insertNewDomain(domain_designation):
    """
    Insert a new domain into django db given it name


    :param domain_designation: name of the domain

    :type domain_designation: str - required   

    :return domain object
    :rtype DomainJson
    """
    domain_obj = DomainJson(designation = domain_designation)
    domain_obj_json = domain_obj.setDomain()
    return domain_obj_json


def insertDDISource(pfam_a, pfam_b, dict_pfams, id_source_information):
    """
    Insert a new DDI into django db.
    1- get the id of both domains if they arn't in the dictionary
    2- verify if the DDI already exists (and the reverse to [A-B - B-A])
    3- if not, insert


    :param pfam_a: name of the domain A
    :param pfam_b: name of the domain B
    :param dict_pfams: dictionary of the domains inserted
    :param id_source_information: source of the information (3did, Ipfam,...)

    :type pfam_a: str - required   
    :type pfam_b: str - required  
    :type dict_pfams: dictionary - required  
    :type id_source_information: Int - required  

    """
    date_day = datetime.datetime.now().date
    if pfam_a not in dict_pfams:
        domain_json = insertNewDomain(pfam_a)
        dict_pfams[domain_json.designation] = domain_json.id

    if pfam_b not in dict_pfams:
        domain_json = insertNewDomain(pfam_b)
        dict_pfams[domain_json.designation] = domain_json.id

    id_ddi_pair = DomainInteractionPairJson.verifyDDIpairExistence(pfam_a,pfam_b)
    if id_ddi_pair == -1:
        id_pfam_a = dict_pfams[pfam_a]
        id_pfam_b = dict_pfams[pfam_b]
        ddi_inserted = DomainInteractionPairJson(domain_a = id_pfam_a, domain_b = id_pfam_b)
        ddi_inserted = ddi_inserted.setDomainInteractionPair()
        id_ddi_pair = ddi_inserted.id

    id_ddi_source_information = DomainInteractionSourceJson.verifyDDIpairSourceExistence(id_ddi_pair, id_source_information)
    if id_ddi_source_information == -1:
        ddi_source_obj = DomainInteractionSourceJson(date_creation = date_day, domain_interaction = id_ddi_pair, information_source = id_source_information)
        print(ddi_source_obj)
        print('B---------------')
        ddi_source_obj = ddi_source_obj.setDomainInteractionSource()

def insertIpfamDDI(list_tuples, dict_elements):
    """
    Insert a new DDI extracted from IpFam

    :param list_tuples: list with all the domains (list of tuples)
    :param dict_elements: dict with all the domains inserted

    :type list_tuples: list(tuple(str, str)) - required   
    :type domain_designation: dictionary - required   

    """
    for dom_a, dom_b in list_tuples:
       insertDDISource(dom_a, dom_b, dict_elements, 1)

def insert3DidPfam(list_tuples, dict_elements):
    """
    Insert a new DDI extracted from 3Did

    :param list_tuples: list with all the domains (list of tuples)
    :param dict_elements: dict with all the domains inserted

    :type list_tuples: list(tuple(str, str)) - required   
    :type domain_designation: dictionary - required   

    """
    for dom_a, dom_b in list_tuples:
        insertDDISource(dom_a, dom_b, dict_elements, 2)

def insertDomineDdi(list_tuples, dict_elements, information_source):
    """
    Insert a new DDI extracted from DOMINE

    :param list_tuples: list with all the domains (list of tuples)
    :param dict_elements: dict with all the domains inserted

    :type list_tuples: list(tuple(str, str)) - required   
    :type domain_designation: dictionary - required   

    """
    for dom_a, dom_b in list_tuples:
        insertDDISource(dom_a, dom_b, dict_elements, information_source)


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()


dict_sources_ddi = getListSourceDomainsDict()

#Verify the DDI existence
dict_values = get_list_domains()
id_ddi_pair = DomainInteractionPairJson.verifyDDIpairExistence('PF00001','PF00001')

#Insert iPfam
list_tuples_iPfam = get_pfam_interactions_iPfam()
#insertIpfamDDI(list_tuples_iPfam, dict_values)


#Insert 3 did
list_tuple_interactions_3did = get_pfam_interactions_3did()
#insert3DidPfam(list_tuple_interactions_3did, dict_values)


list_interactions_DOMINE = get_pfam_interactions_DOMINE()

list_sources_ddi_sources = ['ME','RCDP','Pvalue','Fusion','DPEA','PE','GPE','DIPD','RDFF','KGIDDI','INSITE','DomainGA','PP','SameGO']

#Insert DOMINE
for source_name, id_source_domain in dict_sources_ddi.items():
    list_tuples_by_source = get_ddi_by_DOMINE_technic_name(source_name,list_interactions_DOMINE)
    print(len(list_tuples_by_source))
    insertDomineDdi(list_tuples_by_source, dict_values, id_source_domain)




