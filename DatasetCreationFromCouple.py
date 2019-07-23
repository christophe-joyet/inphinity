from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.CoupleJ import CoupleJson
from objects_API.ProteinJ import ProteinJson
from objects_API.ProteinPFAMJ import ProteinPFAMJson
from objects_API.DomainInteractionPairJ import DomainInteractionPairJson
from objects_API.DomainInteractionSourceJ import DomainInteractionSourceJson

import pickle
import os
import numpy
from numpy import genfromtxt

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()



def getCouplesLevelOne():
    """
    Return all the couples with an interaction of level 1

    :return: list of couples
    :rtype: list(CoupleJson)
    """
    dict_parameters_couple = {}
    dict_parameters_couple['level'] = 1
    list_couples = CoupleJson.getCouplesByFilterParameter(dict_parameters_couple)

    return list_couples


def getProteinsIdsByOrganism(id_organism:int):
    """
    get all the proteins ids given an organism id

    :param id_organism: id of the organism

    :type id_organism: int 


    :return: array with the proteins ids
    :rtype: array[int]

    """

    list_proteins = ProteinJson.getByOrganismID(id_organism)
    list_proteins_id = [protein.id for protein in list_proteins]
    return list_proteins_id


def getDomainsIdsByProteinId(id_protein:int):

    """
    get all the pfam ids in a given protein id

    :param id_protein: id of the protein

    :type id_protein: int 


    :return: array with the domains ids
    :rtype: array[int]

    """
    dict_parameters_pfamProt = {}
    dict_parameters_pfamProt['protein'] = id_protein
    list_domains_protein = ProteinPFAMJson.getProteinPfamByFilterParameter(dict_parameters_pfamProt)
    list_domains_id = [proteinPfam.domain for proteinPfam in list_domains_protein]
    return list_domains_id

def getDictIdsProtsIdDomains(id_organism:int):

    dict_prots_domains = {}
    list_proteins_ids = getProteinsIdsByOrganism(id_organism)
    for id_prot in list_proteins_ids:
        list_domains_ids = getDomainsIdsByProteinId(id_prot)
        dict_prots_domains[id_prot] = list_domains_ids
    return dict_prots_domains

def saveDictIntoPick(path_file:str, dictionary_data:dict):
    """
    save the dictionary into a pickle file

    :param path_file: path of the picke file
    :param dictionary_data: dictionary with the data that you want to save

    :type path_file: str 
    :type dictionary_data: str 

    """
    pickle_out = open(path_file,"wb")
    pickle.dump(dictionary_data, pickle_out)
    pickle_out.close()



def loadDictFromPIckle(path_file:str):
    """
    load the dictionary from a file

    :param path_file: path of the picke file

    :type path_file: str 

    :return: dictionary with the ddis sources ids
    :rtype: dict[(domain_id_a, domain_id_b)] : id sources

    """
    pickle_in = open(path_file,"rb")
    dict_ddi_source_ids = pickle.load(pickle_in)
    return dict_ddi_source_ids




def getDomainsPairIds(path_file:str):
    """
    get all the ddis pair sources. 
    Save the dictionary in a pickle file in case of "down" problems

    :param path_file: path of the picke file

    :type path_file: str 

    :return: dictionary with the ddis sources ids
    :rtype: dict[(domain_id_a, domain_id_b)] : id sources

    """
    dict_ddi_source_ids = {}
    list_domains_pairs_ids = DomainInteractionPairJson.getAllAPI()
    dict_domains_pairs = {(ddi.domain_a, ddi.domain_b) : ddi.id for ddi in list_domains_pairs_ids}
    for key, value in dict_domains_pairs.items():
        dict_parameters_ddi_source = {}
        dict_parameters_ddi_source['domain_interaction'] = value
        list_ddis_test = DomainInteractionSourceJson.getDomainInteractionSourceByFilterParameter(dict_parameters_ddi_source)
        list_ids_sources = [ddi_source.information_source for ddi_source in list_ddis_test]
        dict_ddi_source_ids[key] = list_ids_sources


    saveDictIntoPick(path_file, dict_ddi_source_ids)
    return dict_ddi_source_ids

def loadDictDDIs(path_file:str):
    """
    load the ddi source if the pickle file already existe

    :param path_file: path of the picke file

    :type path_file: str 

    :return: dictionary with the ddis sources ids
    :rtype: dict[(domain_id_a, domain_id_b)] : id sources

    """
    dict_ddis_source_data = {}
    exists = os.path.isfile(path_file)
    if exists:
        dict_ddis_source_data = loadDictFromPIckle(path_file)
    else:
        dict_ddis_source_data = getDomainsPairIds(path_file)
    return dict_ddis_source_data

def loadCSVIdsCouples(path_file:str):
    """
    load the ids of the couples already treated

    :param path_file: path of the picke file

    :type path_file: str 

    :return: numpy array
    :rtype: array[ints]

    """
    ids_couples = []
    exists = os.path.isfile(path_file)
    if exists:
        ids_couples = genfromtxt(path_couples_calulated, delimiter=',')
    else:
        ids_couples = []
    return ids_couples


def searchPFAMmatch(list_pfamId_orga_A:list, list_pfamId_orga_B:list, dict_ddi_sources:dict):
    """
    Combine all the pfam ids and check if they exist in the dictionary of ddi, if yes the tuple was add to an array and returned.
    This method check the id_pfam_A - id_pfam_B AND id_fam_B - id_pfam_A

    :param list_pfamId_orga_A: list of pfam ids in the organism A
    :param list_pfamId_orga_B: list of pfam ids in the organism B
    :param dict_ddi_sources: dictionary that contain all the pairs

    :type list_pfamId_orga_A: list 
    :type list_pfamId_orga_B: list 
    :type dict_ddi_sources: dictionary 

    :return: array with the tuples of id ddi matches
    :rtype: array[(id_domain_a, id_domain_b)]

    """
    array_tuples_match = []
    for id_pfam_orga_a in list_pfamId_orga_A:
        for id_pfam_orga_b in list_pfamId_orga_B:
            if (id_pfam_orga_a, id_pfam_orga_b) in dict_ddi_sources:
                tuple_key = (id_pfam_orga_a, id_pfam_orga_b)
                array_tuples_match.append(tuple_key)
                continue
            elif (id_pfam_orga_b, id_pfam_orga_a) in dict_ddi_sources:
                tuple_key = (id_pfam_orga_b, id_pfam_orga_a)
                array_tuples_match.append(tuple_key)

    return array_tuples_match

def scoreDDICalculation(list_of_ddi_ids_tuples_pair:list, dict_ddi_scores:dict):
    """
    Calculate the PPI score based on sum of the DDI scores
    IF the DDI score comes from iPfam or 3DID I put automaticly 9 which correspond to the maximum of DDI score (combination between all the predictide methods)

    :param list_of_ddi_ids_tuples_pair: list of tuples DDI between two PPI
    :param dict_ddi_scores: dictionary with all the DDI sources

    :type list_of_ddi_ids_tuples_pair: list of tuples 
    :type dict_ddi_scores: dict 


    :return: the PPI score
    :rtype: int

    """
    score_ppi = 0
    score_ppi_aux = 0
    for tuple_ddi in list_of_ddi_ids_tuples_pair:
        list_sources = dict_ddi_scores[tuple_ddi]
        if 1 in list_sources or 2 in list_sources:
            score_ppi_aux = 9
        else:
            score_ppi_aux = len(list_sources)

        score_ppi += score_ppi_aux
        score_ppi_aux = 0
    return score_ppi

def calculatePPIscoresForCouple(dict_proteins_pfam_bact:dict, dict_proteins_pfam_phage:dict, dict_ddi_sources:dict):
    dict_scores_results = {}
    for key_bact, value_bact in dict_proteins_pfam_bact.items():
        for key_phage, value_phage in dict_proteins_pfam_phage.items():
            list_tuples_ddi = searchPFAMmatch(value_bact, value_phage, dict_ddi_sources)
            score_ppi = scoreDDICalculation(list_tuples_ddi, dict_ddi_sources)
            if score_ppi in dict_scores_results:
                dict_scores_results[score_ppi] += 1
            else:
                dict_scores_results[score_ppi] = 1

    return dict_scores_results


dict_protein_pfam_ids = getDictIdsProtsIdDomains(4314)
print('hello')

#load ddi scores
dict_ddis_ids = loadDictDDIs('files_data/picke_ddi_source.pickle')
#load couples
list_couples = getCouplesLevelOne()

#performe the scores
#path file with couples calculated
path_couples_calulated = 'files_data/ids_couples.csv'


#Test load couples dict
#dict_ddis_ids = loadDictFromPIckle('files_data/couple_991.p')


ids_couples_treated = loadCSVIdsCouples(path_couples_calulated)
id_bacterium = 0
id_phage = 0
for couple_obj in list_couples:
    if couple_obj.id not in ids_couples_treated:
        print('I calculate {0}'.format(str(couple_obj.id)))
        if id_bacterium != couple_obj.bacterium:
            id_bacterium = couple_obj.bacterium
            dict_protein_pfam_ids_bacterium = getDictIdsProtsIdDomains(id_bacterium)

        id_phage = couple_obj.bacteriophage
        dict_protein_pfam_ids_bacteriophage = getDictIdsProtsIdDomains(id_phage)

        dict_scores_frequency = calculatePPIscoresForCouple(dict_protein_pfam_ids_bacterium, dict_protein_pfam_ids_bacteriophage, dict_ddis_ids)

        path_save_couple_dict = 'files_data/couple_' + str(couple_obj.id) + '.p'
        saveDictIntoPick(path_save_couple_dict, dict_scores_frequency)
    ids_couples_treated.append(couple_obj.id)
    numpy.savetxt(path_couples_calulated, ids_couples_treated, delimiter=",")


dict_protein_pfam_ids_bacterium = getDictIdsProtsIdDomains(5859)
dict_protein_pfam_ids_phage = getDictIdsProtsIdDomains(5326)
print('hello')

dict_scores_frequency = calculatePPIscoresForCouple(dict_protein_pfam_ids_bacterium, dict_protein_pfam_ids_phage, dict_ddis_ids)
print('Hello')