from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

import json
import numpy as np


from objects_API.CoupleJ import CoupleJson
from objects_API.ProteinJ import ProteinJson
from objects_API.WholeDNAJ import WholeDNAJson
from objects_API.ContigJ import ContigJson


def writeDictProteinError(dict_protein_error:dict, organism_id:int, path_base:str):
    path_base = path_base + 'proteins_'
    writeDictsErrors(dict_protein_error, path_base, organism_id)

def writeDictWholeDNAError(dict_protein_error:dict, organism_id:int, path_base:str):
    path_base = path_base + 'WD_'
    writeDictsErrors(dict_protein_error, path_base, organism_id)

def writeDictContigError(dict_protein_error:dict, organism_id:int, path_base:str):
    path_base = path_base + 'Contig_'
    writeDictsErrors(dict_protein_error, path_base, organism_id)

def writeDictsErrors(dictionaty_error:dict, path_file_write:str, organism_id:int):
    path_write_json = path_file_write + str(organism_id) + '.json'
    with open(path_write_json, 'w') as file:
         file.write(json.dumps(dictionaty_error)) # use `json.loads` to do the reverse


def checkSequenceProteins(list_proteins:ProteinJson):
    """
    check if the sequence exists for a protein given a list of them

    :param list_proteins: list of the proteins

    :type id: array[ProteinJ]

    :return: dictionary with the ids and sequences if they are to short or inexistent
    :rtype: dict[id_prot]:sequence
    """
    dict_proteins_error = {}
    for protein in list_proteins:
        if len(protein.sequence_AA) < 15:
            print('error')
            id_prot = protein.id
            sequence_prot = protein.sequence_AA
            dict_proteins_error[id_prot] = sequence_prot

    return dict_proteins_error

def checkWholeDna(whole_dna):
    """
    check if the sequence exists for a whwhole_dna

    :return: dictionary with the id and sequences if they are to short or inexistent
    :rtype: dict[id_whole_dna]:sequence
    """
    dict_whole_dna_error = {}

    if whole_dna is None:
        dict_whole_dna_error['No whole DNA'] = 'No Whole DNA'

    elif len(whole_dna.sequence_DNA) < 100:
        dict_whole_dna_error[whole_dna.id] = whole_dna.sequence_DNA
    return dict_whole_dna_error

def checkSequenceContigs(list_contigs:ContigJson):
    """
    check if the sequence exists for a contig given a list of them

    :return: dictionary with the ids and sequences if they are to short or inexistent
    :rtype: dict[id_contig]:sequence
    """
    dict_contigs_error = {}

    if list_contigs is None:
        dict_contigs_error['No whole Contig'] = 'No Whole Contig'
    else:
        for contig in list_contigs:
            if len(contig.sequence_DNA) < 15:
                print('error')
                id_contig = contig.id
                contig_sequence = contig.sequence_DNA
                dict_contigs_error[id_contig] = contig_sequence

    return dict_contigs_error





def getAllCouples():
    list_couple = CoupleJson.getAllAPI()
    return list_couple


def getProteinsListByOrganismId(organism_id:int):
    list_proteins = ProteinJson.getByOrganismID(organism_id)
    return list_proteins

def getWholeGenomeByOrganismId(organism_id:int):
    try:
        whole_dna_obj = WholeDNAJson.getByOrganismID(organism_id)
        return whole_dna_obj
    except:
        return None



def getContigsByOrganismId(organism_id:int):
    list_contig = ContigJson.getByOrganismID(organism_id)
    return list_contig


def organismValidation(id_organism:int, path_to_save:str):
    print('start organism {0}'.format(id_organism))

    list_proteins = getProteinsListByOrganismId(id_organism)
    whole_dna_obj = getWholeGenomeByOrganismId(id_organism)
    list_contigs = getContigsByOrganismId(id_organism)
    dict_error_protein = {}
    dict_errors_whole_dna = {}
    dict_errors_contig = {}

    dict_errors_whole_dna = checkWholeDna(whole_dna_obj)
    dict_errors_contig = checkSequenceContigs(list_contigs)
    dict_error_protein = checkSequenceProteins(list_proteins)

    if len(dict_error_protein) > 0:
        writeDictProteinError(dict_error_protein, id_organism, path_write)

    if len(dict_errors_contig) > 0:
        writeDictContigError(dict_errors_contig, id_organism, path_write)

    if len(dict_errors_whole_dna) > 0:
        if 'No whole DNA' not in dict_errors_whole_dna.keys():
            writeDictWholeDNAError(dict_errors_whole_dna, id_organism, path_write)

    print('End organism {0}'.format(id_organism))
    return True

def writeOrganismValidated(path_save:str, list_ids_bacteria:np):
    np.savetxt(path_save, list_ids_bacteria, delimiter=",", fmt='%d')

def readOrganismValidated(path_save:str):
    list_ids_organigms = []

    list_ids_organigms = np.genfromtxt(path_save, delimiter=',')

    return list_ids_organigms

path_organisme_validated = 'error_organisms/organisms_validated.csv'

list_ids_bacteria_validated = readOrganismValidated(path_organisme_validated)

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

list_ids_bacteria = []
list_ids_phage = []

list_couples = getAllCouples()
list_ids_bacteria = [couple_obj.bacterium for couple_obj in list_couples]
list_ids_bacteria = list(set(list_ids_bacteria))

list_ids_phage = [couple_obj.bacteriophage for couple_obj in list_couples]
list_ids_phage = list(set(list_ids_phage))

list_ids_organism_validated = readOrganismValidated(path_organisme_validated)

for phage_id in list_ids_phage:

    path_write = 'error_organisms/phage/'

    if phage_id not in list_ids_organism_validated:
        result_insertion = organismValidation(phage_id, path_write)
        if result_insertion:
            list_ids_organism_validated = np.append(list_ids_organism_validated, phage_id)
            writeOrganismValidated(path_organisme_validated, list_ids_organism_validated)


for bacterium_id in list_ids_bacteria:

    path_write = 'error_organisms/phage/'

    if bacterium_id not in list_ids_bacteria_validated:
        result_insertion = organismValidation(bacterium_id, path_write)
        if result_insertion:
            list_ids_bacteria_validated = np.append(list_ids_bacteria_validated, bacterium_id)
            writeOrganismValidated(path_organisme_validated, list_ids_bacteria_validated)

print('Hello')
print('Hello')