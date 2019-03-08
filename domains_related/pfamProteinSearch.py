import datetime
import tempfile
import subprocess
import re
import os
import pandas as pd
import numpy as np
import csv


from objects_API.DomainJ import DomainJson
from objects_API.SourcePFAMJ import SourcePFAMJson
from objects_API.ProteinJ import ProteinJson
from objects_API.ProteinPFAMJ import ProteinPFAMJson

from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI




def getListDomainsInDB():
    """
    get a dictionary with all the pfams in the API DB

    :return: dictionary{PF00000:id}
    :rtype dictionary

    """
    list_domains = DomainJson.getAllAPI()
    dictionary_pfam = {}
    print(list_domains)
    for domain in list_domains:
        dictionary_pfam[domain.designation] = domain.id
    return dictionary_pfam


def insertPFAMInAPI(designation_pfam):
    """
    Insert a new domain in tha API and return it

    :param designation_pfam: name of the domain

    :type designation_pfam: string 

    :return: the domain inserted
    :rtype domainJ

    """
    domain_json = DomainJson(designation=designation_pfam)
    domain = domain_json.setDomain()
    return domain

    print(id_domain)

def getPFAMId(designation_pfam):
    """
    verify if a given domain exists and return its id, if not it is inserted and return its id

    :param designation_pfam: name of the domain

    :type designation_pfam: string 

    :return: id of the domain
    :rtype int

    """

    id_pfam = 0

    if designation_pfam not in dictionary_domains:
        domain = insertPFAMInAPI(designation_pfam)
        dictionary_domains[designation_pfam] = domain.id
    id_pfam = dictionary_domains[designation_pfam]

    return id_pfam
        
        
def getDictProteinsByOrganismID(id_organism:int):
    """
    get a dictionary with all the protein sequence given an organism id

    :param id_organism: id of the organism

    :type id_organism: int 

    :return: dictionary{id_protein:sequence_AA}
    :rtype dictionary

    """

    dictionary_protein_sequence = {}
    list_proteins = ProteinJson.getByOrganismID(660)
    for protein in list_proteins:
        protein_id = protein['id']
        protein_sequence_AA = protein['sequence_AA']

        dictionary_protein_sequence[protein_id] = protein_sequence_AA
    return dictionary_protein_sequence


def createTemporaryFileWithSequence(sequence_AA, protein_id):
    """
    create file that is used in the domain search

    :param sequence_AA: id of the organism
    :param protein_id: id of the organism

    :type sequence_AA: str 
    :type protein_id: int 

    :return: file name
    :rtype str

    """
    file_name = str(protein_id) + '.fasta'
    with open(os.path.join('tmp_prots', file_name), 'w') as file:
        file.write(">no_head\n")
        file.write(sequence_AA)
        file.close()
    return file_name


def executePFAMSearch(file_aa_sequence):
    """
    start the pfam research base on the sequence

    :param file_aa_sequence: file name where the aa sequence is

    :type file_aa_sequence: str 

    :return: output of the command
    :rtype str

    """
    file_aa_sequence = 'tmp_prots/' + file_aa_sequence
    command_line = 'hmmscan --acc --noali /home/diogo.leite/Scripts_diogo/pfam_search/Pfam-A.hmm ' + file_aa_sequence
    proc = subprocess.Popen([command_line], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    print('-------------------')
    
    #output = os.popen(command_line).read()
    return out.decode("utf-8")



def splitPFAMVersion(PFAM_designation):
    """
    remove the version of the PFAM (PFXXXX.YY) remove YY

    :param PFAM_designation: PFAM name

    :type PFAM_designation: str 

    :return: pfam name
    :rtype str

    """
    pfam_name_version = PFAM_designation.split('.')
    assert len(pfam_name_version) == 2
    return pfam_name_version[0]


def validatePfamEValueProteinSeq(pfam_vec, value_vec):
    """
    recive a vector with pfams and values parsed on the command output
    and return the pfam name an evalue

    :param pfam_vec: PFAM name
    :param value_vec: values obtained on the domain search

    :type pfam_vec: str 
    :type value_vec: str 

    :return: vector with pfam and e-value
    :rtype vector[Pfam,evalue]

    """
    vec_results = []
    if len(pfam_vec) == 1 and len(value_vec) > 8:
        pfam_designation = splitPFAMVersion(pfam_vec[0])
        vec_results = [pfam_designation, value_vec[3]]
    return vec_results


def parseDomainsResults(content:str):
    """
    parse the output of the pfam search command and return a vector with
    all the pfams and e-values found

    :param content: output content of the command

    :type content: str 

    :return: vector with pfam and e-value
    :rtype vector[[Pfam,evalue],[Pfam,evalue],[Pfam,evalue],...]

    """
    pfam_domains_re = re.compile('(PF[0-9]+.[0-9]+)', re.IGNORECASE)
    e_value_re = re.compile('-?[\d.]+(?:e-?\d+)?', re.IGNORECASE)
    pfam_value = None
    e_value_value = None
    array_values = []
    vec_pfam_evalue = []
    domains_exist = []

    for line_text in content.splitlines():
        pfam_value = pfam_domains_re.findall(line_text)
        e_value_value = e_value_re.findall(line_text)
        vec_pfam_evalue = validatePfamEValueProteinSeq(pfam_value, e_value_value)

        if  len(vec_pfam_evalue) == 2 and vec_pfam_evalue[0] not in domains_exist :
            domains_exist.append(vec_pfam_evalue[0])
            array_values.append(vec_pfam_evalue)
        elif len(vec_pfam_evalue) >0:
            print('Duplicated domain {0}'.format(vec_pfam_evalue[0]))
    return array_values
        


def insertPFAMProteins(vec_pfam_prot, protein_id):
    """
    Insert pfam-protein into the API

    :param vec_pfam_prot: vec with all the pfams and e-values scores
    :param protein_id: id of the protein on the API

    :type vec_pfam_prot: vector[[Pfam,evalue],[Pfam,evalue],[Pfam,evalue],...] 
    :type protein_id: int 

    """
    data_day = datetime.datetime.now().date()
    for element in vec_pfam_prot:
        pfam = element[0]
        e_value = element[1]

        id_pfam = getPFAMId(pfam)

        protein_pfam_json = ProteinPFAMJson(domain = id_pfam, person_responsible = 5, protein = protein_id, source = 5, date_creation = data_day, e_value = e_value)
        protein_pfam_json.setProteinPFAM()


def readCSVBacteriumID(file_path):
    dataframe_bac_id = pd.read_csv(file_path, dtype={'organism_ptr_id':int, 'strain_id':int})
    return dataframe_bac_id

def saveProteinVerifiedID(path_file, protein_id_vec):
    np.savetxt(path_file, protein_id_vec, delimiter=",", fmt='%i')

def getProteinsIdVerifiedID(path):
    if os.path.isfile(path) == True:
        vec_ids_proteins = np.genfromtxt(path, delimiter=',')
    else:
        return []
    return vec_ids_proteins

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

domains_Start = 'Scores for complete sequence (score includes all domains):'
dictionary_domains = getListDomainsInDB()

#list_protein_PFAM = ProteinPFAMJson.getAllAPI()
#for element in list_protein_PFAM:
#    print(element)

#protein_pfam_json = ProteinPFAMJson(domain = 1, person_responsible = 5, protein = 8078373, source = 1, date_creation = aaa)

#protein_pfam_json = protein_pfam_json.setProteinPFAM()


#print(protein_pfam_json)
#file = open("restls_fasta_dom.txt", "r") 
#content =  file.read() 
#vec_domains = parseDomainsResults(content)

#Get Bacteriums IDS
path_csv = 'inphinity_orm_bacterium.csv'
datafram_ids_bact = readCSVBacteriumID(path_csv)

#Treatment of each bacterium
for index, row in datafram_ids_bact.iterrows():
    path_file_ids = str(index) + '_bacterium.csv'
    vec_ids_proteins = getProteinsIdVerifiedID(path_file_ids)
    id_bacterium = row['organism_ptr_id']

    #Get List of proteins of this bacterium
    dict_proteins = getDictProteinsByOrganismID(id_bacterium)

    #process each protein
    for key, value in dict_proteins.items():
        #Create temp files
        temp_file_name = createTemporaryFileWithSequence(value, key)
        #Research of PFAM domains
        results_pfam_search = executePFAMSearch(temp_file_name)
        #Parse the content rensults of HMMSCAN
        vec_domains = parseDomainsResults(results_pfam_search)
        #Insert PFAM domains for the protein
        insertPFAMProteins(vec_domains, key)
        vec_ids_proteins.append(key)

        #Each 100 proteins save them
        if len(vec_ids_proteins) % 100 == 0:
            saveProteinVerifiedID(path_file_ids, array_values)

        print(vec_domains)

