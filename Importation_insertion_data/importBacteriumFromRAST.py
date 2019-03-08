from files_treatment_new.xls_gen_bank_rast import Xls_gen_bank
from files_treatment_new.fasta_contigs_RAST import Fasta_contigs_RAST

from files_treatment_new.genbank_file_RAST import Genbank_proteic_RAST

from files_treatment_new.fasta_contigs_RAST import Fasta_contigs_RAST


import glob
import os 


from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from files_treatment_new.xls_gen_bank_rast import Xls_gen_bank
from files_treatment_new.fasta_contigs_RAST import Fasta_contigs_RAST

from objects_new.Contigs_new import Contig

from objects_API.ContigJ import ContigJson
from objects_API.StrainJ import StrainJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.GeneJ import GeneJson
from objects_API.ProteinJ import ProteinJson

from Patric.ImportFiles import ImportFilesPatric

def check_file_exits(file_path):
    """
    This method just verify if a given file exists (it is necessary to give the complete path)

    :param file_path: complete path of the file

    :type file_path: string - required

    :return True or False according the existance
    :rtype boolean

    """
    file_exists = os.path.exists(file_path)
    if file_exists is True:
        return True
    else:
        return False

def get_list_ids_files_in_path(path):
    """
    This method list all files in a given path and return a list with these names

    :param path: path where it necessary to list the files

    :type path: string - required

    :return list with the files paths
    :rtype list(str)

    :note when the start point is smaller than end point (int the contig), it is because the "Strand field int excel file is negative
    """

    current_path = os.getcwd() + path
    list_files = os.listdir(current_path)

    return list_files

def createContig(contigObj, organismID):
    """
    insert a Contig into a REST API

    :param contigObj: Contig DBA object that you want to insert
    :param organismID: ID of the organism which has this wholeDNA
    :param listProtein: List of the proteins of the contig

    :type contigObj: WholeDNA
    :type organismID: int
    :type listProtein: List[int]

    :return: id of the Contig inserted
    :rtype int
    """
    contigObjJson = ContigJson(id_db_online = contigObj.id_contig_db_outside, sequence_DNA= contigObj.sequence, fasta_head = contigObj.head, organism = organismID)
    contigObjJson = contigObjJson.setContig()
    return contigObjJson['id']

def createContigNew(contigObj, bacteriumId):
    head_cnt = ''
    if '>' not in contigObj.head:
        head_cnt = '>' + contigObj.head
    else: 
        head_cnt = contigObj.head

    contigObj = Contig(id_contig_db_outside = contigObj.id_contig_db_outside, head = head_cnt, sequence = contigObj.sequence)
    idContig = createContig(contigObj, bacteriumId)
    return idContig

def createGene(id_bacterium, dna_sequence, start_contig, end_contig, fk_contig, id_db_online, function = None, fasta_head = None):
    """
    insert a Gene into a REST API

    :param id_bacterium: ID of the organisms
    :param dna_sequence: DNA sequence of the gene
    :param start_contig: start position of the gene in the contig
    :param end_contig: end position of the gene in the contig
    :param fk_contig: id of the contig
    :param function: function of the gene
    :param fasta_head: fasta head of the gene

    :type id_bacterium: int
    :type dna_sequence: str
    :type start_contig: int - can be None
    :type end_contig: int - can be None
    :type fk_contig: int - can be None
    :type function: str - can be None
    :type fasta_head: str - can be None

    :return: id of the gene inserted
    :rtype int
    """
    geneObjJson = GeneJson(sequence_DNA = dna_sequence, organism = id_bacterium,  position_start_contig = start_contig, position_end_contig = end_contig, contig = fk_contig, fasta_head = fasta_head, id_db_online = id_db_online)
    geneObjJson = geneObjJson.setGene()
    return geneObjJson.id

def createProtein(id_db_online, fk_organism, fk_gene, sequence_aa, description):
    """
    insert a Protein into a REST API

    :param proteinOBJ: Protein DBA object that you want to insert
    :param fkOrganism: id of the organism
    :param fkGene: id of the gene

    :type proteinOBJ: Protein
    :type fkOrganism: int
    :type fkGene: int

    :return: id of the protein inserted
    :rtype int
    """
    proteinObjJson = ProteinJson(id_db_online = id_db_online, organism = fk_organism, gene = fk_gene, sequence_AA = sequence_aa, description = description)
    proteinObjJson = proteinObjJson.setProtein()
    return proteinObjJson.id

def createAndInsertElements(contig_obj, id_bacterium, xls_genbank_rast_obj):
    #listProts = xls_genbank_patric_obj.get_proteins_information_in_excel()
    list_proteins = xls_genbank_rast_obj.get_proteins_objects_by_contig_id(contig_obj.head)
    list_proteins_ids = xls_genbank_rast_obj.get_proteins_ids_by_contig_id(contig_obj.head)

    contig_id = createContigNew(contig_obj, id_bacterium)

    for protein_obj in list_proteins:
        gene_function = None
        gene_id_db_online = None
        if protein_obj.sequence_prot == None:
            gene_function = protein_obj.description
            gene_id_db_online = protein_obj.id_accession
        fasta_head_gene = '>' + protein_obj.id_accession
        id_gene = createGene(id_bacterium, protein_obj.sequence_dna, protein_obj.start_point_cnt, protein_obj.end_point_cnt, fk_contig = contig_id, function = gene_function, fasta_head = fasta_head_gene, id_db_online = gene_id_db_online)
        if protein_obj.sequence_prot != None and len(protein_obj.sequence_prot) > 5:
            createProtein(id_db_online = protein_obj.id_accession, fk_organism = id_bacterium, fk_gene = id_gene, sequence_aa = protein_obj.sequence_prot, description = protein_obj.designation)



def createStrain(designation, fk_specie):
    #Information for the Strain
    strain_obj = StrainJson(designation = designation, specie = fk_specie)
    id_strain = strain_obj.setStrain()
    return id_strain

def createBacterium(acc_number, person_responsible, source_data, fk_strain):

    bacteriumObjJson = BacteriumJson(acc_number = acc_number, person_responsible = person_responsible, source_data = source_data, strain = fk_strain)
    bacteriumObjJson = bacteriumObjJson.setBacterium()

    return bacteriumObjJson.id
         

#Token connection
conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()
#End token connection


#Uncomment this line for the first insertion
list_files = get_list_ids_files_in_path('/RAST/Xls/')
list_files_error = []
list_files_done = []

cwd = os.getcwd()

dict_names = {}
dict_names['470.1135'] = 'Acinetobacter baumannii MDR-ZJ06'
dict_names['470.1134'] = 'Acinetobacter baumannii MDR-TJ'


for file_name in list_files:
    if file_name[:-4] in dict_names:
        path_file_xls = cwd + '/RAST/Xls/' + file_name[:-3] + 'xls'
        path_file_contig = cwd + '/RAST/CONTIG/' + file_name[:-3] + 'contigs.fa'

        xls_file_exists = check_file_exits(path_file_xls)
        contig_file_exist = check_file_exits(path_file_contig)

        assert xls_file_exists == True and contig_file_exist == True, 'A file is missing'

        xls_obj = Xls_gen_bank(path_file_xls, sheet_name = 'Sheet1')
        contig_fasta_file_patric_obj = Fasta_contigs_RAST(path_file = path_file_contig)

        qty_cntg = contig_fasta_file_patric_obj.get_qty_of_contigs()
        list_cnt = contig_fasta_file_patric_obj.create_contigs_from_file()




        gi_name = "Greg_" + dict_names[file_name[:-4]]
        acc_value = "Greg_" + dict_names[file_name[:-4]]

        person_responsible = 2
        source_data = 3

        fk_specie = 613


        strain = createStrain(dict_names[file_name[:-4]], fk_specie)
        id_strain = strain.id

        id_bacterium = createBacterium(acc_value, person_responsible, source_data, id_strain)


        for contig_old in list_cnt:
            createAndInsertElements(contig_old, id_bacterium, xls_obj)

dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)

path = dir_path + '/Patric/organisms/'
print(path)

import_files_obj = ImportFilesPatric(path, '.contigs.fasta','.xls')
dict_files = import_files_obj.getOrganismsFiles()



#Information for the bacterium





manageOrganismsContent(dict_files)


    #### Test old methods
for contig_id in list_contigs_ids_fasta:
    print(contig_id)
    list_proteins = xls_genbank_patric_obj.get_proteins_objects_by_contig_id(contig_id)
print('fini')