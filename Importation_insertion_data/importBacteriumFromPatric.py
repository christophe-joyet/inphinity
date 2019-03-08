import glob
import os 


from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI
from files_treatment_new.xls_gen_bank_patric import XlsGenBankPatric
from files_treatment_new.fasta_contigs_patric import FastaContigsPatric

from objects_new.Contigs_new import Contig

from objects_API.ContigJ import ContigJson
from objects_API.StrainJ import StrainJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.GeneJ import GeneJson
from objects_API.ProteinJ import ProteinJson

from Patric.ImportFiles import ImportFilesPatric


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

def createAndInsertElements(contig_obj, id_bacterium, xls_genbank_patric_obj):
    #listProts = xls_genbank_patric_obj.get_proteins_information_in_excel()
    list_proteins = xls_genbank_patric_obj.get_proteins_objects_by_contig_id(contig_obj.head)
    list_proteins_ids = xls_genbank_patric_obj.get_proteins_ids_by_contig_id(contig_obj.head)

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

def manageOrganismsContent(dict_files):
    for key, value in dict_files.items():


        contig_path_file = value[0]
        excel_path_file = value[1]



        xls_genbank_patric_obj = XlsGenBankPatric(path_file= excel_path_file, sheet_name = 'Features in Pseudomonas aerugin')
        contig_fasta_file_patric_obj = FastaContigsPatric(path_file = contig_path_file)

        list_contigs_ids_fasta = contig_fasta_file_patric_obj.get_list_contigs_id()

        list_contigs_ids_xls = xls_genbank_patric_obj.get_contigs_id_sorted()



        list_diff = list(set(list_contigs_ids_fasta) - set(list_contigs_ids_xls))
        print(list_diff)
        #assert len(list_contigs_ids_fasta) == len(list_contigs_ids_xls)

        #Strain creation -------------------
        strain_designation = key.replace("Pseudomonas aeruginosa ","")
        fk_specie = 417
        #strain_obj = StrainJson.verifyStrainExistanceDesignationFkSpecie(strain_designation, fk_specie)
        #if strain_obj == None:
        strain = createStrain(strain_designation, fk_specie)
        id_strain = strain.id
        #else:
        #    id_strain = strain_obj.id
        #Bacterium creation --------------------
        acc_number = 'Greg_Patric_' + strain_designation
        source_data = 5
        person_responsible = 2

        id_bacterium = createBacterium(acc_number, person_responsible, source_data, id_strain)


        list_of_contigs = contig_fasta_file_patric_obj.create_contigs_from_file()

        for contig_old in list_of_contigs:
            createAndInsertElements(contig_old, id_bacterium, xls_genbank_patric_obj)
 
            

#Token connection
conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()
#End token connection


#end Test strain

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