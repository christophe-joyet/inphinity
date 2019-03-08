import pandas as pd
import csv
import numpy as np


from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

from objects_API.SourceDataJ import SourceDataJson
from objects_API.PersonResponsibleJ import PersonResponsibleJson
from objects_API.WholeDNAJ import WholeDNAJson
from objects_API.ContigJ import ContigJson
from objects_API.GeneJ import GeneJson
from objects_API.ProteinJ import ProteinJson
from objects_API.BacteriumJ import BacteriumJson
from objects_API.BacteriophageJ import BacteriophageJson

from objects_new.Sources_data_new import Source_data
from objects_new.Sources_new import Source
from objects_new.Contigs_new import Contig
from objects_new.Gene_new import Gene
from objects_new.Proteins_new import Protein
from objects_new.Organisms_new import Organism
from objects_new.WholeDNA_new import WholeDNA
from objects_new.Strains_new import Strain


#===============================================
# Script used to migrate the data (Bacteriophage)
# to the new API
#===============================================

def writeCSVStateInsertion(arrayContentStates):
    np.savetxt('stateInsertionBacteriophage_V2.csv', arrayContentStates.astype(int), fmt='%i', delimiter=",")

def writeCSVProteinNumber(idOrganism, qtyProteins):
    with open('proteinsNumberBacteriophage_V2.csv', 'a') as csvFile:
        writer = csv.writer(csvFile, delimiter = ',', lineterminator='\n')
        writer.writerow([idOrganism, qtyProteins])
    csvFile.close()

def writeCSVProteinContigNotEqual(idOrganism, qtyProteins, qtyProtsInContig):
    with open('contigProteinsNotEqualsBacteriophage_V2.csv', 'a') as csvFile:
        writer = csv.writer(csvFile, delimiter = ',', lineterminator='\n')
        writer.writerow([idOrganism, qtyProteins, qtyProtsInContig])
    csvFile.close()

def treatAccNumber(accNumber, id_strain_bd):
    """
    Used to return null if the acc does not exists

    :param accNumber: gi number

    :type accNumber: string

    :return: None if the acc number is non-existent
    :rtype None or string

    """
    if accNumber == 'NA':
        return 'remove_' + str(id_strain_bd)
    else:
        return accNumber


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
    return contigObjJson.id

def treatGiNumber(giNumber):
    """
    Used to return null if the gi does not exists

    :param giNumber: gi number

    :type giNumber: string

    :return: None if the gi number is non-existent
    :rtype None or string

    """
    if giNumber == 'NA':
        return None
    else:
        return giNumber

def mapProteinsContigs(listProteins, listContigs):
    """
    map the list of proteins with the contig into two dictionaries.

    :param listProteins: List of all the proteins of the bacteriophage
    :param listContigs: List of all the contigs of the bacteriophage

    :type listProteins: List[Protein]
    :type listProteins: List[Contig]

    :return: two dictionaries with the proteins and contigs mapped
    :rtype dictProteins: dictionary{int: list[Proteins]}
    :rtype dictContigs: dictionary{int: Contig}

    """
    dictProteins = {}
    dictContigs = {}
    for element in listProteins:
        if element.fk_id_contig in dictProteins.keys():
            dictProteins[element.fk_id_contig].append(element)
        else:
            dictProteins[element.fk_id_contig] = [element]

    for element in listContigs:
        dictContigs[element.id_contig] = element
    
    return dictProteins, dictContigs

def createWholeDNA(wholeDNAObj, organismID):
    """
    insert a WholeDNA into a REST API

    :param wholeDNAObj: WholeDNA DBA object that you want to insert
    :param organismID: ID of the organism which has this wholeDNA

    :type wholeDNAObj: WholeDNA
    :type organismID: int

    :return: id of the WholeDNA inserted
    :rtype int
    """
    wholeDNAObjJson = WholeDNAJson(id_db_online = wholeDNAObj.head_id, sequence_DNA= wholeDNAObj.sequence, fasta_head = wholeDNAObj.head, organism=organismID)
    wholeDNAObjJson = wholeDNAObjJson.setWholeDNA()
    return wholeDNAObjJson.id

def verifyQtyProtContigs(qtyProteins, dictProteins, dictContigs):
    """
    Verify if all the proteins of the bacteriophage are mapped into the contigs

    :param qtyProteins: List of all the proteins of the bacteriophage
    :param dictProteins: List of all the contigs of the bacteriophage
    :param dictContigs: List of all the contigs of the bacteriophage

    :type qtyProteins: List[Protein]
    :type dictProts: dictionary{int: list[Proteins]}
    :type dictConts: dictionary{int: Contig}

    :return: True if the number of the proteins are exact of false in case of contrary and qty of counted proteins
    :rtype bool and int

    """
    qtyCountProt = 0
    for key, value in dictContigs.items():
        if key in dictProteins.keys():
            qtyCountProt += len(dictProteins[key])

    if qtyCountProt == qtyProteins:
        return True, qtyCountProt
    else:
        return False, qtyCountProt

def createBacteriophage(bacteriophageOBJ, strainDesignation):
    """
    insert a Bacteriophage into a REST API

    :param bacteriophageOBJ: Organism DBA object that you want to insert
    :param strainDesignation: designation of the strain

    :type proteinOBJ: Organism
    :type strainDesignation: String

    :return: id of the bacteriophage inserted
    :rtype int
    """


    strainDesignation = strainDesignation + '_' + bacteriophageOBJ.acc_num
    bacteriophageObjJson = BacteriophageJson(acc_number = bacteriophageOBJ.acc_num, gi_number = bacteriophageOBJ.gi, person_responsible = bacteriophageOBJ.fk_source, source_data = bacteriophageOBJ.fk_source_data, designation = strainDesignation)
    bacteriophageObjJson = bacteriophageObjJson.setBacteriophage()
    return bacteriophageObjJson.id

def createGene(geneObj, start_contig, end_contig, fk_contig):
    """
    insert a Gene into a REST API

    :param geneObj: Gene DBA object that you want to insert
    :param start_contig: start position of the gene in the contig
    :param end_contig: end position of the gene in the contig
    :param fk_contig: id of the contig

    :type geneObj: WholeDNA
    :type start_contig: int - can be None
    :type end_contig: int - can be None
    :type fk_contig: int - can be None

    :return: id of the gene inserted
    :rtype int
    """
    geneObjJson = GeneJson(sequence_DNA = geneObj.dna_sequence, fasta_head = geneObj.dna_head, position_start = geneObj.start_position, position_end = geneObj.end_position, organism = geneObj.FK_id_organism, number_of_seq = geneObj.gene_number, position_start_contig = start_contig, position_end_contig = end_contig, contig = fk_contig)
    geneObjJson = geneObjJson.setGene()
    return geneObjJson.id


def createProtein(proteinOBJ, fkOrganism, fkGene):
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
    #proteinObjJson = ProteinJson(id_db_online = proteinOBJ.id_accession, organism = fkOrganism, gene = fkGene, sequence_AA = proteinOBJ.sequence_prot, fasta_head = proteinOBJ.sequence_prot, description = proteinOBJ.designation, accession_num = proteinOBJ.id_accession, position_start = proteinOBJ.start_point, position_end = proteinOBJ.end_point)

    proteinObjJson = ProteinJson(id_db_online = proteinOBJ.id_accession, organism = fkOrganism, gene = fkGene, sequence_AA = proteinOBJ.sequence_prot, fasta_head = proteinOBJ.sequence_prot, description = 'No description', accession_num = 'No acc', position_start = proteinOBJ.start_point, position_end = proteinOBJ.end_point)

    proteinObjJson = proteinObjJson.setProtein()
    return proteinObjJson.id

def createGeneProt(proteinOBJ, fkOrganism, fkContig = None):
    """
    insert the bacteriophage, the whole genome, the proteins, and the genes

    :param proteinOBJ: Protein object (old DB)
    :param fkOrganism: id of the organism created
    :param fkContig: id of the contig created in case of exists

    :type proteinOBJ: Protein_new
    :type fkOrganism: int
    :type fkContig: int


    """

    dna_head_gene = ''
    if '>' not in proteinOBJ.designation:
        dna_head_gene = '>' + proteinOBJ.designation
    else: 
        dna_head_gene = proteinOBJ.designation

    geneObjInsert = Gene(dna_head = dna_head_gene, 
                     dna_sequence = proteinOBJ.sequence_dna, FK_id_organism = fkOrganism, start_position = proteinOBJ.start_point, end_position = proteinOBJ.end_point)    
    geneId = None
    if fkContig != None:
        geneId = createGene(geneObjInsert, proteinOBJ.start_point_cnt, proteinOBJ.end_point_cnt, fkContig)
    elif proteinOBJ.sequence_dna is not None:
        geneId = createGene(geneObjInsert, None, None, None)
    proteinOBJInsert = Protein(id_accession = proteinOBJ.id_accession, sequence_prot = proteinOBJ.sequence_prot, designation = proteinOBJ.designation)

    if (proteinOBJ.sequence_prot != None):
        proteinId = createProtein(proteinOBJInsert, fkOrganism, geneId)




def createContigNew(contigObj, bacteriophageId):
    head_cnt = ''
    if '>' not in contigObj.head:
        head_cnt = '>' + contigObj.head
    else: 
        head_cnt = contigObj.head

    contigObj = Contig(id_contig_db_outside = contigObj.id_contig_db_outside, head = head_cnt, sequence = contigObj.sequence)
    idContig = createContig(contigObj, bacteriophageId)
    return idContig

def getWholeGenomById(id_whole_genome):
    """
    get a whole genome according to it id

    :param id_whole_genome: Id of the whole genome

    :type id_whole_genome: int

    :return: whole genome Object
    :rtype WholeDNA

    """
    wholeDNAObj = WholeDNA.get_whole_dna_by_id(id_whole_genome)
    return wholeDNAObj

def insertBacteriophageProteinsWholeDNAContigs(bacteriophageObj, idStrainDB, wholeDNAObj, dictProts, dictConts):

    """
    insert the bacteriophage, the whole genome, the contigs, the proteins, and the genes

    :param bacteriophageObj: Organism object (old DB)
    :param idStrainDB: id of the strain on the DB (oldDB)
    :param wholeDNAObj: Whole Genome object (old DB)
    :param dictProts: Dict with the proteins Key: id of the contig - Value: list of proteins (old DB)
    :param dictConts: Dict with the contigs Key: id of the contig - Value: contig object (old DB)

    :type bacteriophageObj: Organism_new
    :type idStrainDB: int
    :type wholeDNAObj: WholeDNA_new
    :type dictProts: dictionary{int: list[Proteins]}
    :type dictConts: dictionary{int: Contig}


    :return: matrix with all the data
    :rtype panda dataframe

    """
    strainObj = Strain.get_strain_by_id(idStrainDB)
    bacteriophageCreatedID = createBacteriophage(bacteriophageObj, strainObj.designation)
    wholeDNACreatedID = createWholeDNA(wholeDNAObj, bacteriophageCreatedID)

    for key,values in dictProts.items():
        contig = dictConts[key]
        listProteins = values
        contigID = createContigNew(contig, bacteriophageCreatedID)
        for protein in listProteins:
            createGeneProt(protein, bacteriophageCreatedID, contigID)

def getProteinsByOrganism(idOrganism):
    """
    get all the proteins of a given Organism

    :param idOrganism: Id of the organism

    :type idOrganism: int

    :return: list of proteins
    :rtype list[Protein]

    """
    listProteins = []
    qtyprots = 0
    listProteins = Protein.get_all_Proteins_by_organism_id(idOrganism)
    assert len(listProteins) > 1
    writeCSVProteinNumber(idOrganism, len(listProteins))
    return listProteins

def getContigsByOrganis(idOrganism):
    """
    get all the contigs of a given Organism

    :param idOrganism: Id of the organism

    :type idOrganism: int

    :return: list of contigs
    :rtype list[Contig]

    """
    listContigs = []
    qtyContigs = 0
    listContigs = Contig.get_all_Contigs_by_organism_id(idOrganism)
    return listContigs

def insertBacteriophageProteinsWholeDNA(bacteriophageObj, idStrainDB, listProteins, wholeDNA):
    """
    insert the bacteriophage, the whole genome, the proteins, and the genes

    :param bacteriophageObj: Organism object (old DB)
    :param idStrainDB: id of the strain on the DB (old DB)
    :param wholeDNA: Whole Genome object (old DB)
    :param listProteins: Array with all the proteins

    :type bacteriophageObj: Organism_new
    :type idStrainDB: int
    :type wholeDNA: WholeDNA_new
    :type listProteins: array[Protein]


    :return: matrix with all the data
    :rtype panda dataframe

    """
    strainObj = Strain.get_strain_by_id(idStrainDB)
    bacteriophageCreatedID = createBacteriophage(bacteriophageObj, strainObj.designation)
    wholeDNACreatedID = createWholeDNA(wholeDNA, bacteriophageCreatedID)
    for protein in listProteins:
        createGeneProt(protein, bacteriophageCreatedID)

def load_get_bacteriophage(pathFile, arrayStates):
    """
    read the CSV with the organism id, strain id and their states. This method call the other for the insertion of the bacteriophage on the database

    :param pathFile: path of teh CSV file

    :type pathFile: String

    :return: matrix with all the data
    :rtype panda dataframe

    """

    arrayProblems = [11548, 11549]
    wholeDNA = None
    dataframe = pd.read_csv(pathFile)

    for index, row in dataframe.iterrows():
        if row['organism_type'] == 4:
            row_verification = [[row['strain_db'], row['strain_api'],2]]
            
            id_strain_API = row['strain_api']
            if id_strain_API not in arrayStates[:, 1] and id_strain_API not in arrayProblems:
                print(row['strain_db'])
                arrayStates = np.append(arrayStates, row_verification, axis=0)
                writeCSVStateInsertion(arrayStates)

                bacteriophageList = Organism.get_organism_by_fk_strain(row['strain_db'])
                assert len(bacteriophageList) == 1



                #if bacteriophageList[0].qty_proteins > bacteriophageList[1].qty_proteins:
                #    bateriophage_take = 0
                #else:
                #    bateriophage_take = 1
                bateriophage_take = 0
                bacteriophage = bacteriophageList[bateriophage_take]

                organism_code = row['strain_db']
                organism_codeAPI = row['strain_api']

                accNumber = treatAccNumber(bacteriophage.acc_num, organism_code)
                bacteriophage.acc_num = accNumber
                source_data = bacteriophage.fk_source_data
                person_responsible = bacteriophage.fk_source
                gi_number = treatGiNumber(bacteriophage.gi)
                wholeDNAObj = getWholeGenomById(bacteriophage.fk_whole_genome)
                strain = bacteriophage.fk_strain
                

                print(strain)



                listProts = getProteinsByOrganism(bacteriophage.id_organism)
                listContigs = getContigsByOrganis(bacteriophage.id_organism)

                if len(listContigs)> 0:
                    dictProts, dictConts = mapProteinsContigs(listProts, listContigs)
                    canInsert, qtyProts = verifyQtyProtContigs(len(listProts),dictProts,dictConts)
                    if canInsert == True:
                        print('Hello Contigs')
                        insertBacteriophageProteinsWholeDNAContigs(bacteriophage, organism_code, wholeDNAObj, dictProts, dictConts)
                    else:
                        writeCSVProteinContigNotEqual(bacteriophage.id_organism, len(listProts),qtyProts)
                        insertBacteriophageProteinsWholeDNA(bacteriophage, organism_code, listProts, wholeDNAObj)
                else:
                    print("Hello")
                    insertBacteriophageProteinsWholeDNA(bacteriophage, organism_code, listProts, wholeDNAObj)

                row_verification = [row['strain_db'], row['strain_api'],1]
                arrayStates[-1] = row_verification
                writeCSVStateInsertion(arrayStates)


def load_CSV_Insertion(pathFile):
    """
    read the CSV with the insertion data:
     - id strain API; id strain DB; id organism; state

    the states can be:
     1 - inserted
     2 - in course

    :param pathFile: path of the CSV file

    :type pathFile: String

    :return: matrix with all the data
    :rtype panda dataframe

    """

    contentCSV = np.genfromtxt(pathFile, delimiter=',', dtype=np.int32)
    return contentCSV


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

dataFrame = pd.read_csv('./correspondenceIDSStrains5.csv')


path = './correspondenceIDSStrains5.csv'
pathInsertion = './stateInsertionBacteriophage_V2.csv'


dataframState = load_CSV_Insertion(pathInsertion)


load_get_bacteriophage(path, dataframState)

