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



#===============================================
# Script used to migrate the data (Bacteria)
# to the new API
#===============================================

def writeCSVStateInsertion(arrayContentStates):
    np.savetxt('stateInsertionBacteria.csv', arrayContentStates.astype(int), fmt='%i', delimiter=",")

def writeCSVProteinNumber(idOrganism, qtyProteins):
    with open('proteinsNumber.csv', 'a') as csvFile:
        writer = csv.writer(csvFile, delimiter = ',', lineterminator='\n')
        writer.writerow([idOrganism, qtyProteins])
    csvFile.close()

def writeCSVProteinContigNotEqual(idOrganism, qtyProteins, qtyProtsInContig):
    with open('contigProteinsNotEquals.csv', 'a') as csvFile:
        writer = csv.writer(csvFile, delimiter = ',', lineterminator='\n')
        writer.writerow([idOrganism, qtyProteins, qtyProtsInContig])
    csvFile.close()

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

#def getContigsbyid(idContig):
#    contigObj = Contig.get_contig_by_id(idContig)
#    return contigObj

#def getContigs(listProteins):
#    fkContig = 0
#    listIdContigs = []
#    qtyprots = 0
#    for proteinObj in listProteins:
#        fkContig = proteinObj.fk_id_contig
#        if fkContig != None and fkContig not in listIdContigs:
#            contigObj = getContigsbyid(fkContig)
#            listIdContigs.append(contigObj.id_contig)
#            qtyprots +=1
#        elif fkContig != None and fkContig in listIdContigs:
#            qtyprots +=1
#    return listIdContigs

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

def treatAccNumber(accNumber):
    """
    Used to return null if the acc does not exists

    :param accNumber: gi number

    :type accNumber: string

    :return: None if the acc number is non-existent
    :rtype None or string

    """
    if accNumber == 'NA':
        return None
    else:
        return accNumber

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

def verifyQtyProtContigs(qtyProteins, dictProteins, dictContigs):
    """
    Verify if all the proteins of the bacterium are mapped into the contigs

    :param qtyProteins: List of all the proteins of the bacterium
    :param dictProteins: List of all the contigs of the bacterium
    :param dictContigs: List of all the contigs of the bacterium

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

def mapProteinsContigs(listProteins, listContigs):
    """
    map the list of proteins with the contig into two dictionaries.

    :param listProteins: List of all the proteins of the bacterium
    :param listContigs: List of all the contigs of the bacterium

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

def createGeneProt(proteinOBJ, fkOrganism, fkContig = None):
    """
    insert the bacterium, the whole genome, the proteins, and the genes

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
    geneId = 0
    if fkContig != None:
        geneId = createGene(geneObjInsert, proteinOBJ.start_point_cnt, proteinOBJ.end_point_cnt, fkContig)
    else:
        geneId = createGene(geneObjInsert, None, None, None)
    proteinOBJInsert = Protein(id_accession = proteinOBJ.id_accession, sequence_prot = proteinOBJ.sequence_prot, designation = proteinOBJ.designation)

    if (proteinOBJ.sequence_prot != None):
        proteinId = createProtein(proteinOBJInsert, fkOrganism, geneId)
    
def createContigNew(contigObj, bacteriumId):
    head_cnt = ''
    if '>' not in contigObj.head:
        head_cnt = '>' + contigObj.head
    else: 
        head_cnt = contigObj.head

    contigObj = Contig(id_contig_db_outside = contigObj.id_contig_db_outside, head = head_cnt, sequence = contigObj.sequence)
    idContig = createContig(contigObj, bacteriumId)
    return idContig



def insertBacteriumProteinsWholeDNA(bacteriumObj, idStrainAPI, listProteins, wholeDNA):
    """
    insert the bacterium, the whole genome, the proteins, and the genes

    :param bacteriumObj: Organism object (old DB)
    :param idStrainAPI: id of the strain on the new DB (API)
    :param wholeDNA: Whole Genome object (old DB)
    :param listProteins: Array with all the proteins

    :type bacteriumObj: Organism_new
    :type idStrainAPI: int
    :type wholeDNA: WholeDNA_new
    :type listProteins: array[Protein]


    :return: matrix with all the data
    :rtype panda dataframe

    """
    bacteriumCreatedID = createBacterium(bacteriumObj, idStrainAPI)
    wholeDNACreatedID = createWholeDNA(wholeDNA, bacteriumCreatedID)

    for protein in listProteins:
        createGeneProt(protein, bacteriumCreatedID)

def testContigContainDNA(list_contig):
    """
    verify tha all the contigs contains a DNA sequence

    :param list_contig: contain a list of contigs that you want to verify

    :type list_contig: list(Contigs)


    :return: true or false according the validation
    :rtype bool

    """
    validation_contig = True
    for contig_obj in list_contig:
        if len(contig_obj.sequence) < 40:
            validation_contig = False
            break
    return validation_contig



def insertBacteriumProteinsWholeDNAContigs(bacteriumObj, idStrainAPI, wholeDNAObj, dictProts, dictConts):

    """
    insert the bacterium, the whole genome, the contigs, the proteins, and the genes

    :param bacteriumObj: Organism object (old DB)
    :param idStrainAPI: id of the strain on the new DB (API)
    :param wholeDNAObj: Whole Genome object (old DB)
    :param dictProts: Dict with the proteins Key: id of the contig - Value: list of proteins (old DB)
    :param dictConts: Dict with the contigs Key: id of the contig - Value: contig object (old DB)

    :type bacteriumObj: Organism_new
    :type idStrainAPI: int
    :type wholeDNAObj: WholeDNA_new
    :type dictProts: dictionary{int: list[Proteins]}
    :type dictConts: dictionary{int: Contig}


    :return: matrix with all the data
    :rtype panda dataframe

    """

    #Test contigs
    list_contigs = Contig.get_all_Contigs_by_organism_id(bacteriumObj.id_organism)
    contigs_validation = testContigContainDNA(list_contigs)
    if contigs_validation == True:
        bacteriumCreatedID = createBacterium(bacteriumObj, idStrainAPI)
        wholeDNACreatedID = createWholeDNA(wholeDNAObj, bacteriumCreatedID)

        for key,values in dictProts.items():
            contig = dictConts[key]
            listProteins = values
            contigID = createContigNew(contig, bacteriumCreatedID)
            for protein in listProteins:
                createGeneProt(protein, bacteriumCreatedID, contigID)

def insertBacteriumProteinsWholeDNAContigsNopositions(bacteriumObj, idStrainAPI, wholeDNAObj, dictProts, dictConts):

    """
    insert the bacterium, the whole genome, the contigs, the proteins, and the genes without protein s contig positions

    :param bacteriumObj: Organism object (old DB)
    :param idStrainAPI: id of the strain on the new DB (API)
    :param wholeDNAObj: Whole Genome object (old DB)
    :param dictProts: Dict with the proteins Key: id of the contig - Value: list of proteins (old DB)
    :param dictConts: Dict with the contigs Key: id of the contig - Value: contig object (old DB)

    :type bacteriumObj: Organism_new
    :type idStrainAPI: int
    :type wholeDNAObj: WholeDNA_new
    :type dictProts: dictionary{int: list[Proteins]}
    :type dictConts: dictionary{int: Contig}


    :return: matrix with all the data
    :rtype panda dataframe

    """

    #Test contigs
    list_contigs = Contig.get_all_Contigs_by_organism_id(bacteriumObj.id_organism)
    contigs_validation = testContigContainDNA(list_contigs)
    if contigs_validation == True:
        bacteriumCreatedID = createBacterium(bacteriumObj, idStrainAPI)
        wholeDNACreatedID = createWholeDNA(wholeDNAObj, bacteriumCreatedID)

        for contig in list_contigs:
            print('Hello')
            contigID = createContigNew(contig, bacteriumCreatedID)
        for key,proteinlist in dictProts.items():
            for protein in proteinlist:
                createProtein(protein, bacteriumCreatedID, None)
                print('Hello')

def checkBacteriumExistsByAcc(acc_number):
    print(acc_number)
    bacterium_existence = BacteriumJson.verifiyBacteriumExistanceByAcc(acc_number)
    return bacterium_existence

def load_get_bacterium(pathFile, arrayStates):
    """
    read the CSV with the organism id, strain id and their states. This method call the other for the insertion of the bacterium on the database

    :param pathFile: path of teh CSV file

    :type pathFile: String

    :return: matrix with all the data
    :rtype panda dataframe

    """
    wholeDNA = None
    dataframe = pd.read_csv(pathFile)
    for index, row in dataframe.iterrows():
        print(row['organism_type'])
        if row['organism_type'] == 1:
            row_verification = [[row['strain_db'], row['strain_api'],2]]

            id_strain_API = row['strain_api']
            if id_strain_API not in arrayStates[:, 1]:
                arrayStates = np.append(arrayStates, row_verification, axis=0)
                writeCSVStateInsertion(arrayStates)

                bacteriumList = Organism.get_organism_by_fk_strain(row['strain_db'])
                assert len(bacteriumList) == 1

                bacterium = bacteriumList[0]
                bacterium.acc_num = bacterium.acc_num + '_V2'

                #Test acc_value
                bacterium_existance = checkBacteriumExistsByAcc(bacterium.acc_num)
                #############
                if bacterium_existance is True:
                    print('This acc {0} already exists'.format(bacterium.acc_num))
                    continue



                accNumber = treatAccNumber(bacterium.acc_num)
                source_data = bacterium.fk_source_data
                person_responsible = bacterium.fk_source
                gi_number = treatGiNumber(bacterium.gi)
                wholeDNAObj = getWholeGenomById(bacterium.fk_whole_genome)
                strain = bacterium.fk_strain

                organism_code = row['strain_db']
                organism_codeAPI = row['strain_api']
                listProts = getProteinsByOrganism(bacterium.id_organism)
                listContigs = getContigsByOrganis(bacterium.id_organism)



                if len(listContigs)> 0 and len(listContigs[0].sequence) > 10:
                    dictProts, dictConts = mapProteinsContigs(listProts, listContigs)
                    canInsert, qtyProts = verifyQtyProtContigs(len(listProts),dictProts,dictConts)
                    if canInsert == True:
                        print('Hello)')
                        insertBacteriumProteinsWholeDNAContigs(bacterium, organism_codeAPI, wholeDNAObj, dictProts, dictConts)
                    elif qtyProts == 0:
                        insertBacteriumProteinsWholeDNAContigsNopositions(bacterium, organism_codeAPI, wholeDNAObj, dictProts, dictConts)
                    else:
                        writeCSVProteinContigNotEqual(bacterium.id_organism, len(listProts),qtyProts)
                else:
                    print("Hello")
                    insertBacteriumProteinsWholeDNA(bacterium, organism_codeAPI, listProts, wholeDNAObj)

                row_verification = [row['strain_db'], row['strain_api'],1]
                arrayStates[-1] = row_verification
                writeCSVStateInsertion(arrayStates)


def load_CSV_Correspondance(pathFile):
    """
    read the CSV with the organism id and strain id

    :param pathFile: path of teh CSV file

    :type pathFile: String

    :return: matrix with all the data
    :rtype panda dataframe

    """

    df = pd.read_csv(pathFile)
    return df

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

def create_sourceData_json(sourceDataDB):
    """
    insert a source Data into a REST API

    :param sourceDataDB: Source data that you want to insert

    :type sourceDataDB: SourceData

    :return: id of the source Data inserted
    :rtype int

    """
    sourceDataREST = SourceDataJson(designation=sourceDataDB.designation)
    sourceDataREST = sourceDataREST.setSourceData()
    return sourceDataREST.id

def get_organism_type_by_strain(fkStrain):
    """
    check the organism type and verify if any errors exists

    :param fkStrain: FK strain of the organism

    :type fkStrain: int

    :return: error code or type of the organism (1 bacterium - 2 phage)
    :rtype int

    :note: -1 - more than one organism with this strain
           -2 - No organism with this strain
           -3 - Other errors
           

    """
    listOrganism = Organism.get_organism_by_fk_strain(fkStrain)
    qtyOrganism = len(listOrganism)
    if qtyOrganism > 1:
        return -1
    elif qtyOrganism == 0:
        return -2
    elif qtyOrganism == 1:
        return listOrganism[0].fk_type
    else:
        return -3


def check_organisms(pandaDatafram):
    """
    check add fields to dataframe and confirm the validity of the organism on the database

    :param pandaDatafram: dataframe previously saved

    :type pandaDatafram: pandas Datafram

    :return: new dataframe with the corrections
    :rtype pandas dataframe

        
    """
    pandaDatafram['error details'] = '--'
    pandaDatafram['organism Type'] = '--'
    for index, row in pandaDatafram.iterrows():
        organism_code = get_organism_type_by_strain(row['strain_db'])
        if organism_code == 1 or organism_code == 2:
            pandaDatafram.loc[index,'organism Type'] = organism_code
        else:
            pandaDatafram.loc[index,'error details'] = organism_code
    return pandaDatafram

def creatPersonResponsible(sourcePersonDB):
    """
    insert a Person Responsible into a REST API

    :param sourcePersonDB: person responsible that you want to insert

    :type sourcePersonDB: Source

    :return: id of the Person Responsible inserted
    :rtype int
    """

    personResponsibleREST = PersonResponsibleJson(name=sourcePersonDB.designation)
    personResponsibleREST = personResponsibleREST.setPersonResponsible()
    return personResponsibleREST.id

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
    proteinObjJson = ProteinJson(id_db_online = proteinOBJ.id_accession, organism = fkOrganism, gene = fkGene, sequence_AA = proteinOBJ.sequence_prot, fasta_head = proteinOBJ.sequence_prot, description = proteinOBJ.designation, accession_num = proteinOBJ.id_accession, position_start = proteinOBJ.start_point, position_end = proteinOBJ.end_point)
    proteinObjJson = proteinObjJson.setProtein()
    return proteinObjJson.id

def createBacterium(bacteriumOBJ, fkStrain):
    """
    insert a Bacterium into a REST API

    :param bacteriumOBJ: Organism DBA object that you want to insert
    :param fkStrain: id of the strain

    :type proteinOBJ: Organism
    :type fkStrain: int

    :return: id of the bacterium inserted
    :rtype int
    """

    bacteriumObjJson = BacteriumJson(acc_number = bacteriumOBJ.acc_num, gi_number = bacteriumOBJ.gi, person_responsible = bacteriumOBJ.fk_source, source_data = bacteriumOBJ.fk_source_data, strain = fkStrain)
    bacteriumObjJson = bacteriumObjJson.setBacterium()
    return bacteriumObjJson.id

#def createBacteriophage(bacteriophageOBJ,  fkPersonResp, fkSourceData, fkBaltimore):

#    bacteriophageObjJson = BacteriophageJson(acc_number = bacteriophageOBJ.acc_num, gi_number = bacteriophageOBJ.gi, person_responsible = fkPersonResp, source_data = fkSourceData, baltimore_classification = fkBaltimore)
#    bacteriophageObjJson = bacteriophageObjJson.setBacteriophage()
#    return bacteriophageObjJson.id



conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()


#Test Source data#################

#sourcesList = SourceDataJson.getAllAPI()

#sourceDataObj = Source_data(designation='Patric - Platform Greg')
#print(sourceDataObj)

#sourceDataNew = create_sourceData_json(sourceDataObj)
#print(sourceDataNew)

#Test Person Responsible#################

#personResponsibleList = PersonResponsibleJson.getAllAPI()
#print(personResponsibleList)


#personResponsibleObj = Source(designation='Gregory Resch')
#print(personResponsibleObj)

#personRespoNew = creatPersonResponsible(personResponsibleObj)
#print(personRespoNew)


#Test Whole DNA#################

#listWholeDNA = WholeDNAJson.getAllAPI()
#print(listWholeDNA)



#wholeDNAObj = WholeDNA(head = '>44454wwww', head_id = 1545454, sequence = 'ACCTTGGG')
#wholeDNAInserted = createWholeDNA(wholeDNAObj, 14)
#print(wholeDNAInserted)

#Test Contig #################
#listContig = ContigJson.getAllAPI()
#print(listContig)

#contigObj = Contig(id_contig_db_outside = '1231qq2DqiogoDD', head = ">asdasd", sequence = "ACCCTTTGGGGAAA")
#contigInserted = createContig(contigObj, 14)
#print(contigInserted)

#Test Gene #################
#listGenes = GeneJson.getAllAPI()
#print(listGenes)

#geneObj = Gene(dna_head = ">124", 
#                 dna_sequence = "ASQWASAS", start_position = 0, 
#                 end_position = 3, FK_id_organism = 14, FK_id_protein = 1, gene_number = 5)
#geneInserted = createGene(geneObj)

#Test Protein ##########
#listProteins = ProteinJson.getAllAPI()
#print(listProteins)

#proteinOBJ = Protein(designation = 'qqqq', sequence_prot = 'AAABBBCDSDSDS', start_point = 1, end_point = 5, start_point_cnt = 56, end_point_cnt = 108, fk_id_contig = 5)

#proteinIdCreated = createProtein(proteinOBJ, 15, 1, 'wqeqwe')
#print(proteinIdCreated)

##Test Bacterium ##########
#listBacterium = BacteriumJson.getAllAPI()
#print(listBacterium)

#bacteriumObj = Organism(gi = 'qqqqqqQwqwDGI', acc_num = 'qQwqwqDwqwACC', fk_source = 10, fk_strain = 15, fk_type = 1, fk_whole_genome = 5, fk_source_data = 9)

#bacteriumCreated = createBacterium(bacteriumObj, 3,3,2)
#print(bacteriumCreated)


##Test Bacteriophage ##########
#listBacteriophage = BacteriophageJson.getAllAPI()
#print(listBacteriophage)

#bacteriophageObj = Organism(gi = 'BacteriophageGI', acc_num = 'qQwqttwqDwqwACC', fk_source = 10, fk_strain = 15, fk_type = 1, fk_whole_genome = 5, fk_source_data = 9)

#bacteriophageCreated = createBacteriophage(bacteriophageObj, 3,3,2)
#print(bacteriophageCreated)

#===============================================
# Update the datafram that contain the strains.
# Verify the vaility of the organisms and strains
#===============================================

#Load CSV
#pathFile = './correspondenceIDSStrains.csv'
#dfPanda = load_CSV_Correspondance(pathFile)
#print(dfPanda)

#dfPanda = check_organisms(dfPanda)
#dfPanda.to_csv('./correspondenceIDSStrains2.csv', sep=',')
#listOrganism = Organism.get_organism_by_fk_strain(4578)
#for element in listOrganism:
#    print(element)

#===============================================
# Load only the bacteria and start the impor-
# tation inside the new db (API rest)
#===============================================

dataFrame = pd.read_csv('./correspondenceIDSStrains4.csv')

path = './correspondenceIDSStrains4.csv'
pathInsertion = './stateInsertionBacteria.csv'

dataframState = load_CSV_Insertion(pathInsertion)
load_get_bacterium(path, dataframState)

print('end')
