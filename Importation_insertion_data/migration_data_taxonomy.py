import csv
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI


from objects_API.FamilyJ import FamilyJson
from objects_API.GenusJ import GenusJson
from objects_API.SpecieJ import SpecieJson
from objects_API.StrainJ import StrainJson

from objects_new.Families_new import Family
from objects_new.Genus_new import Genus
from objects_new.Species_new import Specie
from objects_new.Strains_new import Strain



#===============================================
# Script used to migrate the data (only taxonomy)
# to the new API
#===============================================


def create_family_json(familyBD):
    """
    insert a family into a REST API

    :param familyBD: family that you want to insert

    :type familyBD: Family

    """

    famObjREST = FamilyJson(designation=familyBD.designation)
    familyREST = famObjREST.setFamily()

    listGenusBD = Genus.get_genus_by_family_id(familyBD.id_family)

    idFamilyJson = familyREST.id

    for genusObj in listGenusBD:
        create_genus_json(genusObj, idFamilyJson)
    

def create_genus_json(genusOBJ, idFamJson):
    """
    insert a genus into a REST API

    :param genusOBJ: genus that you want to insert
    :param idFamJson: FK of the genus family

    :type genusOBJ: Genus
    :type idFamJson: int

    """
    genusObjREST = GenusJson(designation=genusOBJ.designation, family=idFamJson)
    genusREST = genusObjREST.setGenus()

    listSpecieBD = Specie.get_all_species_by_genus_id(genusOBJ.id_genus)

    idGenusJson = genusREST.id

    for specieObj in listSpecieBD:
        create_specie_json(specieObj, idGenusJson)


def create_specie_json(specieOBJ, idGenus):
    """
    insert a specie into a REST API

    :param specieOBJ: specie that you want to insert
    :param idGenus: FK of the family specie

    :type specieOBJ: Genus
    :type idGenus: int

    """
    specieObjREST = SpecieJson(designation=specieOBJ.designation, genus=idGenus)
    specieREST = specieObjREST.setSpecie()

    listStrainsBD = Strain.get_strains_by_specie_id(specieOBJ.id_specie)

    idSpecieJson = specieREST.id

    for strainObj in listStrainsBD:
        create_strain_json(strainObj, idSpecieJson)


def create_strain_json(strainOBJ, idSpecie):
    """
    insert a strain into a REST API

    :param strainOBJ: strain that you want to insert
    :param idGenus: FK of the specie strain

    :type strainOBJ: strain
    :type idGenus: int

    """
    print(strainOBJ)
    strainObjREST = StrainJson(designation=strainOBJ.designation, specie=idSpecie)
    strainREST = strainObjREST.setStrain()
    createCSVCorrespondence(strainOBJ.id_strain, strainREST.id)
    

def createCSVCorrespondence(idOldstrain, idNewStrain):
    """
    Save the ids correspondence between the objects created 
    for the API and those in the database into a csv file

    :param idOldstrain: id of the strain in the database
    :param idNewStrain: id of the strain in the API

    :type idOldstrain: int
    :type idNewStrain: int

    """
    with open('correspondenceIDSStrains.csv','a', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows([[idOldstrain,idNewStrain]])


conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()
