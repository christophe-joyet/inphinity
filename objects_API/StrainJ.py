from marshmallow import Schema, fields, post_load
import urllib.request

from rest_client.StrainRest import StrainAPI

class StrainSchema(Schema):
    """
    This class map the json into the object Strain

    ..note:: see marshmallow API
    """
    id = fields.Int()
    designation = fields.Str()
    specie = fields.Int()
    bacteria = fields.List(fields.Url)

    @post_load
    def make_Strain(self, data):
        return StrainJson(**data)

class StrainJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, designation = '', specie = None, bacteria=[]):
        """
        Initialization of the class

        :param id: name of the function
        :param designation: name of the strain
        :param specie: id of the specie
        :param bacteria: vector of species endpoints

        :type id: int
        :type designation: string 
        :type specie: int 
        :type bacteria: vector[string]

        """
        self.id = id
        self.designation = designation
        self.specie = specie
        self.bacteria = bacteria

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} designation {1} FK specie {2} Nb bacteria {3}'.format(self.id, self.designation, self.specie, len(self.bacteria))

    def getAllAPI():

        """
        get all the strains on the database

        :return: list of Strain
        :rtype: vector[StrainJ]
        """
        list_strain = StrainAPI().get_all()
        schema = StrainSchema()
        results = schema.load(list_strain, many=True)
        return results

    def setStrain(self):
        """
        set new strain

        :return: new strain completed with the id
        :rtype: StrainJ
        """
        schema = StrainSchema(only=['designation', 'specie'])
        jsonFam = schema.dump(self)
        resultsCreation = StrainAPI().set_strain(jsonData = jsonFam.data)
        schema = StrainSchema()
        results = schema.load(resultsCreation)
        return results[0]

    def getByID(id_strain:int):

        """
        get a strain given its id

        :param id_strain: id of the strain that it will be returned

        :type id_strain: int

        :return: a json of the strain 
        :rtype: StrainJson
        """
        strain = StrainAPI().get_by_id(id_strain)
        schema = StrainSchema()
        results = schema.load(strain, many=False)
        return results

    def verifyStrainExistanceDesignationFkSpecie(designation:str, fk_specie:int):
        """
        Verify if a strain already exists given its designation and fk_specie

        :param designation: designation of the strain
        :param fk_specie: fk of the specie

        :type designation: string
        :type fk_specie: integer

        :return: None or a strain if exists
        :rtype: StrainJ
        """

        designation = urllib.parse.quote(designation)
        strain_obj = StrainAPI().getByDesignationFkSpecie(designation, fk_specie)
        print(strain_obj['strain_exists'])
        if strain_obj['strain_exists'] == True:
            schema = StrainSchema(only=['designation', 'specie','id'])
            results = schema.load(strain_obj, many=False)
            return results[0]
        else:
            return None