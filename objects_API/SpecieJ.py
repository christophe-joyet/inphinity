from marshmallow import Schema, fields, post_load

from rest_client.SpecieRest import SpecieAPI

class SpecieSchema(Schema):
    """
    This class map the json into the object Specie

    ..note:: see marshmallow API
    """
    id = fields.Int()
    designation = fields.Str()
    genus = fields.Int()
    strains = fields.List(fields.Url)

    @post_load
    def make_Specie(self, data):
        return SpecieJson(**data)


class SpecieJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, designation = '', genus = None, strains=[]):
        """
        Initialization of the class

        :param id: name of the function
        :param designation: name of the specie
        :param genus: id of the genus
        :param strains: vector of strains endpoints

        :type id: int
        :type designation: string 
        :type genus: int 
        :type strains: vector[string]

        """
        self.id = id
        self.designation = designation
        self.genus = genus
        self.strains = strains

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} designation {1} FK Genus {2} Nb strains {3}'.format(self.id, self.designation, self.genus, len(self.strains))

    def getAllAPI():

        """
        get all the species on the database

        :return: list of Specie
        :rtype: vector[SpecieJ]
        """
        list_specie = SpecieAPI().get_all()
        schema = SpecieSchema()
        results = schema.load(list_specie, many=True)
        return results

    def setSpecie(self):
        """
        set new genus

        :return: new specie completed with the id
        :rtype: SpecieJ
        """
        schema = SpecieSchema(only=['designation', 'genus'])
        jsonFam = schema.dump(self)
        resultsCreation = SpecieAPI().set_specie(jsonData = jsonFam.data)
        schema = SpecieSchema()
        results = schema.load(resultsCreation)
        return results[0]


    def getByID(id_specie:int):

        """
        get a specie given its id

        :param id_specie: id of the specie that it will be returned

        :type id_specie: int

        :return: a json of the specie 
        :rtype: SpecieJson
        """
        specie = SpecieAPI().get_by_id(id_specie)
        schema = SpecieSchema()
        results = schema.load(specie, many=False)
        return results[0]


