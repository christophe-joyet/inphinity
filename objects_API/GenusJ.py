from marshmallow import Schema, fields, post_load

from rest_client.GenusRest import GenusAPI

class GenusSchema(Schema):
    """
    This class map the json into the object Genus

    ..note:: see marshmallow API
    """
    id = fields.Int()
    designation = fields.Str()
    family = fields.Int()
    species = fields.List(fields.Url)

    @post_load
    def make_Genus(self, data):
        return GenusJson(**data)

class GenusJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, designation = '', family = None, species=[]):
        """
        Initialization of the class

        :param id: name of the function
        :param designation: name of the genus
        :param family: id of the family
        :param species: vector of species endpoints

        :type id: int
        :type designation: string 
        :type family: int 
        :type species: vector[string]

        """
        self.id = id
        self.designation = designation
        self.family = family
        self.species = species

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} designation {1} FK family {2} Nb species {3}'.format(self.id, self.designation, self.family, len(self.species))

    def getAllAPI():

        """
        get all the genuses on the database

        :return: list of Genus
        :rtype: vector[GenusJ]
        """
        list_genus = GenusAPI().get_all()
        schema = GenusSchema()
        results = schema.load(list_genus, many=True)
        return results

    def setGenus(self):
        """
        set new genus

        :return: new genus completed with the id
        :rtype: GenusJ
        """
        schema = GenusSchema(only=['designation', 'family'])
        jsonFam = schema.dump(self)
        resultsCreation = GenusAPI().set_genus(jsonData = jsonFam.data)
        schema = GenusSchema()
        results = schema.load(resultsCreation)
        return results

    def getByID(id_genus:int):

        """
        get a Genus given its id

        :param id_genus: id of the genus that it will be returned

        :type id_genus: int

        :return: a json of the genus 
        :rtype: GenusJson
        """
        genus = GenusAPI().get_by_id(id_genus)
        schema = GenusSchema()
        results = schema.load(genus, many=False)
        return results[0]