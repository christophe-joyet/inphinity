from marshmallow import Schema, fields, post_load

from rest_client.FamilyRest import FamilyAPI


class FamilySchema(Schema):
    """
    This class map the json into the object Family

    ..note:: see marshmallow API
    """
    id = fields.Int()
    designation = fields.Str()
    genuses = fields.List(fields.Url)

    @post_load
    def make_Family(self, data):
        return FamilyJson(**data)

class FamilyJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, designation = '', genuses=[]):
        """
        Initialization of the class

        :param id: name of the function
        :param designation: name of the family
        :param genuses: vector of genus endpoints

        :type id: int
        :type designation: string 
        :type genuses: vector[string]

        """
        self.id = id
        self.designation = designation
        self.genuses = genuses

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} designation {1}'.format(self.id, self.designation)



    def getAllAPI():

        """
        get all the families on the database

        :return: list of Families
        :rtype: vector[FamilyJ]
        """
        list_family = FamilyAPI().get_all()
        schema = FamilySchema()
        results = schema.load(list_family, many=True)
        return results

    def setFamily(self):
        """
        set new family

        :return: new family completed with the id
        :rtype: FamilyJ
        """
        schema = FamilySchema(only=['designation'])
        jsonFam = schema.dump(self)
        resultsCreation = FamilyAPI().set_family(jsonData = jsonFam.data)
        schema = FamilySchema()
        results = schema.load(resultsCreation)
        return results[0]

    def getByID(id_family:int):

        """
        get a Family given its id

        :param id_family: id of the Family that it will be returned

        :type id_family: int

        :return: a json of the family 
        :rtype: FamilyJson
        """
        family = FamilyAPI().get_by_id(id_family)
        schema = FamilySchema()
        results = schema.load(family, many=False)
        return results[0]