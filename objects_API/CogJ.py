from marshmallow import Schema, fields, post_load

from rest_client.CogRest import CogAPI

class CogsSchema(Schema):
    """
    This class map the json into the object Cog

    ..note:: see marshmallow API
    """
    id = fields.Int()
    designation = fields.Str()

    @post_load
    def make_Cog(self, data):
        return CogJson(**data)

class CogJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, designation = ''):
        """
        Initialization of the class

        :param id: name of the function
        :param designation: name of the cog

        :type id: int
        :type designation: string 

        """
        self.id = id
        self.designation = designation

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} designation {1}'.format(self.id, self.designation)

    def getAllAPI():

        """
        get all the cogs on the database

        :return: list of Cogs
        :rtype: vector[CogJ]
        """
        list_cogs = CogAPI().getAll()
        schema = CogsSchema()
        results = schema.load(list_cogs, many=True)
        return results

    def setCog(self):
        """
        set new cog

        :return: new cog completed with the id
        :rtype: CogJ
        """
        schema = CogsSchema()

        json_cog = schema.dump(self)
        resultsCreation = CogAPI().setCog(jsonData = json_cog)
        schema = CogsSchema()
        results = schema.load(resultsCreation)
        return results


    def getByID(id_cog:int):

        """
        get a Cog given its id

        :param id_cog: id of the cog that it will be returned

        :type id_cog: int

        :return: a json of the cog 
        :rtype: CogJson
        """
        cog = CogAPI().getById(id_cog)
        schema = CogsSchema()
        results = schema.load(cog, many=False)
        return results

    def getCogsByFilterParameter(dict_parameters:dict):
        """
        get a list of cogs given a filters by fields E.G: dict['level']=1
        return all couples with the level at 1

        :param dict_parameters: dictionary that contain the fields and vaules to filter

        :type dict_parameters: dictionary

        :return: a json of the couple
        :rtype: list[CogJson]
        """
        url_parameters = ''
        for key_param in dict_parameters:
            url_parameters += key_param + '=' + str(dict_parameters[key_param]) + '&'

        url_parameters = url_parameters[:-1]
        list_cogs = CogAPI().getCogsByParameters(url_parameters)
        schema = CogsSchema()
        results = schema.load(list_cogs, many=True)
        return results[0]