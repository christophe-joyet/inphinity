from marshmallow import Schema, fields, post_load

from rest_client.CogMethodScoreRest import CogMethodScoreAPI

class CogsMethodScoreSchema(Schema):
    """
    This class map the json into the object Cog

    ..note:: see marshmallow API
    """
    id = fields.Int()
    designation = fields.Str()

    @post_load
    def make_CogMethodScore(self, data):
        return CogMethodScoreJson(**data)


class CogMethodScoreJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, designation = ''):
        """
        Initialization of the class

        :param id: name of the function
        :param designation: name of the cog method score

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
        :rtype: vector[CogMethodScoreJ]
        """
        list_cogs_method_score = CogMethodScoreAPI().getAll()
        schema = CogsMethodScoreSchema()
        results = schema.load(list_cogs_method_score, many=True)
        return results

    def setCogMethodScore(self):
        """
        set new cog method score

        :return: new cog method score completed with the id
        :rtype: CogMethodScoreJ
        """
        schema = CogsMethodScoreSchema()

        json_cog = schema.dump(self)
        resultsCreation = CogMethodScoreAPI().setCogMethodScore(jsonData = json_cog)
        schema = CogsMethodScoreSchema()
        results = schema.load(resultsCreation)
        return results

    def getByID(id_cog:int):

        """
        get a Cog method score given its id

        :param id_cog_method_score: id of the cog method score that it will be returned

        :type id_cog_method_score: int

        :return: a json of the cog  method score
        :rtype: CogMethodScoreJson
        """
        cogMethodScore = CogMethodScoreAPI().getById(id_cog)
        schema = CogsMethodScoreSchema()
        results = schema.load(cogMethodScore, many=False)
        return results

    def getCogsMethodScoreByFilterParameter(dict_parameters:dict):
        """
        get a list of cogs method score given a filters by fields E.G: dict['level']=1
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

        list_cogs_method_Score = CogMethodScoreAPI().getCogsMethodScoreByParameters(url_parameters)
        schema = CogsMethodScoreSchema()
        results = schema.load(list_cogs_method_Score, many=True)
        return results[0]