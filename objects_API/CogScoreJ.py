from marshmallow import Schema, fields, post_load

from rest_client.CogScoreRest import CogScoreAPI

class CogsScoreSchema(Schema):
    """
    This class map the json into the object Cog

    ..note:: see marshmallow API
    """
    id = fields.Int()
    score = fields.Float()
    cog_method_score = fields.Int()
    cog_interaction_source = fields.Int()

    @post_load
    def make_CogScore(self, data):
        return CogScoreJson(**data)

class CogScoreJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, score = 0.0, cog_method_score = -1, cog_interaction_source = -1):
        """
        Initialization of the class

        :param id: id of the cog score
        :param score: score of the cogs interaction
        :param cog_methode_score: id of the cog score
        :param cog_interaction_source: id of the cog interaction source

        :type id: int
        :type score: Float 
        :type cog_methode_score: int 
        :type cog_interaction_source: int 

        """
        self.id = id
        self.score = score
        self.cog_method_score = cog_method_score
        self.cog_interaction_source = cog_interaction_source


    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} score {1} cog method score id {2} cog interaction source id {3}'.format(self.id, self.score, self.cog_method_score, self.cog_interaction_source)

    def getAllAPI():

        """
        get all the cogs scores on the database

        :return: list of Cogs
        :rtype: vector[CogScoreJ]
        """
        list_cogs_score = CogScoreAPI().getAll()
        print(list_cogs_score)
        schema = CogsScoreSchema()
        results = schema.load(list_cogs_score, many=True)
        return results

    def setCogScore(self):
        """
        set new cog score

        :return: new cog score completed with the id
        :rtype: CogScoreJ
        """
        schema = CogsScoreSchema()

        json_cog = schema.dump(self)
        resultsCreation = CogScoreAPI().setCogScore(jsonData = json_cog)
        schema = CogsScoreSchema()
        results = schema.load(resultsCreation)
        return results

    def getByID(id_cog_score:int):

        """
        get a Cog  score given its id

        :param id_cog_score: id of the cog method score that it will be returned

        :type id_cog_score: int

        :return: a json of the cog  method score
        :rtype: CogMethodScoreJson
        """
        cogMethodScore = CogScoreAPI().getById(id_cog_score)
        schema = CogsScoreSchema()
        results = schema.load(cogMethodScore, many=False)
        return results


    def getCogsScoreByFilterParameter(dict_parameters:dict):
        """
        get a list of cogs score given a filters by fields E.G: dict['level']=1
        return all couples with the level at 1

        :param dict_parameters: dictionary that contain the fields and vaules to filter

        :type dict_parameters: dictionary

        :return: a json of the couple
        :rtype: list[CogScoreJson]
        """
        url_parameters = ''
        for key_param in dict_parameters:
            url_parameters += key_param + '=' + str(dict_parameters[key_param]) + '&'

        url_parameters = url_parameters[:-1]

        list_cogs_Score = CogScoreAPI().getCogsScoreByParameters(url_parameters)
        schema = CogsScoreSchema()
        results = schema.load(list_cogs_Score, many=True)
        return results[0]