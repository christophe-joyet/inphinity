from marshmallow import Schema, fields, post_load

from rest_client.CogInteractionPairRest import CogInteractionPairAPI

class CogInteractionPairSchema(Schema):
    """
    This class map the json into the object Cog interaction pair

    ..note:: see marshmallow API
    """
    id = fields.Int()
    cog_a = fields.Int()
    cog_b = fields.Int()

    @post_load
    def make_Cog(self, data):
        return CogInteractionPairJson(**data)

class CogInteractionPairJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self,  id = None, cog_a = -1, cog_b = -1):
        """
        Initialization of the class

        :param id: name of the function
        :param cog_a_id: id of the cog A
        :param cog_b_id: id of the cog B

        :type id: int
        :type cog_a_id: int 
        :type cog_b_id: int 

        """
        self.id = id
        self.cog_a = cog_a
        self.cog_b = cog_b

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} cog A {1} cog B {2}'.format(self.id, self.cog_a, self.cog_b)

    def getAllAPI():

        """
        get all the cogs interaction pair on the database

        :return: list of Cogs interaction pair
        :rtype: vector[CogInteractionPairJ]
        """
        list_cogs_interaction_pair = CogInteractionPairAPI().getAll()
        schema = CogInteractionPairSchema()
        results = schema.load(list_cogs_interaction_pair, many=True)
        return results

    def setCogInteractionPair(self):
        """
        set new cog interaction pair

        :return: new cog interaction pair completed with the id
        :rtype: CogJ
        """
        schema = CogInteractionPairSchema()

        json_cog = schema.dump(self)
        resultsCreation = CogInteractionPairAPI().setCog(jsonData = json_cog)
        schema = CogInteractionPairSchema()
        results = schema.load(resultsCreation)
        return results

    def getByID(id_cog_interaction_pair:int):

        """
        get a Cog interaction pair given its id

        :param id_cog_interaction_pair: id of the cog that it will be returned

        :type id_cog_interaction_pair: int

        :return: a json of the cog interaction pair
        :rtype: CogJson
        """
        cog_interaction_pair = CogInteractionPairAPI().getById(id_cog_interaction_pair)
        schema = CogInteractionPairSchema()
        results = schema.load(cog_interaction_pair, many=False)
        return results

    def getCogsInteractionPairByFilterParameter(dict_parameters:dict):
        """
        get a list of cogs given a filters by fields E.G: dict['level']=1
        return all cogs interaction pair with the level at 1

        :param dict_parameters: dictionary that contain the fields and vaules to filter

        :type dict_parameters: dictionary

        :return: a json of the cog interaction pair
        :rtype: list[CogInteractionPairJson]
        """
        url_parameters = ''
        for key_param in dict_parameters:
            url_parameters += key_param + '=' + str(dict_parameters[key_param]) + '&'

        url_parameters = url_parameters[:-1]
        list_cogs_interaction_pair = CogInteractionPairAPI().getCogsByParameters(url_parameters)
        schema = CogInteractionPairSchema()
        results = schema.load(list_cogs_interaction_pair, many=True)
        return results[0]