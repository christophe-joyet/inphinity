from marshmallow import Schema, fields, post_load
import datetime

from rest_client.CogInteractionSourceRest import CogInteractionSourceAPI


class CogInteractionSourceSchema(Schema):
    """
    This class map the json into the object Domain Interaction Source

    ..note:: see marshmallow API
    """
    id = fields.Int()
    date_creation = fields.Date()
    cogs_interaction = fields.Int()
    information_source = fields.Int()

    @post_load
    def make_CogInteractionSource(self, data):
        return CogInteractionSourceJson(**data)

class CogInteractionSourceJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self,  date_creation:datetime, cogs_interaction:int, information_source:int, id = None):
        """
        Initialization of the class

        :param id: name of the function
        :param date_creation: cog ssocre interaction source date
        :param cogs_interaction: id of the cog interaction
        :param information_source: id of the information source

        :type id: int
        :type date_creation: date 
        :type cogs_interaction: int 
        :type information_source: int 

        """
        self.id = id
        self.date_creation = date_creation
        self.cogs_interaction = cogs_interaction
        self.information_source = information_source

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} creation date {1} cog interaction id {2} source id {3}'.format(self.id, self.date_creation, self.cogs_interaction, self.information_source)

    def getAllAPI():

        """
        get all the cog  interaction source on the database

        :return: list of cog  interaction source
        :rtype: vector[CogInteractionSourceJson]
        """
        list_cog_interaction_source = CogInteractionSourceAPI().getAll()

        schema = CogInteractionSourceSchema()
        results = schema.load(list_cog_interaction_source, many=True)
        return results

    def setCogInteractionSource(self):
        """
        set new cog  interaction source

        :return: new cog  interaction source completed with the id
        :rtype: CogInteractionSourceJ
        """
        schema = CogInteractionSourceSchema()

        json_cog = schema.dump(self)
        resultsCreation = CogInteractionSourceAPI().setCogInterSource(jsonData = json_cog)
        schema = CogInteractionSourceSchema()
        results = schema.load(resultsCreation)
        return results

    def getByID(id_cog__inter_source:int):

        """
        get a cog  interaction source given its id

        :param id_cog__inter_source: id of the cog  interaction source that it will be returned

        :type id_cog__inter_source: int

        :return: a json of the cog  interaction source
        :rtype: CogMethodJson
        """
        cogMethod = CogInteractionSourceAPI().getById(id_cog__inter_source)
        schema = CogInteractionSourceSchema()
        results = schema.load(cogMethod, many=False)
        return results

    def getCogsInteractSourceByFilterParameter(dict_parameters:dict):
        """
        get a list of cog  interaction source given a filters by fields E.G: dict['level']=1
        return all couples with the level at 1

        :param dict_parameters: dictionary that contain the fields and vaules to filter

        :type dict_parameters: dictionary

        :return: a json of the  cog  interaction source
        :rtype: list[CogMethodJson]
        """
        url_parameters = ''
        for key_param in dict_parameters:
            url_parameters += key_param + '=' + str(dict_parameters[key_param]) + '&'

        url_parameters = url_parameters[:-1]

        list_cogs_interact_source = CogInteractionSourceAPI().getCogsInteractionSourceByParameters(url_parameters)
        schema = CogInteractionSourceSchema()
        results = schema.load(list_cogs_interact_source, many=True)
        return results

