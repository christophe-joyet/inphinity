from marshmallow import Schema, fields, post_load

from rest_client.CogSourceInformationRest import CogSourceInformationAPI

class CogSourceInformationSchema(Schema):
    """
    This class map the json into the object Cog Source Information

    ..note:: see marshmallow API
    """
    id = fields.Int()
    designation = fields.Str()

    @post_load
    def make_CogSourceInformation(self, data):
        return CogSourceInformationJson(**data)

class CogSourceInformationJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, designation = ''):
        """
        Initialization of the class

        :param id: name of the function
        :param designation: name of the cog source information

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
        get all the cogs source information on the database

        :return: list of Cogs soruce information
        :rtype: vector[CogSOurceINformationJ]
        """
        list_cogs_source_information = CogSourceInformationAPI().getAll()
        schema = CogSourceInformationSchema()
        results = schema.load(list_cogs_source_information, many=True)
        return results

    def setCogSourceInformation(self):
        """
        set new cog source information

        :return: new cog source information completed with the id
        :rtype: CogSourceInformationJ
        """
        schema = CogSourceInformationSchema()

        json_cog = schema.dump(self)
        resultsCreation = CogSourceInformationAPI().setCogSourceInformation(jsonData = json_cog)
        schema = CogSourceInformationSchema()
        results = schema.load(resultsCreation)
        return results

    def getByID(id_cog_source_information:int):

        """
        get a Cog source information given its id

        :param id_cog_source_information: id of the cog source information that it will be returned

        :type id_cog_source_information: int

        :return: a json of the cog source information
        :rtype: CogSourceInformationJson
        """
        cog_source_information = CogSourceInformationAPI().getById(id_cog_source_information)
        schema = CogSourceInformationSchema()
        results = schema.load(cog_source_information, many=False)
        return results

    def getCogsSourceInformationByFilterParameter(dict_parameters:dict):
        """
        get a list of cogs source information given a filters by fields E.G: dict['level']=1
        return all couples with the level at 1

        :param dict_parameters: dictionary that contain the fields and vaules to filter

        :type dict_parameters: dictionary

        :return: a json of the couple
        :rtype: list[CogSourceInformationJson]
        """
        url_parameters = ''
        for key_param in dict_parameters:
            url_parameters += key_param + '=' + str(dict_parameters[key_param]) + '&'

        url_parameters = url_parameters[:-1]
        list_cogs_source_information = CogSourceInformationAPI().getCogsSourceInformationByParameters(url_parameters)
        schema = CogSourceInformationSchema()
        results = schema.load(list_cogs_source_information, many=True)
        return results[0]