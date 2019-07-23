from marshmallow import Schema, fields, post_load

from rest_client.PPICogScoreRest import PPICogScoreAPI


class PPIScogScoreSchema(Schema):
    """
    This class map the json into the object PPi cog score

    ..note:: see marshmallow API
    """
    id = fields.Int()
    ppi_interaction = fields.Int()
    cog_score = fields.Int()
    date_creation = fields.Date()

    @post_load
    def make_PPICogScore(self, data):
        return PPIScogScoreSchema(**data)

class PPICogScoreJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, ppi_interaction, cog_score, date_creation):
        """
        Initialization of the class

        :param id: id of the ppi cog score
        :param ppi_interaction: id of the ppi interaction
        :param cog_score: id of the cog score
        :param date_creation: creation date (year, month, day)

        :type id: int
        :type ppi_interaction: int
        :type cog_score: int
        :type date_creation: Date 

        """
        self.id = id
        self.ppi_interaction = ppi_interaction
        self.cog_score = cog_score
        self.date_creation = date_creation

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} ppi_interaction {1} cog_score {2} date_creation {3}'.format(self.id, self.ppi_interaction, self.cog_score, self.date_creation)

    def getAllAPI():

        """
        get all the ppi cog scores on the database

        :return: list of PPI cog score
        :rtype: vector[PPICogScoreJ]
        """
        list_ppi_cog_score = PPICogScoreAPI().getAll()
        schema = PPIScogScoreSchema()
        results = schema.load(list_ppi_cog_score, many=True)
        return results

    def setPPICogScore(self):
        """
        set new PPI cog score

        :return: new cog completed with the id
        :rtype: PPICogScoreJ
        """
        schema = PPIScogScoreSchema()

        json_ppi_cog_score = schema.dump(self)
        resultsCreation = PPICogScoreAPI().setPPICogScore(jsonData = json_cog)
        schema = PPIScogScoreSchema()
        results = schema.load(resultsCreation)
        return results


    def getByID(id_ppi_cog_score:int):

        """
        get a PPI cog score given its id

        :param id_ppi_cog_score: id of the cog that it will be returned

        :type id_ppi_cog_score: int

        :return: a json of the PPICogScore 
        :rtype: PPICogScoreJson
        """
        ppi_cog_score = PPICogScoreAPI().getById(id_ppi_cog_score)
        schema = PPIScogScoreSchema()
        results = schema.load(cog, many=False)
        return results


    def getPPICogScoreByFilterParameter(dict_parameters:dict):
        """
        get a list of PPI cog score given a filters by fields E.G: dict['level']=1
        return all couples with the level at 1

        :param dict_parameters: dictionary that contain the fields and vaules to filter

        :type dict_parameters: dictionary

        :return: a json of the PPI cog score
        :rtype: list[PPICogScoreJson]
        """
        url_parameters = ''
        for key_param in dict_parameters:
            url_parameters += key_param + '=' + str(dict_parameters[key_param]) + '&'

        url_parameters = url_parameters[:-1]
        list_ppi_cog_score = PPICogScoreAPI().getPPICogsScoreByParameters(url_parameters)
        schema = PPIScogScoreSchema()
        results = schema.load(list_ppi_cog_score, many=True)
        return results[0]