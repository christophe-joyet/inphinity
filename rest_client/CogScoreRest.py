import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest

class CogScoreAPI(object):

    """
    This class manage the requests for the cogScore objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='cogscore/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function

    def getAll(self):
        """
        get all the cogs Score on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def setCogScore(self, jsonData):
        """
        set new cogs Score in the database

        :return: json file with the last genus created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post

    def getById(self, id_cog_score:int):
        """
        get a Cog  Score given it id

        :param id_cog_score: id of the cog score

        :type id_cog_score: int

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += str(id_genus) + '/'

        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def getCogsScoreByParameters(self, url_parameters:str):
        """
        return a list of cogs Score according the parameters you send

        :param url_parameters: string that contains the parameters values (that design the fields)

        :type url_parameters: str

        :return: json file with all the data
        :rtype: string (json format)
        """


        self.function += '?' + url_parameters

        result_get = GetRest(function = self.function).performRequest()
        return result_get