import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest


class CogInteractionSourceAPI(object):
    """
    This class manage the requests for the cog score interaction source objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='coginteractsource/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function

    def getAll(self):
        """
        get all the cog interaction source on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def setCogInterSource(self, jsonData):
        """
        set new cog interaction source in the database

        :return: json file with the last couple created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post

    def getById(self, id_cog_score_interact_source:int):
        """
        get a cog interaction source given it id

        :param id_cog_interact_source: id of the cog score

        :type id_cog_interact_source: int

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += str(id_genus) + '/'

        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def getCogsInteractionSourceByParameters(self, url_parameters:str):
        """
        return a list of cog interaction source according the parameters you send

        :param url_parameters: string that contains the parameters values (that design the fields)

        :type url_parameters: str

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += '?' + url_parameters

        result_get = GetRest(function = self.function).performRequest()
        return result_get