import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest

class CogAPI(object):
    """
    This class manage the requests for the cog objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='cog/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function

    def getAll(self):
        """
        get all the cogs on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def setCog(self, jsonData):
        """
        set new cogs in the database

        :return: json file with the last genus created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post

    def getById(self, id_cog:int):
        """
        get a cog given it id

        :param id_cog: id of the cog

        :type id_cog: int

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += str(id_cog) + '/'

        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def getCogsByParameters(self, url_parameters:str):
        """
        return a list of cogs according the parameters you send

        :param url_parameters: string that contains the parameters values (that design the fields)

        :type url_parameters: str

        :return: json file with all the data
        :rtype: string (json format)
        """


        self.function += '?' + url_parameters

        result_get = GetRest(function = self.function).performRequest()
        return result_get