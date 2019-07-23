import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest


class CogInteractionPairAPI(object):
    """
    This class manage the requests for the cog interaction pair objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='coginteractpair/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function

    def getAll(self):
        """
        get all the cogs interaction pair on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def setCog(self, jsonData):
        """
        set new cogs interaction pair in the database

        :return: json file with the last genus created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post

    def getById(self, id_cog_interaction_pair:int):
        """
        get a cog interaction pair given it id

        :param id_cog_interaction_pair: id of the cog interaction pair

        :type id_cog_interaction_pair: int

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += str(id_cog_interaction_pair) + '/'

        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def getCogsInteractionsPairsByParameters(self, url_parameters:str):
        """
        return a list of cogs interaction pair according the parameters you send

        :param url_parameters: string that contains the parameters values (that design the fields)

        :type url_parameters: str

        :return: json file with all the data
        :rtype: string (json format)
        """


        self.function += '?' + url_parameters

        result_get = GetRest(function = self.function).performRequest()
        return result_get