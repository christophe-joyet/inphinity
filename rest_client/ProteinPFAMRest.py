import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest

class ProteinPFAMAPI(object):
    """
    This class manage the requests for the ProteinPfam objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='proteinpfam/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function

    def get_all(self):
        """
        get all the proteinPfam on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def setProteinPFAM(self, jsonData):
        """
        set new proteinPFAM in the database

        :return: json file with the last proteinpfam created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post

    def getProteinsPFAMByParameters(self, url_parameters:str):
        """
        return a proteinsPfam according the parameters you send

        :param url_parameters: string that contains the parameters values (that design the fields)

        :type url_parameters: str

        :return: json file with all the data
        :rtype: string (json format)
        """


        self.function += '?' + url_parameters

        result_get = GetRest(function = self.function).performRequest()
        return result_get