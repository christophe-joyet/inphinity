import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest

class ProteinAPI(object):
    """
    This class manage the requests for the Proteins objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='protein/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function

    def get_all(self):
        """
        get all the Proteins on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def set_protein(self, jsonData):
        """
        set new protein in the database

        :return: json file with the last protein created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post

    def getByOrganismID(self, organism_id):
        """
        get all proteins of a given organism

        :param organism_id: organism ID

        :type organism_id: int

        :return: all the proteins of the given organism id
        :rtype: ProteinJson
        """

        self.function += 'organism_id/' + str(organism_id)

        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def getProteinsByParameters(self, url_parameters:str):
        """
        return a proteins according the parameters you send

        :param url_parameters: string that contains the parameters values (that design the fields)

        :type url_parameters: str

        :return: json file with all the data
        :rtype: string (json format)
        """


        self.function += '?' + url_parameters

        result_get = GetRest(function = self.function).performRequest()
        return result_get