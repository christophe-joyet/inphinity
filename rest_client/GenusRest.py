import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest

class GenusAPI(object):
    """
    This class manage the requests for the genus objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='genus/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function

    def get_all(self):
        """
        get all the genus on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def set_genus(self, jsonData):
        """
        set new genus in the database

        :return: json file with the last genus created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post

    def get_by_id(self, id_genus:int):
        """
        get a genus given it id

        :param id_genus: id of the genus

        :type id_genus: int

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += str(id_genus) + '/'

        result_get = GetRest(function = self.function).performRequest()
        return result_get