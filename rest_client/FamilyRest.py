import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest


class FamilyAPI(object):
    """
    This class manage the requests for the family objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='family/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function
        

    def get_all(self):
        """
        get all the families on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def set_family(self, jsonData):
        """
        set new family in the database

        :return: json file with the last family created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post

    def get_by_id(self, id_family:int):
        """
        get a family given it id

        :param id_family: id of the family

        :type id_family: int

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += str(id_family) + '/'

        result_get = GetRest(function = self.function).performRequest()
        return result_get
