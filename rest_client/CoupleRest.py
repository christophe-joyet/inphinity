import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest


class CoupleAPI(object):
    """
    This class manage the requests for the couple objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='couple/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function

    def get_all(self):
        """
        get all the couples on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def set_couple(self, jsonData):
        """
        set new couple in the database

        :return: json file with the last couple created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post

    def getCoupleByBactPhageIds(self, bact_id:int, phage_id:int):
        """
        return a couple according the ids of the bacterium and phage

        :param bact_id: If of the bacterium
        :param phage_id: If of the bacteriophage

        :type bact_id: int
        :type phage_id: int

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += 'organismsid/' + str(bact_id) + '/' + str(phage_id) + '/'

        result_get = GetRest(function = self.function).performRequest()
        return result_get
