import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest


class DomainInteractionSourceAPI(object):
    """
    This class manage the requests for the Domain Interaction Source objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='domaininteractsource/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function


    def get_all(self):
        """
        get all the domains Interaction Sources on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def setDomainInteractionSource(self, jsonData):
        """
        set new domain interaction source  in the database

        :return: json file with the last domain created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post

    def getIdDDISource(self, id_ddi:int, id_source:int):
        """
        get teh id of a ddi source data if already exists in the database

        :param id_ddi: ID of the DDI
        :param id_source: ID of the Source

        :type id_ddi: int
        :type id_source: int

        :return: json that contain the information about the existence and the id or -1
        :rtype: string (json format)
        """


        self.function += 'ddi_info_source_existence/' + str(id_ddi) + '/' + str(id_source)

        result_get = GetRest(function = self.function).performRequest()
        return result_get