import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest

class BacteriophageAPI(object):
    """
    This class manage the requests for the Bacteriophage objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='bacteriophage/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function

    def get_all(self):
        """
        get all the Bacteriophage on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def get_by_id(self, id_bacteriophage:int):
        """
        get a bacteriophage given it id

        :param id_bacteriophage: id of th bacteriophage

        :type id_bacteriophage: int

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += str(id_bacteriophage) + '/'

        result_get = GetRest(function = self.function).performRequest()
        return result_get


    def set_bacteriophage(self, jsonData):
        """
        set new bacteriophage in the database

        :return: json file with the last bacteriophage created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post

    def setBacteriophageExistsByDesignation(self, designation):
        """
        Verify if a bacteriophage exists according a designation

        :param designation: designation name of a bacteriophage

        :type designation: string

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += 'design/' + acc_value + '/exists/'

        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def getBacteriophageByDesignation(self, designation):
        """
        Get a bacteriophage exists according a designation

        :param designation: designation name of a bacteriophage

        :type designation: string

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += 'design/' + designation + '/'

        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def getBacteriophageByAcc(self, acc_value):
        """
        return a bacteriophage according to its acc

        :param acc_value: accession number of the bacteriophage that you want to return

        :type acc_value: string

        :return: json file with all the data
        :rtype: string (json format)
        """

        self.function += 'accnumber/' + acc_value + '/'

        result_get = GetRest(function = self.function).performRequest()
        return result_get