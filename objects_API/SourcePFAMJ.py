from marshmallow import Schema, fields, post_load

from rest_client.SourcePFAMRest import SourcePFAMAPI


class SourcePFAMSchema(Schema):
    """
    This class map the json into the object SourcePFAM

    ..note:: see marshmallow API
    """
    id = fields.Int()
    designation = fields.Str()

    @post_load
    def make_SourcePFAM(self, data):
        return SourcePFAMJson(**data)

class SourcePFAMJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, designation = ''):
        """
        Initialization of the class

        :param id: name of the function
        :param designation: name of the sourcePFAM

        :type id: int
        :type designation: string 

        """
        self.id = id
        self.designation = designation

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} designation {1}'.format(self.id, self.designation)

    def getAllAPI():

        """
        get all the sourcePFAM on the database

        :return: list of SourcePFAM
        :rtype: vector[SourcePFAMJ]
        """
        list_source_pfam = SourcePFAMAPI().get_all()
        schema = SourcePFAMSchema()
        results = schema.load(list_source_pfam, many=True)
        return results

    def setSourcePFAM(self):
        """
        set new domain

        :return: new domain completed with the id
        :rtype: SourcePFAMJ
        """
        schema = SourcePFAMSchema(only=['designation'])
        json_source_pfam = schema.dump(self)
        resultsCreation = SourcePFAMAPI().setSourcePFAM(jsonData = json_source_pfam.data)
        schema = SourcePFAMSchema()
        results = schema.load(resultsCreation)
        return results[0]