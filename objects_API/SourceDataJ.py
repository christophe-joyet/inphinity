from marshmallow import Schema, fields, post_load

from rest_client.SourceDataRest import SourceDataAPI

class SourceDataSchema(Schema):
    """
    This class map the json into the object SourceData

    ..note:: see marshmallow API
    """
    id = fields.Int()
    designation = fields.Str()

    @post_load
    def make_SourceData(self, data):
        return SourceDataJson(**data)

class SourceDataJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, designation = ''):
        """
        Initialization of the class

        :param id: name of the function
        :param designation: name of the Source data

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
        get all the Source Data on the database

        :return: list of Source data
        :rtype: vector[SrouceData]
        """
        list_source_data = SourceDataAPI().get_all()
        schema = SourceDataSchema()
        results = schema.load(list_source_data, many=True)
        return results

    def setSourceData(self):
        """
        set new SourceData

        :return: new Source Data completed with the id
        :rtype: SourceDataJ
        """
        schema = SourceDataSchema(only=['designation'])
        jsonSourData = schema.dump(self)
        resultsCreation = SourceDataAPI().set_sourceData(jsonData = jsonSourData.data)
        schema = SourceDataSchema()
        results = schema.load(resultsCreation)
        return results[0]