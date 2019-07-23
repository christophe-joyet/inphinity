from marshmallow import Schema, fields, post_load

from rest_client.DomainRest import DomainAPI

class DomainSchema(Schema):
    """
    This class map the json into the object Domain

    ..note:: see marshmallow API
    """
    id = fields.Int()
    designation = fields.Str()

    @post_load
    def make_Domain(self, data):
        return DomainJson(**data)

class DomainJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, designation = ''):
        """
        Initialization of the class

        :param id: name of the function
        :param designation: name of the domain

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
        get all the domain on the database

        :return: list of Domain
        :rtype: vector[DomainJ]
        """
        list_domain = DomainAPI().get_all()
        schema = DomainSchema()
        results = schema.load(list_domain, many=True)
        return results

    def setDomain(self):
        """
        set new domain

        :return: new domain completed with the id
        :rtype: DomainJ
        """
        schema = DomainSchema(only=['designation'])
        json_domain = schema.dump(self)
        resultsCreation = DomainAPI().setDomain(jsonData = json_domain.data)
        schema = DomainSchema()
        results = schema.load(resultsCreation)
        return results[0]
