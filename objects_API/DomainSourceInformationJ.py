from marshmallow import Schema, fields, post_load

from rest_client.DomainSourceInformationRest import DomainSourceInformationAPI

class DomainSourceInformationSchema(Schema):
    """
    This class map the json into the object Domain

    ..note:: see marshmallow API
    """
    id = fields.Int()
    designation = fields.Str()

    @post_load
    def make_DomainSourceInformation(self, data):
        return DomainSourceInformationJson(**data)

class DomainSourceInformationJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, designation = ''):
        """
        Initialization of the class

        :param id: id of the domain source information
        :param designation: name of the domain source information

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
        get all the domain source Information on the database

        :return: list of Domain
        :rtype: vector[DomainJ]
        """
        list_domain_source_information = DomainSourceInformationAPI().get_all()
        schema = DomainSourceInformationSchema()
        results = schema.load(list_domain_source_information, many=True)
        return results

    def setDomainSourceInformation(self):
        """
        set new domain

        :return: new domain completed with the id
        :rtype: DomainJ
        """
        schema = DomainSourceInformationSchema(only=['designation'])
        json_domain_source_information = schema.dump(self)
        resultsCreation = DomainSourceInformationAPI().setDomain(jsonData = json_domain_source_information.data)
        schema = DomainSchema()
        results = schema.load(resultsCreation)
        return results[0]