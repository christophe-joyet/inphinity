from marshmallow import Schema, fields, post_load
import datetime


from rest_client.DomainInteractionSourceRest import DomainInteractionSourceAPI

class DomainInteractionSourceSchema(Schema):
    """
    This class map the json into the object Domain Interaction Source

    ..note:: see marshmallow API
    """
    id = fields.Int()
    date_creation = fields.Date()
    domain_interaction = fields.Int()
    information_source = fields.Int()

    @post_load
    def make_DomainInteractionSource(self, data):
        return DomainInteractionSourceJson(**data)

class DomainInteractionSourceJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self,  date_creation:datetime, domain_interaction:int, information_source:int, id = None):
        """
        Initialization of the class

        :param id: name of the function
        :param date_creation: domain interaction source creation date
        :param domain_interaction_id: id of the domain interaction
        :param information_source_id: id of the information source

        :type id: int
        :type date_creation: date 
        :type domain_interaction_id: int 
        :type information_source_id: int 

        """
        self.id = id
        self.date_creation = date_creation
        self.domain_interaction = domain_interaction
        self.information_source = information_source

    def getAllAPI():

        """
        get all the domain interaction source on the database

        :return: list of Domain interaction source
        :rtype: vector[DomainInteractionSourceJson]
        """
        list_domain_interaction_source = DomainInteractionSourceAPI().get_all()
        schema = DomainInteractionSourceSchema()
        results = schema.load(list_domain_interaction_source, many=True)
        return results

    def setDomainInteractionSource(self):
        """
        set new domain interaction source

        :return: new domain interaction source completed with the id
        :rtype: DomainJ
        """
        schema = DomainInteractionSourceSchema(only=['date_creation','domain_interaction','information_source'])
        json_domain = schema.dump(self)
        resultsCreation = DomainInteractionSourceAPI().setDomainInteractionSource(jsonData = json_domain)
        schema = DomainInteractionSourceSchema()
        results = schema.load(resultsCreation)
        return results


    def verifyDDIpairSourceExistence(id_ddi:int, id_source:int):
        """
        verify if the ddi pair source already exists in the database and return it's ID or -1 in case of inexistence

        :param id_ddi: domain a
        :param id_source: domain b

        :type id_ddi: int
        :type id_source: int

        :return: id of the ddi source existence or -1
        :rtype: int
        """

        results_DDI = DomainInteractionSourceAPI().getIdDDISource(id_ddi = id_ddi, id_source = id_source)
        id_ddi = results_DDI['id_ddi_iteract_source']
        return id_ddi


    def getDomainInteractionSourceByFilterParameter(dict_parameters:dict):
        """
        get a list of domains interactions sources given a filters by fields E.G: dict['domain_interaction']=10
        return all domain interaction source for the  domains pair id n 10

        :param dict_parameters: dictionary that contain the fields and values to filter

        :type dict_parameters: dictionary

        :return: a json of the domain interaction source
        :rtype: list[DomainInteractionSourceJ]
        """
        url_parameters = ''
        for key_param in dict_parameters:
            url_parameters += key_param + '=' + str(dict_parameters[key_param]) + '&'

        url_parameters = url_parameters[:-1]
        list_domains_interact_source = DomainInteractionSourceAPI().getIdDDISByParameters(url_parameters)
        schema = DomainInteractionSourceSchema()
        results = schema.load(list_domains_interact_source, many=True)
        return results

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} creation date {1} interaction domain id {2} source id {3}'.format(self.id, self.date_creation, self.domain_interaction, self.information_source)