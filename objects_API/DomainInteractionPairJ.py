from marshmallow import Schema, fields, post_load

from rest_client.DomainInteractionPairRest import DomainInteractionPairAPI


class DomainInteractionPairSchema(Schema):
    """
    This class map the json into the object Domain

    ..note:: see marshmallow API
    """
    id = fields.Int()
    domain_a = fields.Int()
    domain_b = fields.Int()

    @post_load
    def make_DomainInteractionPair(self, data):
        return DomainInteractionPairJson(**data)

class DomainInteractionPairJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self,  domain_a:int, domain_b:int, id = None):
        """
        Initialization of the class

        :param id: name of the function
        :param domain_a_id: name of the domain
        :param domain_b_id: name of the domain

        :type id: int
        :type domain_a_id: int 
        :type domain_b_id: int 

        """
        self.id = id
        self.domain_a = domain_a
        self.domain_b = domain_b

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} id domain A: {1} id domain B: {2}'.format(self.id, self.domain_a, self.domain_b)

    def getAllAPI():

        """
        get all the domain on the database

        :return: list of DomainInteractionPairJ
        :rtype: vector[DomainInteractionPairJ]
        """
        list_domain_interaction_pair = DomainInteractionPairAPI().get_all()
        schema = DomainInteractionPairSchema()
        results = schema.load(list_domain_interaction_pair, many=True)
        return results

    def setDomainInteractionPair(self):
        """
        set new domain

        :return: new domain completed with the id
        :rtype: DomainJ
        """
        schema = DomainInteractionPairSchema(only=['domain_a','domain_b'])
        json_domain = schema.dump(self)
        resultsCreation = DomainInteractionPairAPI().setDomainInteractionPair(jsonData = json_domain.data)
        schema = DomainInteractionPairSchema()
        results = schema.load(resultsCreation)
        return results[0]

    def verifyDDIpairExistence(pfam_a:str, pfam_b:str):
        """
        verify if the ddi pair already exists in the database and return it's ID or -1 in case of inexistence

        :param pfam_a: domain a
        :param pfam_b: domain b

        :type pfam_a: string
        :type pfam_b: string

        :return: id of the ddi or -1
        :rtype: int
        """

        results_DDI = DomainInteractionPairAPI().getIdDDI(pfam_a = pfam_a, pfam_b = pfam_b)
        id_ddi = results_DDI['id_ddi_interaction']
        return id_ddi