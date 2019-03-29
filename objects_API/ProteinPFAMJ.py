from marshmallow import Schema, fields, post_load

from rest_client.ProteinPFAMRest import ProteinPFAMAPI

import datetime

class ProteinPFAMSchema(Schema):
    """
    This class map the json into the object Domain

    ..note:: see marshmallow API
    """
    id = fields.Int()
    date_creation = fields.Date()
    domain = fields.Int()
    person_responsible = fields.Int()
    protein = fields.Int()
    source = fields.Int()
    e_value = fields.Float()

    @post_load
    def make_ProteinPFAM(self, data):
        return ProteinPFAMJson(**data)

class ProteinPFAMJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, date_creation:datetime, domain:int, person_responsible:int, protein:int, source:int, e_value:float, id = None):
        """
        Initialization of the class

        :param id: id of the proteinPFAM
        :param domain_id: id of the PFAM domain
        :param person_responsible_id: id of the person_responsible
        :param protein_id: id of the protein
        :param source_id: id of the source pfam
        :param date_creation: date of the creation
        :param e-value: e-value of the match

        :type id: int
        :type domain_id: int
        :type person_responsible_id: int
        :type protein_id: int
        :type source_id: int
        :type date_creation: datetime
        :type e-value: float

        """
        self.id = id
        self.domain = domain
        self.person_responsible = person_responsible
        self.protein = protein
        self.source = source
        self.date_creation = date_creation
        self.e_value = e_value

    def __str__(self):
        """
        override the Str function 

        """
        date_creation = str(self.date_creation)

        return 'id: {0} domain_id: {1} creation date: {2}'.format(self.id, self.domain, self.date_creation)

    def getAllAPI():

        """
        get all the proteinPFAM on the database

        :return: list of proteinPFAM
        :rtype: vector[proteinPFAMJ]
        """
        list_protein_pfam = ProteinPFAMAPI().get_all()
        schema = ProteinPFAMSchema()
        results = schema.load(list_protein_pfam, many=True)
        return results

    def setProteinPFAM(self):
        """
        set new domain

        :return: new domain completed with the id
        :rtype: DomainJ
        """
        schema = ProteinPFAMSchema(only=['date_creation','domain','person_responsible','protein','source','e_value'])
        json_protein_Pfam = schema.dump(self)
        resultsCreation = ProteinPFAMAPI().setProteinPFAM(jsonData = json_protein_Pfam.data)
        schema = ProteinPFAMSchema()
        results = schema.load(resultsCreation)
        return results[0]

    def getProteinPfamByFilterParameter(dict_parameters:dict):
        """
        get a list of proteinsPfam given a filters by fields E.G: dict['level']=1
        return all proteinPfam with the level at 1

        :param dict_parameters: dictionary that contain the fields and values to filter

        :type dict_parameters: dictionary

        :return: a json of the proteinPfam
        :rtype: list[proteinPfamJson]
        """
        url_parameters = ''
        for key_param in dict_parameters:
            url_parameters += key_param + '=' + str(dict_parameters[key_param]) + '&'

        url_parameters = url_parameters[:-1]
        list_proteinPfam = ProteinPFAMAPI().getProteinsPFAMByParameters(url_parameters)
        schema = ProteinPFAMSchema()
        results = schema.load(list_proteinPfam, many=True)
        return results[0]
