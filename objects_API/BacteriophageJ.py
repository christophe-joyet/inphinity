from marshmallow import Schema, fields, post_load

from rest_client.BacteriophageRest import BacteriophageAPI
from objects_API.OrganismJ import *

class BacteriophageSchema(OrganismSchema):
    """
    This class map the json into the object Bacteriophage

    ..note:: see marshmallow API
    """
    
    #organism = fields.Nested('OrganismSchema', many = False, required = True)
    id = fields.Int()
    baltimore_classification = fields.Int(required=False, allow_none=True)
    designation = fields.String()


    @post_load
    def make_Bacteriophage(self, data):
        return BacteriophageJson(**data)

class BacteriophageJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, acc_number = None, gi_number= None, person_responsible = None, source_data = None, organism_gene = None, organism_contig = None, organism_wholeDNA = None, protein_organism = None, baltimore_classification = None, designation = None):
        """
        Initialization of the class

        :param id: name of the function
        :param baltimore_classification: name of the function

        :type id: int
        :type baltimore_classification: int 

        """
        OrganismJson.__init__(self,acc_number, gi_number, source_data, person_responsible, organism_gene, organism_contig, organism_wholeDNA, protein_organism)
        self.id = id
        self.baltimore_classification = baltimore_classification
        self.designation = designation

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} baltimore class {1}'.format(self.id, self.baltimore_classification)

    def getAllAPI():

        """
        get all the Bacteria on the database

        :return: list of Bacteriophage
        :rtype: vector[BacteriophageJson]
        """
        list_bacteriophages = BacteriophageAPI().get_all()
        schema = BacteriophageSchema()
        results = schema.load(list_bacteriophages, many=True)
        return results[0]

    def getByID(id_bacteriophage:int):

        """
        get a bacteriophage given its id

        :param id_bacteriophage: id of the bacteriophage that it will be returned

        :type id_bacteriophage: int

        :return: a json of the bacteriophage
        :rtype: BacteriophageJson
        """
        bacteriophage = BacteriophageAPI().get_by_id(id_bacteriophage)
        schema = BacteriophageSchema()
        results = schema.load(bacteriophage, many=False)
        return results[0]

    def setBacteriophage(self):
        """
        set new bacteriophage

        :return: new bacteriophage completed with the id
        :rtype: BacteriophageJson

        """

        schema = BacteriophageSchema(only=['acc_number','gi_number','source_data','person_responsible','baltimore_classification', 'designation'])
        jsonBacteriophage = schema.dump(self)
        resultsCreation = BacteriophageAPI().set_bacteriophage(jsonData = jsonBacteriophage.data)
        schema = BacteriophageSchema()
        results = schema.load(resultsCreation)
        return results[0]


    def verifiyBacteriophageExistanceByDesignation(designation):
        """
        Verify if a bacteriophage exists according a designation

        :param designation: accession number of the bacteriophage that you want to check the existence

        :type designation: string

        :return: return True or False according the existence
        :rtype: boolean
        """

        resultsCreation = BacteriophageAPI().getBacteriophageByAcc(acc_value = acc_value)
        bacteriophage_existence = resultsCreation['value']['bacteriophage_exists']
        return bacteriophage_existence

    def getBydesignation(designation:str):

        """
        get a bacteriophage given its designation

        :param designation: designation of the bacteriophage that it will be returned

        :type designation: String

        :return: a json of the bacteriophage
        :rtype: BacteriophageJson
        """
        bacteriophage = BacteriophageAPI().getBacteriophageByDesignation(designation)
        schema = BacteriophageSchema()
        results = schema.load(bacteriophage, many=False)
        return results[0]

    def getByAccnumber(acc_number:str):

        """
        get a bacteriophage given its acc

        :param acc_number: acc of the bacteriophage that it will be returned

        :type acc_number: String

        :return: a json of the bacteriophage
        :rtype: BacteriophageJson
        """
        bacteriophage = BacteriophageAPI().getBacteriophageByAcc(acc_number)
        schema = BacteriophageSchema()
        results = schema.load(bacteriophage, many=False)
        return results[0]