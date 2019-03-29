from marshmallow import Schema, fields, post_load
import json

from rest_client.BacteriumRest import BacteriumAPI
from objects_API.OrganismJ import *

class BacteriumSchema(OrganismSchema):
    """
    This class map the json into the object Bacterium

    ..note:: see marshmallow API
    """
    
    #organism = fields.Nested('OrganismSchema', many = False, required = True)
    id = fields.Int()
    strain = fields.Int()


    @post_load
    def make_Bacterium(self, data):

        return BacteriumJson(**data)

class BacteriumJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, acc_number = None, gi_number= None, person_responsible = None, source_data = None, organism_gene = None, organism_contig = None, organism_wholeDNA = None, protein_organism = None, strain = None):
        """
        Initialization of the class

        :param id: name of the function
        :param strain: name of the function

        :type id: int
        :type strain: string 


        """

        OrganismJson.__init__(self,acc_number, gi_number, source_data, person_responsible, organism_gene, organism_contig, organism_wholeDNA, protein_organism)
        self.id = id
        self.strain = strain

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} strain {1}'.format(self.id, self.strain)

    def getAllAPI():

        """
        get all the Bacteria on the database

        :return: list of Bacteria
        :rtype: vector[BacteriuJson]
        """
        list_bacteria = BacteriumAPI().get_all()
        schema = BacteriumSchema()
        results = schema.load(list_bacteria, many=True)
        return results

    def setBacterium(self):
        """
        set new bacterium

        :return: new bacterium completed with the id
        :rtype: BacteriumJson
        """

        schema = BacteriumSchema(only=['acc_number','gi_number','source_data','person_responsible','strain'])
        jsonBacterium = schema.dump(self)
        resultsCreation = BacteriumAPI().set_bacterium(jsonData = jsonBacterium.data)
        schema = BacteriumSchema()
        results = schema.load(resultsCreation)
        return results[0]

    def getByID(id_bacterium:int):

        """
        get a bacterium given its id

        :param id_bacterium: id of the bacterium that it will be returned

        :type id_bacteriophage: int

        :return: a json of the bacterium 
        :rtype: BacteriumJson
        """
        bacterium = BacteriumAPI().get_by_id(id_bacterium)
        schema = BacteriumSchema()
        results = schema.load(bacterium, many=False)
        return results

    def verifiyBacteriumExistanceByAcc(acc_value):
        """
        Verify if a bacterium exists according an accValue

        :param acc_value: accession number of the bacterium that you want to check the existence

        :type acc_value: string

        :return: return True or False according the existence
        :rtype: boolean
        """

        resultsCreation = BacteriumAPI().setBacteriumExistsByAcc(acc_value = acc_value)
        bacterium_existence = resultsCreation['value']['bacterium_exists']
        return bacterium_existence

    def getByAccnumber(acc_number:str):

        """
        get a bacterium given its acc

        :param acc_number: acc of the bacterium that it will be returned

        :type acc_number: String

        :return: a json of the bacterium
        :rtype: BacteriumJson
        """
        bacterium = BacteriumAPI().getBacteriumByAcc(acc_number)
        schema = BacteriumSchema()
        results = schema.load(bacterium, many=False)
        return results[0]

