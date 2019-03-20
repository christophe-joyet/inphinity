from marshmallow import Schema, fields, post_load

from rest_client.WholeDNARest import WholeDNAAPI

class WholeDNASchema(Schema):
    """
    This class map the json into the object WholeDNA

    ..note:: see marshmallow API
    """
    id = fields.Int()
    id_db_online = fields.Str(required=False, allow_none=True)
    sequence_DNA = fields.Str()
    fasta_head = fields.Str()
    organism = fields.Int()

    @post_load
    def make_WholeDNA(self, data):
        return WholeDNAJson(**data)

class WholeDNAJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, id_db_online = None, sequence_DNA= None, fasta_head = None, organism= None):
        """
        Initialization of the class

        :param id: name of the function
        :param id_db_online: id of the database where the dna comes
        :param sequence_DNA: DNA sequence
        :param fasta_head: head of the fasta from where the DNA was extracted
        :param organism: id of the organism

        :type id: int
        :type id_db_online: string 
        :type sequence_DNA: string
        :type fasta_head: string
        :type organism: int


        """
        self.id = id
        self.id_db_online = id_db_online
        self.sequence_DNA = sequence_DNA
        self.fasta_head = fasta_head
        self.organism = organism

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} fasta head {1}'.format(self.id, self.fasta_head)

    def getAllAPI():

        """
        get all the WholeDNAs on the database

        :return: list of WholeDNAs
        :rtype: vector[WholeDNAJson]
        """
        list_wholeDNA = WholeDNAAPI().get_all()
        schema = WholeDNASchema()
        results = schema.load(list_wholeDNA, many=True)
        return results

    def setWholeDNA(self):
        """
        set new wholeDNA

        :return: new wholeDNA completed with the id
        :rtype: WholeDNAJson
        """

        schema = WholeDNASchema(only=['id_db_online','sequence_DNA','fasta_head','organism'])
        jsonWholeDNA = schema.dump(self)
        resultsCreation = WholeDNAAPI().set_wholeDNA(jsonData = jsonWholeDNA.data)
        schema = WholeDNASchema()
        results = schema.load(resultsCreation)
        return results[0]

    def getByOrganismID(organism_id:int):
        """
        get the whole_dna of a given organism

        :param organism_id: organism ID

        :type organism_id: int

        :return: the whole_dna of the given organism id
        :rtype: WholeDna
        """

        resultsJson = WholeDNAAPI().getByOrganismID(organism_id = organism_id)
        schema = WholeDNASchema()
        results = schema.load(resultsJson, many=False)
        return results[0]