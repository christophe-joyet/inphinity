from marshmallow import Schema, fields, post_load

from rest_client.ContigRest import ContigAPI

class ContigSchema(Schema):
    """
    This class map the json into the object WholeDNA

    ..note:: see marshmallow API
    """
    id = fields.Int()
    id_db_online = fields.Str(required=False, allow_none=True)
    sequence_DNA = fields.Str()
    fasta_head = fields.Str()
    organism = fields.Int()
    protein_contig = fields.List(fields.Url)

    @post_load
    def make_Contig(self, data):
        return ContigJson(**data)

class ContigJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, id_db_online = None, sequence_DNA= None, fasta_head = None, organism= None, protein_contig = None):
        """
        Initialization of the class

        :param id: name of the function
        :param id_db_online: id of the database where the dna comes
        :param sequence_DNA: DNA sequence
        :param fasta_head: head of the fasta from where the DNA was extracted
        :param organism: id of the organism
        :param protein_contig: list of the proteins in the contig

        :type id: int
        :type id_db_online: string 
        :type sequence_DNA: string
        :type fasta_head: string
        :type organism: int
        :type protein_contig: List(int)

        """
        self.id = id
        self.id_db_online = id_db_online
        self.sequence_DNA = sequence_DNA
        self.fasta_head = fasta_head
        self.organism = organism
        self.protein_contig = protein_contig

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} fasta head {2}'.format(self.id, self.fasta_head)

    def getAllAPI():

        """
        get all the Contig on the database

        :return: list of Contig
        :rtype: vector[ContigJson]
        """
        list_contigs = ContigAPI().get_all()
        schema = ContigSchema()
        results = schema.load(list_contigs, many=True)
        return results

    def setContig(self):
        """
        set new contig

        :return: new contig completed with the id
        :rtype: ContigJson
        """

        schema = ContigSchema(only=['id_db_online','sequence_DNA','fasta_head','organism'])
        jsonContig = schema.dump(self)
        resultsCreation = ContigAPI().set_contig(jsonData = jsonContig.data)
        schema = ContigSchema()
        results = schema.load(resultsCreation)
        return results[0]

    def getByOrganismID(organism_id:int):
        """
        get all contigs of a given organism

        :param organism_id: organism ID

        :type organism_id: int

        :return: all the contigs of the given organism id
        :rtype: list[ContigJson]
        """

        resultsCreation = ContigAPI().getByOrganismID(organism_id = organism_id)
        schema = ContigSchema()
        results = schema.load(resultsCreation, many=True)
        return results
