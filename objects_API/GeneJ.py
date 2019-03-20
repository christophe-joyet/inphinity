from marshmallow import Schema, fields, post_load

from rest_client.GeneRest import GeneAPI

class GeneSchema(Schema):
    """
    This class map the json into the object Gene

    ..note:: see marshmallow API
    """
    id = fields.Int()
    id_db_online = fields.Str(required=False, allow_none=True)
    sequence_DNA = fields.Str()
    fasta_head = fields.Str()
    position_start = fields.Int(required=False, allow_none=True)
    position_end = fields.Int(required=False, allow_none=True)
    number_of_seq = fields.Int(required=False, allow_none=True)
    organism = fields.Int()
    protein_gene = fields.List(fields.Url)
    description = fields.Str(required=False, allow_none=True)

    contig_gene = fields.Int(allow_none = True)
    position_start_contig = fields.Int(allow_none = True)
    position_end_contig = fields.Int(allow_none = True)

    @post_load
    def make_GeneDNA(self, data):
        return GeneJson(**data)

class GeneJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, id_db_online = None, sequence_DNA= None, fasta_head = None, position_start = None, position_end = None, number_of_seq=None, organism= None, protein_gene = None, position_start_contig = None, position_end_contig = None, contig = None):
        """
        Initialization of the class

        :param id: name of the function
        :param id_db_online: id of the database where the dna comes
        :param sequence_DNA: DNA sequence
        :param fasta_head: head of the fasta from where the DNA was extracted
        :param position_start: Start position of the gene in the wholeDNA
        :param position_end: End position of the gene in the wholeDNA
        :param number_of_seq: Number of the sequenced gene (first, second,...)
        :param organism: id of the organism
        :param protein_gene: List of the proteins in this gene
        :param contig: contig of the gene
        :param position_start_contig: Start position of the gene on the contig
        :param position_end_contig: End position of the gene on the contig

        :type id: int
        :type id_db_online: string 
        :type sequence_DNA: string
        :type fasta_head: string
        :type position_start: int
        :type position_end: int
        :type number_of_seq: int
        :type organism: int
        :type protein_gene: Array[URL]
        :type contig: int
        :type position_start_contig: int
        :type position_end_contig: int


        """
        self.id = id
        self.id_db_online = id_db_online
        self.sequence_DNA = sequence_DNA
        self.fasta_head = fasta_head
        self.position_start = position_start
        self.position_end = position_end
        self.number_of_seq = number_of_seq
        self.organism = organism
        self.contig = contig
        self.position_start_contig = position_start_contig
        self.position_end_contig = position_end_contig

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} fasta head {2}'.format(self.id, self.fasta_head)

    def getAllAPI():

        """
        get all the Genes on the database

        :return: list of Genes
        :rtype: vector[GeneJson]
        """
        list_genes = GeneAPI().get_all()
        schema = GeneSchema()
        results = schema.load(list_genes, many=True)
        return results

    def setGene(self):
        """
        set new gene

        :return: new gene completed with the id
        :rtype: GeneJson
        """

        schema = GeneSchema(only=['id_db_online','sequence_DNA','fasta_head','position_start','position_end','organism','contig','position_start_contig','position_end_contig'])
        jsonGene = schema.dump(self)
        resultsCreation = GeneAPI().set_gene(jsonData = jsonGene.data)
        schema = GeneSchema()
        results = schema.load(resultsCreation)
        return results[0]