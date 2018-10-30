from marshmallow import Schema, fields, post_load

from rest_client.ProteinRest import ProteinAPI

class ProteinSchema(Schema):
    """
    This class map the json into the object Protein

    ..note:: see marshmallow API
    """
    id = fields.Int()
    id_db_online = fields.Str(allow_none = True)
    organism = fields.Int()
    gene = fields.Int(allow_none = True)
    contig = fields.Int(allow_none = True)
    sequence_AA = fields.Str()
    fasta_head = fields.Str(allow_none = True)
    description = fields.Str(allow_none = True)
    accession_num = fields.Str(allow_none = True)
    position_start = fields.Int()
    position_end = fields.Int()
    position_start_contig = fields.Int(allow_none = True)
    position_end_contig = fields.Int(allow_none = True)



    @post_load
    def make_ProteinJson(self, data):
        return ProteinJson(**data)

class ProteinJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, id_db_online = None, organism= None, gene = None, contig= None, sequence_AA = None, fasta_head = None, description = None, accession_num = None, position_start = None, position_end = None, position_start_contig = None, position_end_contig = None):
        """
        Initialization of the class

        :param id: of the protein
        :param id_db_online: id of the database where the dna comes
        :param organism: FK of the organism
        :param gene: FK of the gene
        :param contig: FK of the contig
        :param sequence_AA: Sequence AA of the protein
        :param fasta_head: Head of the fasta line
        :param description: description of the protein function
        :param accession_num: ACC of the protein
        :param position_start: position start on the gene
        :param position_end: position end on the gene
        :param position_start_contig: position end on the contig
        :param position_end_contig: position end on the contig

        :type id: int
        :type id_db_online: string 
        :type organism: int
        :type gene: int
        :type contig: int
        :type sequence_AA: string
        :type fasta_head: string
        :type description: string
        :type accession_num: string
        :type position_start: int
        :type position_end: int
        :type position_start_contig: int
        :type position_end_contig: int


        """
        self.id = id
        self.id_db_online = id_db_online
        self.organism = organism
        self.gene = gene
        self.contig = contig
        self.sequence_AA = sequence_AA
        self.fasta_head = fasta_head
        self.description = description
        self.accession_num = accession_num
        self.position_start = position_start
        self.position_end = position_end
        self.position_start_contig = position_start_contig
        self.position_end_contig = position_end_contig

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} fasta head {1} organism {2}'.format(self.id, self.fasta_head, self.organism)

    def getAllAPI():

        """
        get all the Proteins on the database

        :return: list of Proteins
        :rtype: vector[ProteinJson]
        """
        list_protein = ProteinAPI().get_all()
        schema = ProteinSchema()
        results = schema.load(list_protein, many=True)
        return results[0]

    def setProtein(self):
        """
        set new protein

        :return: new protein completed with the id
        :rtype: ProteinJson
        """

        schema = ProteinSchema(only=['id_db_online','accession_num','organism','gene','contig','sequence_AA','fasta_head', 'description','position_start','position_end','position_start_contig','position_end_contig'])
        jsonProtein = schema.dump(self)
        resultsCreation = ProteinAPI().set_protein(jsonData = jsonProtein.data)
        schema = ProteinSchema()
        results = schema.load(resultsCreation)
        return results[0]