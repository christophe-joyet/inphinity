from marshmallow import Schema, fields, post_load

from rest_client.CoupleRest import CoupleAPI


class CoupleSchema(Schema):
    """
    This class map the json into the object Family

    ..note:: see marshmallow API
    """

    id = fields.Int()
    bacterium = fields.Int()
    bacteriophage = fields.Int()
    interaction_type = fields.Bool()
    level = fields.Int()
    lysis = fields.Int(required=False, allow_none=True)
    source_data = fields.Int()
    person_responsible = fields.Int()
    validity = fields.Int()

    @post_load
    def make_Couple(self, data):
        return CoupleJson(**data)

class CoupleJson(object):
    """
    This class manage the object and is used to map them into json format
    """
    def __init__(self, id:int = None, interaction_type:int = None, bacteriophage:int = None, bacterium:int = None, level:int = None, lysis:int = None, person_responsible:int = None, source_data:int = None, validity:int = None):

        """
        Initialization of the class

        :param id: id of the couple
        :param interaction_type: type of the interaction (True or False)
        :param bacteriophage_id: ID of the bacteriophage
        :param bacterium_id: ID of the bacterium
        :param level_id: level of the interaction (in terms of taxonomy) (between 0 and 4)
        :param person_responsible_id: ID of the person who has inserted the couple
        :param source_data_id: ID of the source where you know the couple
        :param validity_id: ID of the couple validity (verified, in vivo, generated,...)

        :type id: int
        :type interaction_type: int 
        :type bacteriophage_id: int
        :type bacterium_id: int
        :type level_id: int
        :type person_responsible_id: int
        :type source_data_id: int
        :type validity_id: int

        """
        self.id = id
        self.interaction_type = interaction_type
        self.bacteriophage = bacteriophage
        self.bacterium = bacterium
        self.level = level
        self.lysis = lysis
        self.person_responsible = person_responsible
        self.source_data = source_data
        self.validity = validity

    def getAllAPI():

        """
        get all the Couples on the database

        :return: list of couple
        :rtype: vector[CoupleJson]
        """
        list_couple = CoupleAPI().get_all()
        schema = CoupleSchema()
        results = schema.load(list_couple, many=True)
        return results

    def getByBacteriumPhageIds(bacterium_id:int, phage_id:int):

        """
        get a couple given the id of the bacterium and phage

        :param bact_id: If of the bacterium
        :param phage_id: If of the bacteriophage

        :type bact_id: int
        :type phage_id: int

        :return: a json of the couple
        :rtype: CoupleJson
        """
        couple = CoupleAPI().getCoupleByBactPhageIds(bacterium_id, phage_id)
        schema = CoupleSchema()
        results = schema.load(couple, many=False)
        return results[0]


    def setCouple(self):
        schema = CoupleSchema(only=['interaction_type','bacteriophage','bacterium','level','person_responsible','source_data','validity'])
        jsonCouple = schema.dump(self)
        resultsCouple = CoupleAPI().set_couple(jsonData = jsonCouple.data)
        schema = CoupleSchema()
        results = schema.load(resultsCouple)
        return results[0]

    def setCoupleWithLysis(self):
        schema = CoupleSchema(only=['interaction_type','bacteriophage','bacterium','level','person_responsible','source_data','validity','lysis'])
        jsonCouple = schema.dump(self)
        resultsCouple = CoupleAPI().set_couple(jsonData = jsonCouple.data)
        schema = CoupleSchema()
        results = schema.load(resultsCouple)
        return results[0]