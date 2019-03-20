from marshmallow import Schema, fields, post_load

from rest_client.PersonResponsibleRest import PersonResponsibleAPI

class PersonResponsibleSchema(Schema):
    """
    This class map the json into the object Person Responsible

    ..note:: see marshmallow API
    """
    id = fields.Int()
    name = fields.Str()

    @post_load
    def make_PersonResponsible(self, data):
        return PersonResponsibleJson(**data)

class PersonResponsibleJson(object):
    """
    This class manage the object and is used to map them into json format
    """

    def __init__(self, id = None, name = ''):
        """
        Initialization of the class

        :param id: name of the function
        :param name: name of the Person Responsible

        :type id: int
        :type name: string 

        """
        self.id = id
        self.name = name

    def __str__(self):
        """
        override the Str function 

        """
        return 'id: {0} name {1}'.format(self.id, self.name)

    def getAllAPI():

        """
        get all the Person Responsible on the database

        :return: list of Person Responsible
        :rtype: vector[PersonResponsibleJson]
        """
        list_person_responsible = PersonResponsibleAPI().get_all()
        schema = PersonResponsibleSchema()
        results = schema.load(list_person_responsible, many=True)
        return results

    def setPersonResponsible(self):
        """
        set new Person Responsible

        :return: new Person Responsible completed with the id
        :rtype: PersonResponsibleJson
        """
        schema = PersonResponsibleSchema(only=['name'])
        jsonPerResp = schema.dump(self)
        resultsCreation = PersonResponsibleAPI().set_person_responsible(jsonData = jsonPerResp.data)
        schema = PersonResponsibleSchema()
        results = schema.load(resultsCreation)
        return results[0]