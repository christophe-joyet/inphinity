# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, '../inphinity')

from objects_API.CoupleJ import CoupleJson
from objects_API.StrainJ import StrainJson
from objects_API.GenusJ import GenusJson
from objects_API.FamilyJ import FamilyJson
from objects_API.SourceDataJ import SourceDataJson
from configuration.configuration_api import ConfigurationAPI
from rest_client.AuthenticationRest import AuthenticationAPI

conf_obj = ConfigurationAPI()
conf_obj.load_data_from_ini()
AuthenticationAPI().createAutenthicationToken()

def getAllOfCouples():
    """
    | Get all couples from DB Inphinity

    :return: all couples from DB Inphinity
    :rtype: list 
    """
    list_couples = CoupleJson.getAllAPI()
    return list_couples

def getAllOfStrain():
    """
    | Get all Strains from DB Inphinity

    :return: all strains from DB Inphinity
    :rtype: list 
    """
    list_strain = StrainJson.getAllAPI()
    return list_strain
    
def getAllOfGenus():
    """
    | Get all Genus from DB Inphinity

    :return: all genus from DB Inphinity
    :rtype: list 
    """
    list_genus = GenusJson.getAllAPI()
    return list_genus

def getAllOfFamily():
    """
    | Get all families from DB Inphinity

    :return: all families from DB Inphinity
    :rtype: list 
    """
    list_family = FamilyJson.getAllAPI()
    return list_family

def getAllOfSourceData():
    """
    | Get all source data from DB Inphinity

    :return: all source data from DB Inphinity
    :rtype: list 
    """
    list_source_data = SourceDataJson.getAllAPI()
    return list_source_data

def getCouplesLysis(lysis_type):
    """
    | Get all couples of DB Inphinity according with the lysis type

    :param lysis_type: type of lysis

    :type lysis_type: int or int []

    :return: list of couples
    :rtype: list
    """
    # Check if lysis is an int 
    if isinstance(lysis_type, int):
        tmp = lysis_type
        # Change lysis_type in list
        lysis_type = []
        lysis_type.append(tmp)

    list_couples = getAllOfCouples()
    list_couples_lysis_type = []

    # Get all couple with lysis type
    for couple in list_couples:
        if couple.lysis in lysis_type:
            list_couples_lysis_type.append(couple)

    return list_couples_lysis_type