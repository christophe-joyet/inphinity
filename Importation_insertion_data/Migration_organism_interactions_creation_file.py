import pandas as pd
from objects_new.Couples_new import Couple
from objects_new.Organisms_new import Organism
from objects_new.Strains_new import Strain
from objects_new.Species_new import Specie
from objects_new.Genus_new import Genus


def get_interactions_by_level_old_db(level_id):
    """
    return all the couple in the old db given the level id and ignore those obtain by Greg

    :return List with the couples
    :rtype List(Couple)

    """
    list_interact_level = []
    list_couples_old_db = Couple.get_all_couples()
    for interaction in list_couples_old_db:
        if interaction.fk_level_interact == 4 and interaction.fk_source_data != 3 and interaction.fk_source_data != 5:
            list_interact_level.append(interaction)

    return list_interact_level

def get_acc_by_organism(id_organism:int):
    """
    return all the couple in the old db given the level id and ignore those obtain by Greg

    :return List with the couples
    :rtype List(Couple)

    """
    organism_obj = Organism.get_organism_by_id(id_organism)
    acc_number = organism_obj.acc_num
    return acc_number

def get_phage_name_by_id(phage_id:int):
    """
    return the name of the phage according to it id

    :return phage name
    :rtype str

    """

    organism_obj = Organism.get_organism_by_id(phage_id)
    fk_strain = organism_obj.fk_strain

    strain_obj = Strain.get_strain_by_id(fk_strain)
    strain_designation = strain_obj.designation

    return strain_designation

def get_bacterium_taxonomy_designation_by_id_bacterium(bacterium_id:int):
    """
    return the strain, specie, and genus of the bacterium given it id

    :return list with the respective taxonomy
    :rtype list[str,str,str]

    """

    organism_obj = Organism.get_organism_by_id(bacterium_id)
    fk_strain = organism_obj.fk_strain


    strain_obj = Strain.get_strain_by_id(fk_strain)
    fk_specie = strain_obj.fk_specie
    strain_designation = strain_obj.designation

    specie_obj = Specie.get_specie_by_id(fk_specie)
    fk_genus = specie_obj.fk_genus
    specie_designation = specie_obj.designation

    genus_id = specie_obj.fk_genus
    genus_obj = Genus.get_genus_by_id(fk_genus)
    genus_designation = genus_obj.designation

    list_taxonomy = [strain_designation, specie_designation, genus_designation]

    return list_taxonomy


def get_details_data_interactions(list_interactions):

    data_list = []

    for interaction in list_interactions:
        id_bacterium = interaction.fk_bacteria
        id_phage = interaction.fk_phage

        acc_bacterium = get_acc_by_organism(id_bacterium)
        acc_phage = get_acc_by_organism(id_phage)

        id_interaction = interaction.id_couple
        interaction_type = interaction.interact_pn

        phage_name = get_phage_name_by_id(id_phage)
        bacterium_taxonomy = get_bacterium_taxonomy_designation_by_id_bacterium(id_bacterium)

        data_list.append([id_interaction, interaction_type, id_bacterium, acc_bacterium, bacterium_taxonomy[0], bacterium_taxonomy[1], bacterium_taxonomy[2], id_phage, acc_phage, phage_name])

    dataframe_data = pd.DataFrame(data_list, columns=['interaction_id','interaction_type','id_bacterium','acc_bacterium', 'bact_strain', 'bact_specie','bact_genus','id_phage','acc_phage','phage_name'])

    return dataframe_data

list_interactions = get_interactions_by_level_old_db(4)
print(len(list_interactions))
datafram_values = get_details_data_interactions(list_interactions)
print(datafram_values)

datafram_values.to_csv('data_interactions.csv',sep=',')
print(datafram_values)