from Bio import Entrez
import csv

from objects_new.Organisms_new import *
from objects_new.Strains_new import *
from objects_new.Species_new import *
from objects_new.Genus_new import *
from objects_new.Families_new import *


#Test

def get_name_organism(acc_num):
    try:
        Entrez.email = 'diogo1790@hotmail.com'
        handle = Entrez.efetch(db="nucleotide", id=acc_num, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        organims_name = record[0]["GBSeq_definition"]
    except:
        organims_name = 'Nothing error server or key'
    return organims_name



def write_file(array_data):
    with open("results_taxonomy.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(array_data)

#Test
aaa = get_name_organism('NC_010572.1')
print(aaa)



list_all_organisms = Organism.get_all_Organisms()

aux = 0
values_array = [['id_organism','Family', 'Genus', 'Specie', 'Strain', 'protein number', 'acc_number', 'name NCBI']]
for organism in list_all_organisms:
    acc_num = organism.acc_num
    qty_proteins = organism.qty_proteins
    strain_obj = Strain.get_strain_by_id(organism.fk_strain)
    specie_obj = Specie.get_specie_by_id(strain_obj.fk_specie)
    genus_obj = Genus.get_genus_by_id(specie_obj.fk_genus)
    family_obj = Family.get_family_by_id(genus_obj.fk_family)
    print("hello" + str(aux))
    aux += 1
    name_NCBI = get_name_organism(acc_num)
    array = [organism.id_organism, family_obj.designation, genus_obj.designation, specie_obj.designation, strain_obj.designation, qty_proteins, acc_num, get_name_organism(acc_num)]
    values_array.append(array)
    write_file(values_array)
print(len(list_all_organisms))


