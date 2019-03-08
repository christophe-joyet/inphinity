#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 14:46:42 2017

@author: xavierbrochet
"""

import csv
import codecs
import os
import time
import sys, datetime
import re
from multiprocessing import Pool

from objects_new.Families_new import *
from objects_new.Genus_new import *
from objects_new.Species_new import *
from objects_new.Strains_new import *

from objects_new.WholeDNA_new import *
from objects_new.Organisms_new import *
from objects_new.Contigs_new import *

from objects_new.Proteins_new import *
from objects_new.Gene_new import *

from files_treatment.fasta_parsing import *
from files_treatment.fasta_parsing2 import *

from objects_new.Couples_new import *

#dirData = '/Users/xavierbrochet/Documents/projets/INPHINITY/Database/TestloadDatabase'
dirData = '/home/diogo.leite/Scripts_xavier/data'
couplePhageBact = dirData+'/InteractionFound.csv'
dirBac = dirData+'/bacteria/'
dirPhage = dirData+'/phages2/'
GENOME = '_COMPLETE.fasta'
WGS    = '_WHOLESHOTGUN.fasta'
CDS    = '_DNA_CS.fasta'
PROT   = '_PROTEIN_CS.fasta'
level = {'strain': 4, 'specie': 3}

missingData = dirData+'/missingData.txt'

def parse_organism(taxon, dic_genome, dict_prot, dict_nn, source, orga, AC_Orga):
    
    #%% Taxonomy insertion insertion dans la base de données
    family_obj = Family(designation=taxon['family'])
    id_family = family_obj.create_family()
    genus_obj  = Genus(designation=taxon['genus'], fk_family=id_family)
    id_genus  = genus_obj.create_genus()
    specie_obj = Specie(designation=taxon['specie'], fk_genus=id_genus)
    id_specie = specie_obj.create_specie()
    strain_obj = Strain(designation=taxon['strain'], fk_specie=id_specie)
    id_strain = strain_obj.create_strain()

    #%% Whole_Genome
    if len(dic_genome) == 1:
        for genome_key in dic_genome.keys():
            whole_genome_obj = WholeDNA(head=str(dic_genome[genome_key][1]), head_id=AC_Orga, sequence=str(dic_genome[genome_key][2]))
            id_whole_genome = whole_genome_obj.create_whole_dna()
    else:
        compteur = 0;
        for genome_key in dic_genome.keys():
            if compteur == 0 :
                head = taxon['specie']+'_'+taxon['strain']
                whole_genome_obj = WholeDNA(head = head, head_id = AC_Orga, sequence="NA")
                id_whole_genome = whole_genome_obj.create_whole_dna()
            contig_obj = Contig(id_contig_db_outside = str(dic_genome[genome_key][0]), head=str(dic_genome[genome_key][1]), sequence = str(dic_genome[genome_key][2]), fk_id_whole_genome=id_whole_genome)
            id_contig = contig_obj.create_contig()
            compteur+=1
    
    #%% Organism     
    #source phagedb= 2; NCBI=1
    if len(dic_genome) == 1: #(assemble...)
        for genome_key in dic_genome.keys():
            organism_obj = Organism(id_organism=-1, gi='NA', acc_num=AC_Orga, qty_proteins=len(dict_prot.keys()), assembled=1, qty_contig=0, fk_source=1, fk_strain=id_strain, fk_type = orga, fk_whole_genome = id_whole_genome, fk_source_data=source)
            id_organism = organism_obj.create_organism()
    else:
        organism_obj = Organism(id_organism=-1, gi='NA', acc_num=AC_Orga, qty_proteins=len(dict_prot.keys()), assembled=0, qty_contig=len(dic_genome), fk_source =1, fk_strain=id_strain, fk_type=orga, fk_whole_genome=id_whole_genome, fk_source_data=source)
        id_organism = organism_obj.create_organism()

#%% Protein treatment    
    list_of_proteins = []

    aux = 0
    for protein_key in dict_prot.keys():
        if dict_nn is not None:
            try:
                seq_nc = dict_nn[protein_key.replace('prot','cds')]
                protein_obj = Protein(id_protein = -1, id_accession = dict_prot[protein_key][0],  designation = dict_prot[protein_key][4], sequence_prot = str(dict_prot[protein_key][5]), sequence_dna = str(seq_nc[5]), start_point = dict_prot[protein_key][1], end_point = dict_prot[protein_key][2], start_point_cnt = -1, end_point_cnt = -1, fk_id_contig = -1)
            except KeyError:
                protein_obj = Protein(id_protein = -1, id_accession = dict_prot[protein_key][0], designation = dict_prot[protein_key][4], sequence_prot = str(dict_prot[protein_key][5]), sequence_dna = "", start_point = dict_prot[protein_key][1], end_point = dict_prot[protein_key][2], start_point_cnt = -1, end_point_cnt = -1, fk_id_contig = -1)
        else:
            protein_obj = Protein(id_protein = -1, id_accession = dict_prot[protein_key][0], designation = dict_prot[protein_key][4], sequence_prot = str(dict_prot[protein_key][5]), sequence_dna = "", start_point = dict_prot[protein_key][1], end_point = dict_prot[protein_key][2], start_point_cnt = -1, end_point_cnt = -1, fk_id_contig = -1) #peut etreje ne sais pas trop les valeurs... ici a voir
        
        list_of_proteins.append(protein_obj)
    for protein in list_of_proteins:
        id_protein = protein.create_protein()
        gene_obj = Gene(FK_id_organism=id_organism, FK_id_protein=id_protein)
        value_id_gene = gene_obj.create_gene()
        aux += 1
        if aux % 1500 == 0:
            time.sleep(3)
            
    return id_organism

def setPhage(row, source):
    
    taxon = {}
    phageName = row[0]
    phage_acc = row[1]

    id_phage = Organism.get_organism_id_by_acc_or_designation(phage_acc, row[0].strip())
    if id_phage == -1:
            
        if row[1]!='NA':
            if source==2:
                phage = row[1]+"_"+row[0]
            elif source==1:
                phage = row[1]
        else:
            phage = row[0]

        if row[1]!='NA':
            if source==2:
                filePhageGen = dirPhage+row[1]+"_"+row[0]+GENOME        
                filePhageWGS = dirPhage+row[1]+"_"+row[0]+WGS
            elif source == 1:
                filePhageGen = dirPhage+row[1]+GENOME        
                filePhageWGS = dirPhage+row[1]+WGS
            filePhageCDS = dirPhage+row[1]+CDS
            filePhagePROT = dirPhage+row[1]+PROT
        else:
            filePhageGen = dirPhage+row[0]+GENOME
            filePhageWGS = dirPhage+row[0]+WGS
            filePhageCDS = dirPhage+row[0]+CDS
            filePhagePROT = dirPhage+row[0]+PROT
        
        if os.path.isfile(filePhageGen):
            fasta_file_genome = Fasta_parsing2(filePhageGen)
            dict_values_genome = fasta_file_genome.parse_fasta_genome()
        elif os.path.isfile(filePhageWGS):
            fasta_file_genome = Fasta_parsing2(filePhageWGS)
            dict_values_genome = fasta_file_genome.parse_fasta_WGS()
        else:
            with open(missingData,'a') as file:
                file.write(str(row)+': Phage -> manque le fichier GENOME')
            return None
        
        if os.path.isfile(filePhageCDS):
            fasta_file_cds = Fasta_parsing2(filePhageCDS)
            if row[1]!='NA':
                dict_values_cds = fasta_file_cds.parse_fasta_ncbi()
            else:
                dict_values_cds = fasta_file_cds.parse_fasta_GeneMark()
        else:
            dict_values_cds = None
        
        if os.path.isfile(filePhagePROT):
            fasta_file_protein = Fasta_parsing2(filePhagePROT)
            if row[1]!='NA':
                dict_values_proteins = fasta_file_protein.parse_fasta_ncbi()
            else:
                dict_values_proteins = fasta_file_protein.parse_fasta_GeneMark()
        else:
            with open(missingData,'a') as file:
                file.write(str(row)+': Phage -> manque le fichier PROTEIN')
            return None
        
        taxon = {'family': 'Phage no family', 'genus': 'Phage no genuse', 'specie': 'Phage no Specie', 'strain': row[0].strip()}
        id_phage = parse_organism(taxon, dict_values_genome, dict_values_proteins, dict_values_cds, source, 2, phage_acc)
    
    return id_phage
    
def SetBacterium(row, source):
    
    bacterium = row[0]
    
    id_bact = Organism.get_id_organism_by_acc(bacterium)
    if id_bact == -1:
        taxon = {}
        fileBactGen = dirBac+bacterium+GENOME
        fileBactWGS = dirBac+bacterium+WGS
        fileBactCDS = dirBac+bacterium+CDS
        fileBactPROT = dirBac+bacterium+PROT
        
        if os.path.isfile(fileBactGen):
            fasta_file_genome = Fasta_parsing2(fileBactGen)
            dict_values_genome = fasta_file_genome.parse_fasta_genome()
        elif os.path.isfile(fileBactWGS):
            fasta_file_genome = Fasta_parsing2(fileBactWGS)
            dict_values_genome = fasta_file_genome.parse_fasta_WGS()
        else:
            #ecrire dans un fichier... manque le fichier Genome ou WGS pour la bacterie truc
            with open(missingData,'a') as file:
                file.write(str(row[0])+': Bacterium -> manque le fichier Genome')
            return None
            
        if os.path.isfile(fileBactCDS):
            fasta_file_cds = Fasta_parsing2(fileBactCDS)
            dict_values_cds = fasta_file_cds.parse_fasta_ncbi()
        else:
            dict_values_cds = None
        
        if os.path.isfile(fileBactPROT):
            fasta_file_protein = Fasta_parsing2(fileBactPROT)
            with open(fileBactPROT,'r') as f:
                firstligne = f.readline()
                matchObj = re.search('^>gene_', firstligne)
                if matchObj:
                    dict_values_proteins = fasta_file_protein.parse_fasta_GeneMark()
                else:
                    dict_values_proteins = fasta_file_protein.parse_fasta_ncbi()
        else:
            with open(missingData,'a') as file:
                file.write(str(row[0])+': Bacterium -> manque le fichier PROTEIN')
            return None
        
        #preparation de la taxonomie...
        posL = len(row)
        strain= row[posL-1]
        if 'group' in row[posL-3]:
            specie = row[posL-4]+' '+row[posL-2]
            genus = row[posL-4]
            family = row[posL-5]
        else:
            specie = row[posL-3]+' '+row[posL-2]
            genus = row[posL-3]
            family = row[posL-4]
       
        taxon = {'family': family, 'genus': genus, 'specie': specie.strip(), 'strain': strain}
        id_bact = parse_organism(taxon, dict_values_genome, dict_values_proteins, dict_values_cds, source, 1, bacterium)
        
    return id_bact

#%% Parse fichier Couple
def parseCSV(row):

    print('######## COUPLE:########')
    print('row:'+str(row))
    if row[4]=='phagedb':
        source = 2
    else:
        source = 1
    id_phage = setPhage(row[0:3], source)
    
    id_bact = SetBacterium(row[3:len(row)], 1)
    typeOfInteraction = row[2]
    
    
    #si une des deux valeur son None pas de fichier pour cet organism alors il faut écrire dans un fichier de valeur manquantes et ne pas charger le couple....
    if id_phage is not None and id_bact is not None:    
        ####### Couples ######
        couple_obj = Couple(id_couple = -1, interact_pn = 1, fk_bacteria=id_bact, fk_phage=id_phage, fk_type_inter=2, fk_level_interact=level[row[2]], fk_data_source=source)
        id_couple = couple_obj.create_couple()
    else:
        with open(missingData,'a') as file:
                file.write('---> PB (ci-dessus) avec le couple couple: '+str(row))

#%% Main
if __name__ == '__main__':
    start = time.time()

    file = open(couplePhageBact, "rb")
    decode = codecs.getreader('utf-8')
    cr = csv.reader(decode(file), delimiter=';')
    coupleL = list(cr)
    for val in coupleL:
        parseCSV(val)
    
    #multiprocessing 
#    with Pool(processes=None) as pool:
#        pool.map(parseCSV, coupleL)
    
    end = time.time()
    print(str(datetime.timedelta(seconds=int(end-start))))
