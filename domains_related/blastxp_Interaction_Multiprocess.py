#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 14:46:42 2017

@author: xavierbrochet
"""

# Programme qui blast les interactions Phage-Bact présent dans la base de données
# Sauvegarde des résultats du blast dans la base
# nouvelles featrures résultats des blasts

from Bio import SeqIO
import os, time, datetime
#import sys
#import pymysql.cursors
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from  Bio.Alphabet  import  generic_dna, generic_protein
import tempfile
from gestionBD import GestionBD
from multiprocessing import Pool

from objects_new.F_score_blast_P_new import *
from objects_new.PPI_couples_new import *

##pour les blastX décommenter les doubles ##

#bd = GestionBD("phage_bact", "xavier", "xavier", "localhost")
bd = GestionBD("inphinityDB_proj", "diogo_all", "Miguel1", "localhost")

def getInteraction(table_interaction):
    bd.connectDb()	
    bd.executerReq("SELECT id_couple_CP, FK_id_organism_bact_OR_CP, FK_id_organism_phage_OR_CP FROM "+table_interaction+" LIMIT 10")#+" LIMIT 10"
    interactionsL = bd.resultatReq()
    bd.close()
    return interactionsL

def blast_compute(table_interaction):

    try:
        interaction_id = str(table_interaction[0])
        bacterium_id = str(table_interaction[1])
        phage_id = str(table_interaction[2])
        
        bd.connectDb()
        bd.executerReq("SELECT id_protein_PT, sequence_PT FROM PROTEINS,GENES,ORGANISMS WHERE id_protein_PT=FK_id_protein_PT_GE and FK_id_organism_OR_GE=id_organism_OR AND id_organism_OR ="+bacterium_id) #id_protein_BD_online_PT,
        protBactL = bd.resultatReq()
        bd.close()
        bactProt_output_tmp = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
        #with open("bactProt_seq.fasta", "w") as output_handle:
        for prot in protBactL:
            idProt = prot[0]
            seq = prot[1]
            recordBacProt = SeqRecord(Seq(seq,generic_protein), name = "ProtBac", id = str(idProt), description = "") #bacterium_id , description = str(idProt)
            SeqIO.write(recordBacProt, bactProt_output_tmp, "fasta")
            bactProt_output_tmp.flush()
            #SeqIO.write(recordBacProt, output_handle, "fasta")
        
        #bactProt_output_tmp.seek(0)
        #print(bactProt_output_tmp.read())
        #bactProt_output_tmp.close()
                
        #################### pour le blastX!!! ####################
        #Phage genome.fna 
#        bd.connectDb()
#        bd.executerReq("select designation_GE, designation_SP, designation_ST, head_WD, sequence_WD from ORGANISMS, STRAINS, SPECIES, GENUSES, WHOLE_DNA WHERE FK_id_strain_ST_OR = id_strain_ST AND FK_id_specie_SP_ST = id_specie_SP AND id_dna_WD = FK_id_whole_DNA_DNA_OR AND id_organism_OR ="+phage_id)
#        genomePhage = bd.resultatReqOne()
#        bd.close()
#        if genomePhage != None:
#            #print("+ Phage_id: "+phage_id)
#            #print ("	Name: "+genomePhage[2])
#            seq = genomePhage[4]
#            if seq!='NA' :
#                namePhage = genomePhage[2]
#                head_description = genomePhage[3] 
#                recordPhage = SeqRecord(Seq(seq,generic_dna), name = namePhage, id = phage_id, description = head_description) #creation d'un objet SeqRecord (sequence fasta)
#                phage_output_tmp = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
#                SeqIO.write(recordPhage, phage_output_tmp, "fasta")
#                phage_output_tmp.flush()
##                with open(phage_id+"_phageSeq.fasta", "w") as output_handle2:
##                    SeqIO.write(recordPhage, output_handle2, "fasta")            
#                #phage_output_tmp.seek(0)
#                #print(phage_output_tmp.read())
#            else: #CONTIGS
#                bd.connectDb()
#                bd.executerReq("SELECT id_contig_CT, sequence_CT FROM ORGANISMS, WHOLE_DNA, CONTIGS WHERE FK_id_whole_DNA_DNA_OR = id_dna_WD AND FK_id_whole_genome_WD_CT = id_dna_WD AND id_organism_OR = "+phage_id)
#                contigPhageL = bd.resultatReq()
#                bd.close()
#                phage_output_tmp = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
#                for contig in contigPhageL:
#                    idContig = contig[0]
#                    seq = contig[1]
#                    recordPhageContig = SeqRecord(Seq(seq,generic_dna), name = "contigPhage", id = str(idContig), description = "")
#                    SeqIO.write(recordPhageContig, phage_output_tmp, "fasta")
#                    phage_output_tmp.flush()
        ################################################################################
            
        #Phage protein.faa (multi-fasta)
        bd.connectDb()
        bd.executerReq("SELECT id_protein_PT, sequence_PT FROM PROTEINS,GENES,ORGANISMS WHERE id_protein_PT=FK_id_protein_PT_GE and FK_id_organism_OR_GE=id_organism_OR AND id_organism_OR ="+phage_id) #id_protein_BD_online_PT,
        protPhageL = bd.resultatReq()
        bd.close()
        phageProt_output_tmp = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
#        with open(phage_id+"_phageProtSeq.fasta", "w") as output_handle3:
#            for prot in protPhageL:
#                idProt = prot[0]
#                seq = prot[1]
#                recordPhageProt = SeqRecord(Seq(seq,generic_protein), name = "ProtPhage", id = str(idProt), description = "") #bacterium_id , description = str(idProt)
#                SeqIO.write(recordPhageProt, output_handle3, "fasta")  
        
        for prot in protPhageL:
            idProt = prot[0]
            seq = prot[1]
            recordPhageProt = SeqRecord(Seq(seq,generic_protein), name = "ProtPhage", id = str(idProt), description = "") #bacterium_id , description = str(idProt)
            SeqIO.write(recordPhageProt, phageProt_output_tmp, "fasta")            
            phageProt_output_tmp.flush()
                #SeqIO.write(recordPhageProt, output_handle3, "fasta")
        
       # phageProt_output_tmp.seek(0)
       # print(phageProt_output_tmp.read())
        #print("***** BLAST *****")

        ##pour le test ne pas faire les blastX
        ##fileBlastxResult = "/home/diogo.leite/Scripts_xavier/blastScripts/"+interaction_id+"_blastx.txt"
        #fileBlastxResult2 = "/home/diogo.leite/Scripts_xavier/blastScripts/"+interaction_id+"_blastx2.txt"

        #os.system("blastx -subject bactProt_seq.fasta -query phageSeq.fasta -out "+fileBlastxResult2+" -evalue 0.01 -outfmt '6 std'")	
        #pour le test ne fait pas le blastX voir comment traiter les CONTIG et les GENOME COMPLET...
        ##os.system("blastx -subject "+str(bactProt_output_tmp.name)+" -query "+str(phage_output_tmp.name)+" -out "+fileBlastxResult+" -evalue 0.01 -outfmt '6 std'")
        #print("FIN du BLASTx")
        
        fileBlastpResult = "/home/diogo.leite/Scripts_xavier/blastScripts/"+interaction_id+"_blastp.txt"
        #print(("blastp -subject "+str(bactProt_output_tmp.name)+" -query "+str(phageProt_output_tmp.name)+" -out "+fileBlastpResult+" -evalue 0.01 -outfmt '6 std'"))
        os.system("blastp -subject "+str(bactProt_output_tmp.name)+" -query "+str(phageProt_output_tmp.name)+" -out "+fileBlastpResult+" -evalue 0.001 -outfmt '6 std'")
        #os.system("blastp -subject bactProt_seq.fasta -query phageProtSeq.fasta -out "+fileBlastpResult2+" -evalue 0.01 -outfmt '6 std'")
        bactProt_output_tmp.close()
        phageProt_output_tmp.close()
        
        #print("***** LOAD DATA BLAST *****")
##        bd.connectDb()
##        if os.path.isfile(fileBlastxResult):
##            with open(fileBlastxResult, "r") as fileP:
##                for ligne in fileP:
##                    ligne = ligne.rstrip()
##                    valueL = ligne.split("\t")                
                    #print("INSERT INTO F_SCORE_BLAST_X(FK_id_couple_CP_FSB_X, FK_id_protein_PT_FSB_X, pident_FSB_X, length_FSB_X, mismatch_FSB_X, gapopen_FSB_X, pstart_FSB_X, pend_FSB_X, bstart_FSB_X, bend_FSB_X, evalue_FSB_X, bitscore_FSB_X) VALUES ("+interaction_id+","+valueL[1]+","+valueL[2]+","+valueL[3]+","+valueL[4]+","+valueL[5]+","+valueL[6]+","+valueL[7]+","+valueL[8]+","+valueL[9]+","+valueL[10]+","+valueL[11]+")")
##                    bd.executerReq("INSERT F_SCORE_BLAST_X(FK_id_couple_CP_FSB_X, FK_id_protein_PT_FSB_X, pident_FSB_X, length_FSB_X, mismatch_FSB_X, gapopen_FSB_X, pstart_FSB_X, pend_FSB_X, bstart_FSB_X, bend_FSB_X, evalue_FSB_X, bitscore_FSB_X) VALUES ("+interaction_id+","+valueL[1]+","+valueL[2]+","+valueL[3]+","+valueL[4]+","+valueL[5]+","+valueL[6]+","+valueL[7]+","+valueL[8]+","+valueL[9]+","+valueL[10]+","+valueL[11]+")")
##                bd.commit()
##        bd.close()
##        os.remove(fileBlastxResult)
        
        if os.path.isfile(fileBlastpResult):
            list_scor_blast_p = []
            with open(fileBlastpResult, "r") as fileP:
                for ligne in fileP:
                    ligne = ligne.rstrip()
                    valueL = ligne.split("\t")

                    #add dans la table PPI_couples
                    PPI_couples_obj = PPI_couple(FK_prot_bact = valueL[1], FK_prot_phage = valueL[0], FK_couple = interaction_id)
                    ppi_couple_id   = PPI_couples_obj.create_ppi_couple_if_not_exist()

                    F_score_blast_P_obj = F_score_blast_P(id_f_score_blast_P = -1, pident = valueL[2], length = valueL[3], mismatch = valueL[4], gapopen = valueL[5], pstart = valueL[6], pend = valueL[7], bstart = valueL[8], bend = valueL[9], evalue = valueL[10], bitscore = valueL[11], plen = -1, blen = -1, FK_ppi_couple = ppi_couple_id)
                    list_scor_blast_p.append(F_score_blast_P_obj)
            for element in list_scor_blast_p:
                element.create_f_score_blast_p_no_verification()
        os.remove(fileBlastpResult)

    except Exception as e:
        print(e)
        print(table_interaction)

if __name__ == '__main__':
    start = time.time()
    listInteraction = getInteraction("COUPLES")
    #tList = len(listInteraction);
    
    with Pool(processes=None) as pool:
        pool.map(blast_compute, listInteraction)

    end = time.time()
    print(str(datetime.timedelta(seconds=int(end-start))))