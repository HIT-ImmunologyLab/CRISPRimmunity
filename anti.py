#-*-coding:utf-8-*-

from Bio import SeqIO
import os,subprocess
import sys
import json
import re
import threading
from collections import OrderedDict
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import Entrez
import xml.etree.ElementTree as ET
import re,argparse
import os
from collections import Counter
import traceback
import queue
import requests
from collections import OrderedDict
import copy
import http
import math
import argparse,time
from Bio.Blast import NCBIWWW
import time
from xml.etree.ElementTree import iterparse
from lxml import etree
import shutil
import pandas as pd
import numpy as np
from  scipy.stats import chi2_contingency
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import Levenshtein

# queue_prot_filed = queue.Queue()
# queue_bac_filed = queue.Queue()
# queue_from_prot_to_bac = queue.Queue()

def merge_intervals(intervals):
    """
    :type intervals: List[Interval]
    :rtype: List[Interval]
    """
    length = 0
    intervals_sorted = sorted(intervals, key=lambda x : x[0])
    result = []
    for interval in intervals_sorted:
        # check overlap
        if result and result[-1][1] >= interval[0]:
            # merge
            result[-1][1] = max(result[-1][1], interval[1])
        else:
            result.append(interval)
    #return result
    for item in result:
        length += abs(item[1]-item[0]) + 1
    return length

def predict_known_anti_hth(protein_file,outdir,identity,coverage,prog_dir):
    script = "%s/predict_known_anti.py"%(prog_dir)
    command = "/home/qiusuo/miniconda3/bin/python %s --protein %s --output %s --identity %s --coverage %s --prog_dir %s"%(script,protein_file,outdir,str(identity),str(coverage),prog_dir)
    os.system(command)


def prophage_find_hth(prophage_file,faa_file,faa_type,pro_in_phage_faa_file,pro_in_phage_info_file):
    prophage_dic = {}
    with open(prophage_file) as fi:
        lines = fi.readlines()
        if len(lines) <= 1:
            print('No prophage containing self-targeting found !')
        else:
            for line in lines[1:]:
                if line.strip()=='':
                    continue
                prophage_info = line.strip().split('\t')
                prophage_loc = [prophage_info[4],prophage_info[5]]
                contig_id = prophage_info[1]
                if contig_id not in prophage_dic.keys():
                    prophage_dic.update({contig_id:[]})
                prophage_dic[contig_id].append(prophage_info)
    
    protein_dic = {}
    with open(faa_file) as fi:
        lines = fi.readlines()
        for index,line in enumerate(lines):
            if line.startswith('>'):
                contig_id = line.split('.')[0].replace('>','')
                if contig_id in prophage_dic.keys() and lines[index+1] != 'unknown':
                    if contig_id in protein_dic.keys():
                        protein_dic[contig_id].append(line)
                        protein_dic[contig_id].append(lines[index+1])
                    else:
                        protein_dic[contig_id] = [line,lines[index+1]]
    
    faa_ou = open(pro_in_phage_faa_file,'w')
    f_prophage_protein_info = open(pro_in_phage_info_file,'w')
    f_prophage_protein_info.write('protein_id\tprotein_location\tprotein_size\tprophage_contig_id\tprophage_start\tprophage_end\tmethod\tkey_proteins\tbest_hit_species\tCDS_number\tattl_region\tattr_region\n')
    for contig_id in protein_dic.keys():
        for index,line in enumerate(protein_dic[contig_id]):
            if line.startswith('>'):
                if '|' in line:
                    if line.count('|')>=3:
                        faa_type = 'genbank'
                    else:
                        faa_type = 'fasta'
                else:
                    faa_type = 'fasta'
                if faa_type == 'fasta':
                    protein_loc = line.strip().split('_')
                    protein_start = int(protein_loc[-3])
                    protein_end = int(protein_loc[-2])
                    direction = protein_loc[-1]
                elif faa_type == 'genbank':
                    protein_loc = line.strip().split('|')[2].split('_')
                    protein_start = int(protein_loc[0])
                    protein_end = int(protein_loc[1])
                    direction = protein_loc[-1]
                protein_size = len(protein_dic[contig_id][index+1].strip())
                protein_id = line.strip().strip('>')
                if contig_id in prophage_dic.keys():
                    for prophage_info in prophage_dic[contig_id]:
                        phage_start = int(prophage_info[4])
                        phage_end = int(prophage_info[5])
                        if protein_start >= phage_start and protein_end <= phage_end:
                            faa_ou.write(line+protein_dic[contig_id][index+1])
                            f_prophage_protein_info.write(protein_id+'\t'+str(protein_start)+'_'+str(protein_end)+'_'+direction+'\t'+str(protein_size)+'\t'+contig_id+'\t'+'\t'.join(prophage_info[4:])+'\n')
                            f_prophage_protein_info.flush()
                            break                
    faa_ou.close()
    f_prophage_protein_info.close()

def ParseResultrpsblast(inFileName,outFileName):
    text = open(inFileName).read()
    text = re.sub(u"[\x00-\x08\x0b-\x0c\x0e-\x1f]+", u" ", text)
    #print(text)
    root = ET.fromstring(text)
    BlastOutput_iterations = root.find("BlastOutput_iterations")
    f = open(outFileName,'w')
    f.write('query_id\thit_id\thit_description\tquery_from\tquery_to\thit_from\thit_to\te_value\tidentity\thit_coverage\tquery_coverage\n')
    for Iteration in BlastOutput_iterations.findall("Iteration"):
        strTemp = str(Iteration.find("Iteration_query-def").text)
        Iteration_hits = Iteration.find("Iteration_hits")
        query_length = float(Iteration.find("Iteration_query-len").text)
        for Hit in Iteration_hits.findall("Hit"):
            strDef = str(Hit.find("Hit_def").text)
            Hit_len =  Hit.find("Hit_len").text
            Hit_ID = Hit.find("Hit_id").text
            Hsp = Hit.find("Hit_hsps").find("Hsp")
            query_from = Hsp.find("Hsp_query-from").text
            query_to = Hsp.find("Hsp_query-to").text
            Hit_from = Hsp.find("Hsp_hit-from").text
            Hit_to = Hsp.find("Hsp_hit-to").text
            Hsp_evalue = Hsp.find("Hsp_evalue").text
            Hsp_score = Hsp.find("Hsp_score").text
            hit_length = abs(int(Hit_to)-int(Hit_from))+1
            query_len = abs(int(query_to)-int(query_from))+1
            hit_coverage = str(float(hit_length)/ float(Hit_len))
            query_coverage = str(float(query_len)/ query_length)
            #identity = str(Hsp.find("Hsp_identity").text)
            identity = str(float(Hsp.find("Hsp_identity").text) / float(Hsp.find("Hsp_align-len").text))
            f.write(strTemp+'\t'+Hit_ID+'\t'+strDef+'\t'+query_from+'\t'+query_to+'\t'+Hit_from+'\t'+Hit_to+'\t'+Hsp_evalue+'\t'+identity+'\t'+hit_coverage+'\t'+query_coverage+'\n')
            f.flush()
    f.close()


def get_pro_up_hth(faa_file,faa_type,hth_file,prophage_file,prophage_protein_file,prophage_protein_info_file,up_faa_file,info_file,up_num=3,prot_size=400):
    info_ou = open(info_file,'w')
    info_ou.write('protein_id\tprotein_location\tprotein_size\tin_prophage\tdistance_hth\tneighbor_hth_id\tneighbor_hth_domain\tprophage_contig_id\tprophage_start\tprophage_end\tmethod\tkey_proteins\tbest_hit_species\tCDS_number\tattl_region\tattr_region\n')
    
    header_arr_dict = {}
    seq_arr_dict = {}
    
    with open(faa_file) as fi:
        lines = fi.readlines()
        for index,line in enumerate(lines):
            if line.startswith('>') and lines[index+1] != 'unknown':
                contig_id = '_'.join(line.split('_')[0:-3]).strip('>').strip()
                if contig_id not in header_arr_dict.keys():
                    header_arr_dict.update({contig_id:[]})
                header_arr_dict[contig_id].append(line.strip().replace('>',''))
                seq_arr_dict[contig_id].append(lines[index+1].strip())
    
    prophage_protein_info_dict = {}
    with open(prophage_protein_info_file) as f:
        prophage_contents = f.readlines()
    for prophage_pro_info in prophage_contents[1:]:
        prophage_pro_info = prophage_pro_info.strip().split('\t')
        protein_id = prophage_pro_info[0]
        prophage_protein_info_dict.update({protein_id:prophage_pro_info})      
    
    ou = open(up_faa_file,'w')
    write_protein_ids = []
    with open(hth_file) as fi:
        lines = fi.readlines()
        #info_ou.write(lines[0].strip()+'\tup_protein_info\n')
        for line in lines[1:]:
            #info_ou.write(line.strip())           
            arr = line.strip().split('\t')
            hth_protein_id = arr[0]
            hit_contig_id = '_'.join(hth_protein_id.split('_')[0:-3])
            hth_domain_id = arr[1]
            hth_domain_def = arr[2].split('.')[0].strip()
            if hth_protein_id not in prophage_protein_info_dict.keys():
                continue
            protein_prophage_info = prophage_protein_info_dict[protein_id]
            if '|' in arr[0]:
                if arr[0].count('|')>=3:
                    faa_type = 'genbank'
                else:
                    faa_type = 'fasta'
            else:
                faa_type = 'fasta'
            if faa_type == 'fasta':
                query_info = arr[0].split('_')
                query_contig = '_'.join(query_info[0:-3])
                query_start = int(query_info[-3])
                query_end = int(query_info[-2])
                strand = query_info[-1]
            else:
                query_info = arr[0].split('|')
                strand = query_info[2].split('_')[-1]
                query_contig = query_info[0].replace('>','')
                query_start = query_info[2].split('_')[0]
                query_end =  query_info[2].split('_')[1]
            
            hth_index = header_arr_dict[hit_contig_id].index(arr[0])
            pro_str = ''
            if strand == '+':
                up_3 = 0
                pro_index = hth_index
                while up_3 < int(up_num):
                    pro_index = pro_index - 1
                    if pro_index < 0:
                        break
                    else:
                        pro_contig = ''
                        if '|' in header_arr_dict[hit_contig_id][pro_index]:
                            if header_arr_dict[hit_contig_id][pro_index].count('|')>=3:
                                faa_type = 'genbank'
                            else:
                                faa_type = 'fasta'
                        else:
                            faa_type = 'fasta'
                        if faa_type == 'fasta' and header_arr_dict[hit_contig_id][pro_index].endswith('+'):
                            pro_info = header_arr_dict[hit_contig_id][pro_index].split('_')
                            pro_contig = '_'.join(pro_info[0:-3])
                            pro_start = int(pro_info[-3])
                            pro_end = int(pro_info[-2])
                            pro_strand = pro_info[-1]
                        elif faa_type == 'genbank' and header_arr_dict[hit_contig_id][pro_index].split('|')[2].split('_')[-1] == '+':
                            pro_info = header_arr_dict[hit_contig_id][pro_index].split('|')
                            pro_contig = pro_info[0].replace('>','')
                            pro_start = int(pro_info[2].split('_')[0])
                            pro_end = int(pro_info[2].split('_')[1])
                            pro_strand = pro_info[2].split('_')[-1]                       
                        
                        if pro_contig == query_contig and pro_strand == '+':
                            up_3 += 1
                            if len(seq_arr[pro_index]) < int(prot_size):
                                ou.write('>'+header_arr_dict[hit_contig_id][pro_index]+'\n'+seq_arr_dict[hit_contig_id][pro_index]+'\n')
                                if pro_str == '':
                                    pro_str = header_arr_dict[hit_contig_id][pro_index]+'_'+str(len(seq_arr_dict[hit_contig_id][pro_index].strip()))
                                else:
                                    pro_str = pro_str + '::' + header_arr_dict[hit_contig_id][pro_index]+'_'+str(len(seq_arr_dict[hit_contig_id][pro_index].strip()))                        
                                c_att_pro_sizr = len(seq_arr_dict[hit_contig_id][pro_index])
                                if header_arr_dict[hit_contig_id][pro_index] in prophage_protein_info_dict.keys():
                                    in_prophage_flag = 'yes'
                                else:
                                    in_prophage_flag = 'no'
                                if header_arr_dict[hit_contig_id][pro_index] not in write_protein_ids:
                                    write_protein_ids.append(header_arr_dict[hit_contig_id][pro_index])
                                    info_ou.write(header_arr_dict[hit_contig_id][pro_index]+'\t'+str(pro_start)+'_'+str(pro_end)+'_'+pro_strand+'\t'+str(c_att_pro_sizr)+'\t'+in_prophage_flag+'\t'+str(hth_index-pro_index)+'\t'+hth_protein_id+'\t'+hth_domain_def+'\t'+'\t'.join(prophage_protein_info_dict[hth_protein_id][3:])+'\n')
                                    info_ou.flush()
            else:
                up_3 = 0
                pro_index = hth_index
                while up_3 < int(up_num):
                    pro_index = pro_index + 1
                    if pro_index >= len(header_arr_dict[hit_contig_id]):
                        break
                    else:
                        pro_contig = ''
                        if '|' in header_arr_dict[hit_contig_id][pro_index]:
                            if header_arr_dict[hit_contig_id][pro_index].count('|')>=3:
                                faa_type = 'genbank'
                            else:
                                faa_type = 'fasta'
                        else:
                            faa_type = 'fasta'
                        if faa_type == 'fasta' and header_arr_dict[hit_contig_id][pro_index].endswith('-'):
                            pro_info = header_arr_dict[hit_contig_id][pro_index].split('_')
                            pro_contig = '_'.join(pro_info[0:-3])
                            pro_start = int(pro_info[-3])
                            pro_end = int(pro_info[-2])
                            pro_strand = pro_info[-1]
                        elif faa_type == 'genbank' and header_arr_dict[hit_contig_id][pro_index].split('|')[2].split('_')[-1] == '-':
                            pro_info = header_arr_dict[hit_contig_id][pro_index].split('|')
                            pro_contig = pro_info[0].replace('>','')
                            pro_start = int(pro_info[2].split('_')[0])
                            pro_end = int(pro_info[2].split('_')[1])
                            pro_strand = pro_info[2].split('_')[-1]
                        if pro_contig == query_contig and pro_strand == '-':
                            up_3 += 1
                            if len(seq_arr_dict[hit_contig_id][pro_index]) < int(prot_size):
                                ou.write('>'+header_arr_dict[hit_contig_id][pro_index]+'\n'+seq_arr_dict[hit_contig_id][pro_index]+'\n')
                                if pro_str == '':
                                    pro_str = header_arr_dict[hit_contig_id][pro_index]+'_'+str(len(seq_arr_dict[hit_contig_id][pro_index].strip()))
                                else:
                                    pro_str = pro_str + '::' + header_arr_dict[hit_contig_id][pro_index]+'_'+str(len(seq_arr_dict[hit_contig_id][pro_index].strip()))
                                
                                c_att_pro_sizr = len(seq_arr_dict[hit_contig_id][pro_index])
                                if header_arr_dict[hit_contig_id][pro_index] in prophage_protein_info_dict.keys():
                                    in_prophage_flag = 'yes'
                                else:
                                    in_prophage_flag = 'no'
                                if header_arr_dict[hit_contig_id][pro_index] not in write_protein_ids:
                                    write_protein_ids.append(header_arr_dict[hit_contig_id][pro_index])
                                    info_ou.write(header_arr_dict[hit_contig_id][pro_index]+'\t'+str(pro_start)+'_'+str(pro_end)+'_'+pro_strand+'\t'+str(c_att_pro_sizr)+'\t'+in_prophage_flag+'\t'+str(pro_index-hth_index)+'\t'+hth_protein_id+'\t'+hth_domain_def+'\t'+'\t'.join(prophage_protein_info_dict[hth_protein_id][3:])+'\n')
                                    info_ou.flush()
            
    ou.close()
    info_ou.close()


def rpsblastp_protein(pro_up_hth_faa,out_file_prefix,prog_dir):
    pro_up_domain_xml = out_file_prefix+'_domain.xml'
    pro_up_domain_parse = out_file_prefix+'_domain.txt'
    blast_cline = "%s/software/rpsblast -query %s -comp_based_stats 0 -evalue 0.01 -seg no -max_target_seqs 1 -outfmt 5 -num_threads 20 -db %s/database/CDD_rps/Cdd -out %s" % (prog_dir,pro_up_hth_faa,prog_dir,pro_up_domain_xml)
    os.system(blast_cline)
    ParseResultrpsblast(pro_up_domain_xml,pro_up_domain_parse)


def pro_no_func(faa_file,up_faa_file,up_protein_info_file,pro_up_domain_parse,prophage_protein_info_file,anti_align_file,candi_anti_prefix):
    candi_anti_protein_file = candi_anti_prefix+'_protein.faa'
    f_protein = open(candi_anti_protein_file,'w')
    
    candi_anti_protein_info_file = candi_anti_prefix+'_protein_info.txt'
    f_info = open(candi_anti_protein_info_file,'w')
    f_info.write('protein_id\tprotein_location\tprotein_size\tin_prophage\tdistance_hth\tneighbor_hth_id\tneighbor_hth_domain\tprophage_contig_id\tprophage_start\tprophage_end\tmethod\tkey_proteins\tbest_hit_species\tCDS_number\tattl_region\tattr_region\tknown_anti\tidentity\tquery_coverage\thit_coverage\n')
    func_dic = {}
    with open(pro_up_domain_parse) as fi:
        lines = fi.readlines()
        for line in lines[1:]:
            arr = line.split('\t')
            protein_id = arr[0]
            if protein_id not in func_dic.keys():
                func_dic.update({protein_id:arr})
    
    #known anti
    anti_dic = {}
    with open(anti_align_file) as fi:
        lines = fi.readlines()
        for line in lines[1:]:
            arr = line.strip().split('\t')
            protein_id = arr[0]
            if protein_id not in anti_dic.keys():
                anti_dic.update({protein_id:['','','','']})
            write_info = [arr[1]]+arr[-3:]
            for index,item in enumerate(anti_dic[protein_id]):
                anti_dic[protein_id][index] = anti_dic[protein_id][index]+','+write_info[index]
            for index,item in enumerate(anti_dic[protein_id]):
                anti_dic[protein_id][index] = anti_dic[protein_id][index].lstrip(',')
    
    ##up protein prophage info file
    up_protein_info_dict = {}
    with open(up_protein_info_file) as f:
        contents = f.readlines()
    for line in contents[1:]:
        line = line.strip().split('\t')
        protein_id = line[0]
        up_protein_info_dict.update({protein_id:line})
    
    #predicted anti proteins
    predict_protein_ids = []
    with open(up_faa_file) as fi:
        lines = fi.readlines()
        for index,line in enumerate(lines):
            if line.startswith('>'):
                protein_id = line.strip().replace('>','')
                if protein_id in func_dic.keys():
                    protein_domain_info = func_dic[protein_id]
                if protein_id in anti_dic.keys():
                    protein_anti_info =  anti_dic[protein_id]
                protein_prophage_info = up_protein_info_dict[protein_id]
                if protein_id in anti_dic.keys():
                    predict_protein_ids.append(protein_id)
                    f_info.write('\t'.join(up_protein_info_dict[protein_id]).strip()+'\t'+'\t'.join(anti_dic[protein_id])+'\n')
                    f_info.flush()
                    f_protein.write(line+lines[index+1])
                else:
                    if protein_id not in func_dic.keys():
                        predict_protein_ids.append(protein_id)
                        f_info.write('\t'.join(up_protein_info_dict[protein_id]).strip()+'\t'+'\t'.join(['NA']*4)+'\n')
                        f_info.flush()
                        f_protein.write(line+lines[index+1])
    
    #prophage protein info
    prophage_protein_info_dict = {}
    with open(prophage_protein_info_file) as f:
        prophage_contents = f.readlines()
    for prophage_pro_info in prophage_contents[1:]:
        prophage_pro_info = prophage_pro_info.strip().split('\t')
        protein_id = prophage_pro_info[0]
        prophage_protein_info_dict.update({protein_id:prophage_pro_info})  
    
    #all protein file
    with open(faa_file) as f:
        faa_contents = f.read().strip()
    protein_sequence_dict = {}
    for pro in faa_contents.split('>')[1:]:
        pro_id = pro.split('\n')[0].strip()
        protein_sequence_dict.update({pro_id:pro})
    
    #predicted known anti-proteins
    for protein_id in anti_dic.keys():
        protein_sequence = ''.join(protein_sequence_dict[protein_id].split('\n')[1:]).strip()
        protein_size = len(protein_sequence)        
        faa_type = 'fasta'    
        if protein_id not in predict_protein_ids:
            if protein_id in prophage_protein_info_dict.keys():
                protein_inf_prophage = 'yes'
                prophage_info = prophage_protein_info_dict[protein_id][3:]
            else:
                protein_inf_prophage = 'no'
                prophage_info = ['NA']*9
            if '|' in protein_id:
                if protein_id.count('|')>=3:
                    faa_type = 'genbank'
                else:
                    faa_type = 'fasta'
            else:
                faa_type = 'fasta'
            if faa_type == 'fasta':
                pro_info = protein_id.split('_')
                pro_contig = '_'.join(pro_info[0:-3])
                pro_start = int(pro_info[-3])
                pro_end = int(pro_info[-2])
                pro_strand = pro_info[-1]
           
            elif faa_type == 'genbank':
                pro_info = protein_id.split('|')
                pro_contig = pro_info[0].replace('>','')
                pro_start = int(pro_info[2].split('_')[0])
                pro_end = int(pro_info[2].split('_')[1])
                pro_strand = pro_info[2].split('_')[-1]

            f_info.write(protein_id+'\t'+str(pro_start)+'_'+str(pro_end)+'_'+pro_strand+'\t'+str(protein_size)+'\t'+protein_inf_prophage+'\tNA\tNA\tNA'+'\t'+'\t'.join(prophage_info).strip()+'\t'+'\t'.join(anti_dic[protein_id])+'\n')
            f_info.flush()
            f_protein.write('>'+protein_id+'\n'+protein_sequence.strip()+'\n')
    f_protein.close()
    f_info.close()

def add_acr_homology_to_candidate_dict(known_anti_file,acr_protein_candidate_dict):
    with open(known_anti_file,'r')as fin:
        lines = fin.readlines()
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_proten_id = content[0]
            if cur_proten_id not in acr_protein_candidate_dict:
                acr_protein_candidate_dict[cur_proten_id] = content[1:]

def get_upstream_or_downstream_prot_file(known_aca_file,prot_key_index_value_dict, index_key_prot_info_value_dict,aca_neighbor_prot_seq_file,aca_neighbor_prot_info_file,distance_prot_num=3,acr_prot_size=400):
    fout = open(aca_neighbor_prot_seq_file,'w')
    finfo = open(aca_neighbor_prot_info_file,'w')

    with open(known_aca_file,'r')as fin:
        lines = fin.readlines()
        header_list = ['neighbor_prote_id','distance_with_aca']+lines[0].strip('\n').split('\t')
        finfo.write('\t'.join(header_list)+'\n')
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_aca_protein_info = content[0]
            cur_contig_id = '_'.join(cur_aca_protein_info.split('_')[0:-3])
            cur_aca_protein_direction = cur_aca_protein_info.split('_')[-1]
            cur_aca_index = prot_key_index_value_dict[cur_contig_id][cur_aca_protein_info]
            if cur_aca_protein_direction == '+':
                for i in range(max(int(cur_aca_index)-int(distance_prot_num),1),int(cur_aca_index)):
                    seq_len = len(index_key_prot_info_value_dict[cur_contig_id][str(i)].split('\n')[1])
                    title_info = index_key_prot_info_value_dict[cur_contig_id][str(i)].split('\n')[0].replace('>','')
                    if int(seq_len) <= int(acr_prot_size):
                        fout.write(index_key_prot_info_value_dict[cur_contig_id][str(i)])
                        distance_with_aca = abs(int(i)-int(cur_aca_index))-1
                        neighbor_prot_info_list = [title_info,str(distance_with_aca)]+content
                        finfo.write('\t'.join(neighbor_prot_info_list)+'\n')
            elif cur_aca_protein_direction == '-':
                for i in range(int(cur_aca_index)+1,min(int(cur_aca_index)+int(distance_prot_num)+1,len(prot_key_index_value_dict[cur_contig_id]))):
                    seq_len = len(index_key_prot_info_value_dict[cur_contig_id][str(i)].split('\n')[1])
                    title_info = index_key_prot_info_value_dict[cur_contig_id][str(i)].split('\n')[0].replace('>', '')
                    if int(seq_len) <= int(acr_prot_size):
                        fout.write(index_key_prot_info_value_dict[cur_contig_id][str(i)])
                        distance_with_aca = abs(int(i) - int(cur_aca_index))-1
                        neighbor_prot_info_list = [title_info,str(distance_with_aca)] + content
                        finfo.write('\t'.join(neighbor_prot_info_list) + '\n')
    fout.close()
    finfo.close()

def get_prote_domain_dict(prot_domain_file,prot_domain_dict):
    with open(prot_domain_file,'r')as fin:
        lines = fin.readlines()
        lines[0].strip('\n').split('\t')[1:]
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            if content[0] not in prot_domain_dict:
                prot_domain_dict[content[0]] = content[1:]

def predict_novel_acr_by_aca(known_aca_file,prot_key_index_value_dict,index_key_prot_info_value_dict,outdir,prog_dir,distance_prot_num=3,acr_prot_size=400):
    # get upstream or downstream protein of the Aca protein according to the Aca protein direction
    # filter the protein size
    aca_neighbor_prot_seq_file = os.path.join(outdir, 'aca_neighbor_%s_proteins_filter_prot_size_%s_aa_seq.faa' % (str(distance_prot_num),str(acr_prot_size)))
    aca_neighbor_prot_info_file = os.path.join(outdir, 'aca_neighbor_%s_proteins_filter_prot_size_%s_aa_info.txt' % (str(distance_prot_num),str(acr_prot_size)))
    get_upstream_or_downstream_prot_file(known_aca_file,prot_key_index_value_dict, index_key_prot_info_value_dict,aca_neighbor_prot_seq_file,aca_neighbor_prot_info_file,distance_prot_num,acr_prot_size)

    if os.path.getsize(aca_neighbor_prot_seq_file) > 0:
        # annotate the domain of the upstream or downstream proteins
        up_protein_domain_prefix = os.path.join(outdir, 'aca_neighbor_%s_filter_prot_size_%s_aa_rpsblast' % (str(distance_prot_num),str(acr_prot_size)))
        rpsblastp_protein(aca_neighbor_prot_seq_file, up_protein_domain_prefix,prog_dir)

        # get the Acr candidate dict by Aca
        prot_domain_file = '%s_domain.txt'%(up_protein_domain_prefix)
        prot_domain_dict = {}
        get_prote_domain_dict(prot_domain_file, prot_domain_dict)

        acr_candidate_by_aca_file = os.path.join(outdir, 'aca_neighbor_%s_filter_prot_size_%s_aa.txt' % (str(distance_prot_num),str(acr_prot_size)))
        fout = open(acr_candidate_by_aca_file,'w')
        with open(aca_neighbor_prot_info_file,'r')as fin:
            lines = fin.readlines()
            header_list = lines[0].strip('\n').split('\t')
            fout.write('\t'.join(header_list)+'\n')
            for line in lines[1:]:
                content = line.strip('\n').split('\t')
                cur_prot_id = content[0]
                if cur_prot_id not in prot_domain_dict:
                    strwrite_list = content
                    fout.write('\t'.join(strwrite_list)+'\n')
        fout.close()
    else:
        acr_candidate_by_aca_file = os.path.join(outdir, 'aca_neighbor_%s_filter_prot_size_%s_aa.txt' % (
        str(distance_prot_num), str(acr_prot_size)))
        fout = open(acr_candidate_by_aca_file, 'w')
        with open(aca_neighbor_prot_info_file,'r')as fin:
            lines = fin.readlines()
            header_list = lines[0].strip('\n').split('\t')
            fout.write('\t'.join(header_list)+'\n')
        fout.close()

def get_prote_dict(protein_file,prot_key_index_value_dict,index_key_prot_info_value_dict):
    with open(protein_file,'r')as fin:
        content = fin.read()
        elems = content.split('>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            cur_title_line = elem.split('\n')[0]
            cur_contig_id = '_'.join(cur_title_line.split('_')[0:-3])
            if cur_contig_id not in prot_key_index_value_dict:
                prot_key_index_value_dict[cur_contig_id] = {}
            if cur_title_line not in prot_key_index_value_dict[cur_contig_id]:
                cur_index = str(len(prot_key_index_value_dict[cur_contig_id])+1)
                prot_key_index_value_dict[cur_contig_id][cur_title_line] = cur_index
            if cur_contig_id not in index_key_prot_info_value_dict:
                index_key_prot_info_value_dict[cur_contig_id] = {}
            if cur_index not in index_key_prot_info_value_dict[cur_contig_id]:
                index_key_prot_info_value_dict[cur_contig_id][cur_index] = '>'+elem

def get_protospacer_location_dict(st_info,protospacer_location_dict):
    ori_protospacer_location_dict = {}
    with open(st_info,'r')as fin:
        lines = fin.readlines()
        header_list = lines[0].strip('\n').split('\t')
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_protospacer_info_list = content[header_list.index("protospacer_location")].split('|')
            for cur_protospacer_info_item in cur_protospacer_info_list:
                cur_protospacer_info_content = cur_protospacer_info_item.split('_')
                if cur_protospacer_info_item=='NA':
                    continue
                
                cur_protospacer_contig_id = '_'.join(cur_protospacer_info_content[0:-1]).split('.')[0]
                cur_protospacer_loc_start = cur_protospacer_info_content[-1].split('-')[0]
                cur_protospacer_loc_end = cur_protospacer_info_content[-1].split('-')[1]
                if cur_protospacer_contig_id not in ori_protospacer_location_dict:
                    ori_protospacer_location_dict[cur_protospacer_contig_id] = {}
                if cur_protospacer_loc_start not in ori_protospacer_location_dict[cur_protospacer_contig_id]:
                    ori_protospacer_location_dict[cur_protospacer_contig_id][cur_protospacer_loc_start] = '%s-%s'%(cur_protospacer_loc_start,cur_protospacer_loc_end)

                for contig_key in ori_protospacer_location_dict:
                    protospacer_location_dict[contig_key] = {}
                    sorted_key_list = sorted(ori_protospacer_location_dict[contig_key])
                    for sorted_key in sorted_key_list:
                        protospacer_location_dict[contig_key][sorted_key] = ori_protospacer_location_dict[contig_key][sorted_key]

def check_prophage_contain_protospacer(protospacer_location_dict,contig_id,prophage_start,prophage_end):
    prophage_contain_protospacer_flag = 'No'
    cur_contig_protospacer_dict = protospacer_location_dict[contig_id]
    for protospacer_start_key in cur_contig_protospacer_dict:
        if int(prophage_end) < int(protospacer_start_key):
            break
        else:
            protospacer_start_loc = cur_contig_protospacer_dict[protospacer_start_key].split('-')[0]
            protospacer_end_loc = cur_contig_protospacer_dict[protospacer_start_key].split('-')[1]
            if int(protospacer_start_loc) >= int(prophage_start) and int(protospacer_end_loc) <= int(
                    prophage_end):
                prophage_contain_protospacer_flag = 'Yes'
                break
            elif int(protospacer_start_loc) >= int(prophage_start) and int(protospacer_end_loc) >= int(
                    prophage_end):
                prophage_contain_protospacer_flag = 'Yes'
                break
            elif int(protospacer_start_loc) <= int(prophage_start) and int(protospacer_end_loc) <= int(
                    prophage_end) and int(protospacer_end_loc) >= int(prophage_start):
                prophage_contain_protospacer_flag = 'Yes'
                break
    return prophage_contain_protospacer_flag

def get_prophage_contain_protospacer(prophage_file,prophage_contain_protospacer_file,st_info):
    if os.path.getsize(prophage_file) > 0:
        # step1 get protospacer location dict order by start location
        protospacer_location_dict = OrderedDict()
        get_protospacer_location_dict(st_info, protospacer_location_dict)
        # step2 filter the prophage don't contain the protospacer
        fout = open(prophage_contain_protospacer_file,'w')
        with open(prophage_file,'r')as fin:
            lines = fin.readlines()
            header_list = lines[0].strip('\n').split('\t')
            fout.write('\t'.join(header_list)+'\n')
            for line in lines[1:]:
                content = line.strip('\n').split('\t')
                cur_contig_id = content[header_list.index("bacteria_id")]
                cur_prophage_start = content[header_list.index("prophage_start")]
                cur_prophage_end = content[header_list.index("prophage_end")]
                if cur_contig_id not in protospacer_location_dict:
                    continue
                # check whether the prophage contains the protospacer
                prophage_contain_protospacer_flag = check_prophage_contain_protospacer(protospacer_location_dict,cur_contig_id,cur_prophage_start,cur_prophage_end)
                if prophage_contain_protospacer_flag == 'Yes':
                    fout.write('\t'.join(content)+'\n')
        fout.close()
    else:
        fout = open(prophage_contain_protospacer_file, 'w')
        header_list = ["prefix","bacteria_id","bacteria_def","genome_size","prophage_start","prophage_end",
                       "method","key_proteins","best_hit_species","CDS_number","attl_region","attr_region"]
        fout.write('\t'.join(header_list)+'\n')
        fout.close()

def acRanker(prog_dir,protein_faa_file,result_file):
    cmd_acRanker = '/usr/bin/python %s/software/acranker.py %s %s'%(protein_faa_file,result_file,prog_dir)
    os.system(cmd_acRanker)

def rpsblastpHTHproteins(prog_dir,fileName,outFileName,cdd_db,evalue):
    blast_cline = "%s/software/rpsblast -query %s -comp_based_stats 0 -max_target_seqs 1 -evalue %s -seg no -outfmt 5 -num_threads 20 -db %s -out %s" % (prog_dir,fileName,str(evalue),cdd_db,outFileName)
    os.system(blast_cline)

def merge_intervals(intervals):
    """
    :type intervals: List[Interval]
    :rtype: List[Interval]
    """
    length = 0
    intervals_sorted = sorted(intervals, key=lambda x : x[0])
    result = []
    for interval in intervals_sorted:
        if result and result[-1][1] >= interval[0]:
            result[-1][1] = max(result[-1][1], interval[1])
        else:
            result.append(interval)
    #return result
    for item in result:
        length += abs(item[1]-item[0]) + 1
    return length,result[0][0],result[-1][1]

def ParseResult(blast_file, parse_file, filter_file, identity, coverage, filter_flag='yes'):
    root = etree.parse(blast_file)
    BlastOutput_iterations = root.find("BlastOutput_iterations")
    f = open(parse_file, 'w')
    f.write(
        'protein_id\thit_id\thit_def\tq_start\tq_end\ts_start\ts_end\tevalue\tidentity\tquery_coverage\thit_coverage\n')
    if filter_flag != 'no':
        f_filter = open(filter_file, 'w')
        f_filter.write(
            'protein_id\thit_id\thit_def\tq_start\tq_end\ts_start\ts_end\tevalue\tidentity\tquery_coverage\thit_coverage\n')
    for Iteration in BlastOutput_iterations.findall("Iteration"):
        strTemp = str(Iteration.find("Iteration_query-def").text)
        query_len = str(Iteration.find("Iteration_query-len").text)
        Iteration_hits = Iteration.find("Iteration_hits")
        for Hit in Iteration_hits.findall("Hit"):
            strDef = str(Hit.find("Hit_def").text)
            Hit_len = Hit.find("Hit_len").text
            Hit_ID = Hit.find("Hit_id").text
            Hsps = Hit.find("Hit_hsps").findall("Hsp")
            Hsp_identity = 0
            Hsp_align_len = 0
            query_arr = []
            hit_arr = []
            min_Hsp_evalue = 10
            for Hsp in Hsps:
                query_from = Hsp.find("Hsp_query-from").text
                query_to = Hsp.find("Hsp_query-to").text
                Hit_from = Hsp.find("Hsp_hit-from").text
                Hit_to = Hsp.find("Hsp_hit-to").text
                query_arr.append([int(query_from), int(query_to)])
                hit_arr.append([int(Hit_from), int(Hit_to)])
                Hsp_identity = Hsp_identity + int(Hsp.find("Hsp_identity").text)
                Hsp_align_len = Hsp_align_len + int(Hsp.find("Hsp_align-len").text)
                Hsp_evalue = Hsp.find("Hsp_evalue").text
                if float(Hsp_evalue) < min_Hsp_evalue:
                    min_Hsp_evalue = float(Hsp_evalue)
            hit_length, hit_start, hit_end = merge_intervals(hit_arr)
            query_length, query_start, query_end = merge_intervals(query_arr)
            hit_coverage = str(float(hit_length) / float(Hit_len))
            query_coverage = str(float(query_length) / float(query_len))
            hit_identity = str(float(Hsp_identity) / Hsp_align_len)

            if float(query_coverage) > 1:
                query_coverage = str(1)
            if float(hit_coverage) > 1:
                hit_coverage = str(1)

            f.write(
                strTemp + '\t' + Hit_ID + '\t' + strDef + '\t' + str(query_start) + '\t' + str(query_end) + '\t' + str(
                    hit_start) + '\t' + str(hit_end) + '\t' + str(min_Hsp_evalue) + '\t' + hit_identity + '\t' + str(
                    query_coverage) + '\t' + str(hit_coverage) + '\n')
            f.flush()
            if filter_flag == 'yes':
                if (float(query_coverage) >= float(coverage)) and (float(hit_coverage) >= float(coverage)) and (
                        (float(hit_identity)) >= float(identity)):
                    f_filter.write(strTemp + '\t' + Hit_ID + '\t' + strDef + '\t' + str(query_start) + '\t' + str(
                        query_end) + '\t' + str(hit_start) + '\t' + str(hit_end) + '\t' + str(
                        min_Hsp_evalue) + '\t' + hit_identity + '\t' + str(query_coverage) + '\t' + str(
                        hit_coverage) + '\n')
                    f_filter.flush()
            if filter_flag == 'hit':
                if (float(hit_coverage) >= float(coverage)) and ((float(hit_identity)) >= float(identity)):
                    f_filter.write(strTemp + '\t' + Hit_ID + '\t' + strDef + '\t' + str(query_start) + '\t' + str(
                        query_end) + '\t' + str(hit_start) + '\t' + str(hit_end) + '\t' + str(
                        min_Hsp_evalue) + '\t' + hit_identity + '\t' + str(query_coverage) + '\t' + str(
                        hit_coverage) + '\n')
                    f_filter.flush()

    f.close()
    if filter_flag != 'no':
        f_filter.close()

def annotate_HTH(protein_file, outdir,prog_dir, identity=0.4, coverage=0.7):
    evalue = 0.01
    hth_cdd_db = "%s/database/HTH_rps/hth"%(prog_dir)
    hth_rpsblast_file = os.path.join(outdir, 'prophage_contain_protospacer_protein_rpsblast_hth.txt')
    rpsblastpHTHproteins(prog_dir, protein_file, hth_rpsblast_file, hth_cdd_db, evalue)

    protein_hth_file = os.path.join(outdir, 'prophage_contain_protospacer_protein_hth_info.txt')
    protein_hth_filter_file = os.path.join(outdir, 'prophage_contain_protospacer_protein_hth_filter_info.txt')
    ParseResult(hth_rpsblast_file, protein_hth_file, protein_hth_filter_file, str(identity), str(coverage), 'hit')

def predict_novel_acr_by_HTH(prophage_protein_file,hth_file,outdir,prot_key_index_value_dict, index_key_prot_info_value_dict,prog_dir,identity=0.4, coverage=0.7,distance_prot_num=3,acr_prot_size=400):
    # step1 annote HTH domain
    annotate_HTH(prophage_protein_file,outdir,prog_dir, identity, coverage)
    # step2 get novel Acr protein
    # get upstream or downstream protein of the Aca protein according to the Aca protein direction
    # filter the protein size
    with open(hth_file,'r')as fin:
        lines = fin.readlines()
    if len(lines) == 1:
        acr_candidate_by_hth_file = os.path.join(outdir,'prophage_contain_st_region_HTH_neighbor_%s_filter_prot_size_%s_aa.txt' % (
                                                     str(distance_prot_num), str(acr_prot_size)))
        header_list = ['neighbor_prote_id', 'distance_with_HTH'] + lines[0].strip('\n').split('\t')
        fout = open(acr_candidate_by_hth_file, 'w')
        fout.write('\t'.join(header_list)+'\n')
        fout.close()
    else:
        hth_neighbor_prot_seq_file = os.path.join(outdir, 'prophage_contain_st_region_HTH_neighbor_%s_proteins_filter_prot_size_%s_aa_seq.faa' % (
        str(distance_prot_num), str(acr_prot_size)))
        hth_neighbor_prot_info_file = os.path.join(outdir, 'prophage_contain_st_region_HTH_neighbor_%s_proteins_filter_prot_size_%s_aa_info.faa' % (
        str(distance_prot_num), str(acr_prot_size)))
        get_upstream_or_downstream_prot_file(hth_file, prot_key_index_value_dict, index_key_prot_info_value_dict,
                                             hth_neighbor_prot_seq_file, hth_neighbor_prot_info_file, distance_prot_num,
                                             acr_prot_size)
        prot_domain_dict = {}
        if os.path.getsize(hth_neighbor_prot_seq_file) > 0:
            # annotate the domain of the upstream or downstream proteins
            up_protein_domain_prefix = os.path.join(outdir, 'prophage_contain_st_region_HTH_neighbor_%s_filter_prot_size_%s_aa_rpsblast' % (
            str(distance_prot_num), str(acr_prot_size)))
            rpsblastp_protein(hth_neighbor_prot_seq_file, up_protein_domain_prefix,prog_dir)
            # get the Acr candidate dict by HTH
            prot_domain_file = '%s_domain.txt' % (up_protein_domain_prefix)
            get_prote_domain_dict(prot_domain_file, prot_domain_dict)

        acr_candidate_by_hth_file = os.path.join(outdir, 'prophage_contain_st_region_HTH_neighbor_%s_filter_prot_size_%s_aa.txt' % (
        str(distance_prot_num), str(acr_prot_size)))
        fout = open(acr_candidate_by_hth_file, 'w')
        with open(hth_neighbor_prot_info_file, 'r')as fin:
            lines = fin.readlines()
            header_list = lines[0].strip('\n').split('\t')
            fout.write('\t'.join(header_list) + '\n')
            for line in lines[1:]:
                content = line.strip('\n').split('\t')
                cur_prot_id = content[0]
                if cur_prot_id not in prot_domain_dict:
                    strwrite_list = content
                    fout.write('\t'.join(strwrite_list) + '\n')
        fout.close()

def predict_novel_acr_by_st_in_prophage(prophage_protein_file,hth_file,prot_key_index_value_dict, index_key_prot_info_value_dict,outdir,prog_dir,distance_prot_num=3,acr_prot_size=400):
    if os.path.getsize(prophage_protein_file) > 0:
        # 1. AcRanker
        acRanker_result_file_prefix = os.path.join(outdir, 'prophage_protein_AcRanker_result')
        mThread_acRanker = threading.Thread(target=acRanker, args=(prog_dir,prophage_protein_file, acRanker_result_file_prefix,))
        mThread_acRanker.start()
        # 2. HTH analysis
        predict_novel_acr_by_HTH(prophage_protein_file,hth_file,outdir, prot_key_index_value_dict, index_key_prot_info_value_dict,prog_dir,
                                 identity, coverage, distance_prot_num, acr_prot_size)
        mThread_acRanker.join()

        acr_candidate_by_hth_file = os.path.join(outdir,'prophage_contain_st_region_HTH_neighbor_%s_filter_prot_size_%s_aa.txt' % (
                                                     str(distance_prot_num), str(acr_prot_size)))
        acRanker_result_file = '%s.csv'%(acRanker_result_file_prefix)

        acr_predicted_by_prophage_contain_st = os.path.join(outdir,'acr_predict_by_prophage_contain_self-targeting.txt')
        fout = open(acr_predicted_by_prophage_contain_st,'w')
        with open(acr_candidate_by_hth_file,'r')as fin:
            lines = fin.readlines()
            header_list = lines[0].strip('\n').split('\t')
            fout.write('\t'.join(header_list) + '\n')
            if len(lines) == 1:
                if os.path.exists(acRanker_result_file):
                    with open(acRanker_result_file,'r')as fin:
                        acRanker_lines = fin.readlines()
                        acr_candidate_id = acRanker_lines[1].split(',')[0]
                        strwrite_list = [acr_candidate_id,'AcRanker'] + ['NA']*(len(header_list)-2)
                        fout.write('\t'.join(strwrite_list)+'\n')
            else:
                for line in lines[1:]:
                    fout.write(line)
        fout.close()
    else:
        acr_predicted_by_prophage_contain_st = os.path.join(outdir,'acr_predict_by_prophage_contain_self-targeting.txt')
        fout = open(acr_predicted_by_prophage_contain_st, 'w')
        header_list = ["neighbor_prote_id","distance_with_HTH","protein_id","hit_id","hit_def","q_start","q_end",
                       "s_start","s_end","evalue","identity","query_coverage","hit_coverage"]
        fout.write('\t'.join(header_list)+'\n')
        fout.close()

def get_integrate_anti_dict(acr_candidate_file_one_method,integrate_anti_dict,method):
    if not os.path.exists(acr_candidate_file_one_method):
        return 0
    with open(acr_candidate_file_one_method,'r')as fin:
        lines = fin.readlines()
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_acr_candidate = content[0]
            if cur_acr_candidate not in integrate_anti_dict:
                integrate_anti_dict[cur_acr_candidate] = {}
                integrate_anti_dict[cur_acr_candidate]['homolog_acr'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['homolog_acr_evalue'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['homolog_acr_identity'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['homolog_acr_query_coverage'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['homolog_acr_hit_coverage'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['neighbor_aca'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['homolog_aca'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['distance_aca'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['homolog_aca_evalue'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['homolog_aca_identity'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['homolog_aca_query_coverage'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['homolog_aca_hit_coverage'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['neighbor_prophage_st_hth'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['distance_prophage_st_hth'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['prophage_st_hth_evalue'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['prophage_st_hth_identity'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['prophage_st_hth_query_coverage'] = 'NA'
                integrate_anti_dict[cur_acr_candidate]['prophage_st_hth_hit_coverage'] = 'NA'
            if method == 'acr_homology':
                integrate_anti_dict[cur_acr_candidate]['homolog_acr'] = content[1]
                integrate_anti_dict[cur_acr_candidate]['homolog_acr_evalue'] = content[-4]
                integrate_anti_dict[cur_acr_candidate]['homolog_acr_identity'] = content[-3]
                integrate_anti_dict[cur_acr_candidate]['homolog_acr_query_coverage'] = content[-2]
                integrate_anti_dict[cur_acr_candidate]['homolog_acr_hit_coverage'] = content[-1]
            elif method == 'aca_homology':
                integrate_anti_dict[cur_acr_candidate]['neighbor_aca'] = content[2]
                integrate_anti_dict[cur_acr_candidate]['homolog_aca'] = content[3]
                integrate_anti_dict[cur_acr_candidate]['distance_aca'] = content[1]
                integrate_anti_dict[cur_acr_candidate]['homolog_aca_evalue'] = content[-4]
                integrate_anti_dict[cur_acr_candidate]['homolog_aca_identity'] = content[-3]
                integrate_anti_dict[cur_acr_candidate]['homolog_aca_query_coverage'] = content[-2]
                integrate_anti_dict[cur_acr_candidate]['homolog_aca_hit_coverage'] = content[-1]
            elif method == 'prophage_contain_st':
                integrate_anti_dict[cur_acr_candidate]['neighbor_prophage_st_hth'] = content[2]
                integrate_anti_dict[cur_acr_candidate]['distance_prophage_st_hth'] = content[1]
                integrate_anti_dict[cur_acr_candidate]['prophage_st_hth_evalue'] = content[-4]
                integrate_anti_dict[cur_acr_candidate]['prophage_st_hth_identity'] = content[-3]
                integrate_anti_dict[cur_acr_candidate]['prophage_st_hth_query_coverage'] = content[-2]
                integrate_anti_dict[cur_acr_candidate]['prophage_st_hth_hit_coverage'] = content[-1]

def integrate_anti_module_result(integrated_anti_module_file,acr_homology_with_known_anti_file,acr_predict_by_aca_file,acr_predict_by_prophage_contain_st_file,prot_key_index_value_dict, index_key_prot_info_value_dict):
    fout = open(integrated_anti_module_file,'w')
    header_list = ['acr_candidate_id','acr_candidate_location','acr_prot_direction','acr_candidate_protein_size','distance_with_Aca',
                   'neighbor_aca_protein','homology_with_known_aca','homology_with_known_aca_evalue','homology_with_known_aca_identity',
                   'homology_with_known_aca_query_coverage','homology_with_known_aca_hit_coverage',
                   'homology_with_known_acr','homology_with_known_acr_evalue','homology_with_known_acr_identity',
                   'homology_with_known_acr_query_coverage','homology_with_known_acr_hit_coverage',
                   'distance_with_HTH','HTH_protein','domain_evalue','domain_identity','domain_query_coverage','domain_hit_coverage']
    integrate_anti_dict = {}
    get_integrate_anti_dict(acr_homology_with_known_anti_file, integrate_anti_dict, 'acr_homology')
    get_integrate_anti_dict(acr_predict_by_aca_file, integrate_anti_dict, 'aca_homology')
    get_integrate_anti_dict(acr_predict_by_prophage_contain_st_file, integrate_anti_dict, 'prophage_contain_st')
    fout.write('\t'.join(header_list)+'\n')
    for key in integrate_anti_dict:
        cur_contig_id = '_'.join(key.split('_')[0:-3])
        cur_prot_location = '-'.join(index_key_prot_info_value_dict[cur_contig_id][prot_key_index_value_dict[cur_contig_id][key]].split('\n')[0].split('_')[-3:-1])
        cur_prot_direction = index_key_prot_info_value_dict[cur_contig_id][prot_key_index_value_dict[cur_contig_id][key]].split('\n')[0].split('_')[-1]
        cur_protein_size = str(len(index_key_prot_info_value_dict[cur_contig_id][prot_key_index_value_dict[cur_contig_id][key]].split('\n')[1]))
        strwrite_list = [key,
                         cur_prot_location,
                         cur_prot_direction,
                         cur_protein_size+' aa',
                         integrate_anti_dict[key]['distance_aca'],
                         integrate_anti_dict[key]['neighbor_aca'],
                         integrate_anti_dict[key]['homolog_aca'],
                         integrate_anti_dict[key]['homolog_aca_evalue'],
                         integrate_anti_dict[key]['homolog_aca_identity'],
                         integrate_anti_dict[key]['homolog_aca_query_coverage'],
                         integrate_anti_dict[key]['homolog_aca_hit_coverage'],
                         integrate_anti_dict[key]['homolog_acr'],
                         integrate_anti_dict[key]['homolog_acr_evalue'],
                         integrate_anti_dict[key]['homolog_acr_identity'],
                         integrate_anti_dict[key]['homolog_acr_query_coverage'],
                         integrate_anti_dict[key]['homolog_acr_hit_coverage'],
                         integrate_anti_dict[key]['distance_prophage_st_hth'],
                         integrate_anti_dict[key]['neighbor_prophage_st_hth'],
                         integrate_anti_dict[key]['prophage_st_hth_evalue'],
                         integrate_anti_dict[key]['prophage_st_hth_identity'],
                         integrate_anti_dict[key]['prophage_st_hth_query_coverage'],
                         integrate_anti_dict[key]['prophage_st_hth_hit_coverage']
                         ]
        fout.write('\t'.join(strwrite_list)+'\n')
    fout.close()

def get_prophage_region_dict(prophage_file,prophage_region_dict):
    if os.path.getsize(prophage_file) > 0:
        with open(prophage_file, 'r')as fin:
            lines = fin.readlines()
            header_list = lines[0].strip('\n').split('\t')
            for line in lines[1:]:
                content = line.strip('\n').split('\t')
                cur_contig_id = content[header_list.index("bacteria_id")]
                cur_prophage_start = content[header_list.index("prophage_start")]
                cur_prophage_end = content[header_list.index("prophage_end")]
                if cur_contig_id not in prophage_region_dict:
                    prophage_region_dict[cur_contig_id] = []
                prophage_region_dict[cur_contig_id].append('%s-%s'%(cur_prophage_start,cur_prophage_end))

def check_whether_in_prophage(candidate_info,prophage_region_dict):
    cur_contig = '_'.join(candidate_info.split('_')[0:-3]).split('.')[0]
    cur_candidate_start = candidate_info.split('_')[-3]
    cur_candodate_end = candidate_info.split('_')[-2]
    candidate_in_prophage_flag = 'No'
    candidate_prophage_region = 'No'
    if cur_contig in prophage_region_dict:
        cur_contig_prophage_region_list = prophage_region_dict[cur_contig]
        for cur_contig_prophage_region_item in cur_contig_prophage_region_list:
            cur_prophage_start = cur_contig_prophage_region_item.split('-')[0]
            cur_prophage_end = cur_contig_prophage_region_item.split('-')[1]
            if int(cur_candidate_start) >= int(cur_prophage_start) and int(cur_candodate_end) <= int(cur_prophage_end):
                candidate_in_prophage_flag = 'Yes'
                candidate_prophage_region = cur_prophage_start+'-'+cur_prophage_end
                break
            elif int(cur_candidate_start) <= int(cur_prophage_start) and int(cur_candodate_end) <= int(cur_prophage_end) and int(cur_candodate_end) > int(cur_prophage_start):
                candidate_in_prophage_flag = 'Yes'
                candidate_prophage_region = cur_prophage_start+'-'+cur_prophage_end
                break
            elif int(cur_candidate_start) >= int(cur_prophage_start) and int(cur_candidate_start) <= int(cur_prophage_end) and int(cur_candodate_end) > int(cur_prophage_end):
                candidate_in_prophage_flag = 'Yes'
                candidate_prophage_region = cur_prophage_start+'-'+cur_prophage_end
                break
    return candidate_prophage_region

def check_candodatae_protein_in_prophage(integrated_anti_module_file,integrated_anti_module_add_prophage_info_file,prophage_file):
    prophage_region_dict = {}
    get_prophage_region_dict(prophage_file, prophage_region_dict)
    fout = open(integrated_anti_module_add_prophage_info_file,'w')
    with open(integrated_anti_module_file,'r')as fin:
        lines = fin.readlines()
        header_list = lines[0].strip('\n').split('\t')[0:4]+['Acr_candidate_in_prophage']+lines[0].strip('\n').split('\t')[4:]
        fout.write('\t'.join(header_list)+'\n')
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_candidate_info = content[0]
            whether_in_prophage_flag = check_whether_in_prophage(cur_candidate_info,prophage_region_dict)
            strwrite_list = content[0:4]+[whether_in_prophage_flag]+content[4:]
            fout.write('\t'.join(strwrite_list)+'\n')
    fout.close()

def get_upstream_and_downstream_prot_file(integrated_anti_module_add_prophage_info_file,prot_key_index_value_dict, index_key_prot_info_value_dict,candidate_neighbor_prot_seq_file,candidate_neighbor_prot_info_file,neighbor_num=5):
    fout = open(candidate_neighbor_prot_seq_file,'w')
    finfo = open(candidate_neighbor_prot_info_file,'w')

    with open(integrated_anti_module_add_prophage_info_file,'r')as fin:
        lines = fin.readlines()
        header_list = ['neighbor_prote_id','distance_with_candidate']+lines[0].strip('\n').split('\t')
        finfo.write('\t'.join(header_list)+'\n')
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_aca_protein_info = content[0]
            cur_contig_id = '_'.join(cur_aca_protein_info.split('_')[0:-3])
            cur_aca_index = prot_key_index_value_dict[cur_contig_id][cur_aca_protein_info]

            for i in range(max(int(cur_aca_index)-int(neighbor_num),1),int(cur_aca_index)):
                seq_len = len(index_key_prot_info_value_dict[cur_contig_id][str(i)].split('\n')[1])
                title_info = index_key_prot_info_value_dict[cur_contig_id][str(i)].split('\n')[0].replace('>','')
                fout.write(index_key_prot_info_value_dict[cur_contig_id][str(i)])
                distance_with_aca = abs(int(i)-int(cur_aca_index))-1
                neighbor_prot_info_list = [title_info,str(distance_with_aca)]+content
                finfo.write('\t'.join(neighbor_prot_info_list)+'\n')
            for i in range(int(cur_aca_index)+1,min(int(cur_aca_index)+int(neighbor_num)+1,len(prot_key_index_value_dict[cur_contig_id]))):
                seq_len = len(index_key_prot_info_value_dict[cur_contig_id][str(i)].split('\n')[1])
                title_info = index_key_prot_info_value_dict[cur_contig_id][str(i)].split('\n')[0].replace('>', '')
                fout.write(index_key_prot_info_value_dict[cur_contig_id][str(i)])
                distance_with_aca = abs(int(i) - int(cur_aca_index))-1
                neighbor_prot_info_list = [title_info,str(distance_with_aca)] + content
                finfo.write('\t'.join(neighbor_prot_info_list) + '\n')
    fout.close()
    finfo.close()

def annotate_HTH_acr_candodate_neighbor(protein_file, outdir,prog_dir, identity=0.4, coverage=0.7,neighbor_num=5):
    evalue = 0.01
    hth_cdd_db = "%s/database/HTH_rps/hth"%(prog_dir)
    hth_rpsblast_file = os.path.join(outdir, 'acr_candodate_neighbor_%s_protein_hth_info.xml'%(neighbor_num))
    rpsblastpHTHproteins(prog_dir, protein_file, hth_rpsblast_file, hth_cdd_db, evalue)

    protein_hth_file = os.path.join(outdir, 'acr_candodate_neighbor_%s_protein_hth_info.txt'%(neighbor_num))
    protein_hth_filter_file = os.path.join(outdir, 'acr_candodate_neighbor_%s_protein_hth_filter_info.txt'%(neighbor_num))
    ParseResult(hth_rpsblast_file, protein_hth_file, protein_hth_filter_file, str(identity), str(coverage), 'hit')

def anno_candidate_neighbor_hth(integrated_anti_module_add_prophage_info_file,prot_key_index_value_dict,index_key_prot_info_value_dict,outdir,prog_dir,neighbor_num=5,identity=0.4, coverage=0.7):
    # get upstream and downstream proteins of the acr candidate protein
    candidate_neighbor_prot_seq_file = os.path.join(outdir, 'acr_candidate_neighbor_%s_proteins_seq.faa' % (str(neighbor_num)))
    candidate_neighbor_prot_info_file = os.path.join(outdir, 'acr_candidate_neighbor_%s_proteins_info.txt' % (str(neighbor_num)))
    get_upstream_and_downstream_prot_file(integrated_anti_module_add_prophage_info_file, prot_key_index_value_dict,
                                          index_key_prot_info_value_dict, candidate_neighbor_prot_seq_file,
                                          candidate_neighbor_prot_info_file, neighbor_num)
    # annotate the hth domain
    if os.path.exists(candidate_neighbor_prot_seq_file):
        if os.path.getsize(candidate_neighbor_prot_seq_file)>0:
            annotate_HTH_acr_candodate_neighbor(candidate_neighbor_prot_seq_file, outdir, prog_dir,identity, coverage,neighbor_num)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_protein', help='protein_file of genome.\n')
    parser.add_argument('--input_st_info', help='file containing self-targeting infomation.\n')
    parser.add_argument('--input_prophage', help='the out dir.\n')
    parser.add_argument('--output', help='the out dir.\n')
    parser.add_argument('--identity', help='the out dir.\n')
    parser.add_argument('--coverage', help='the out dir.\n')
    parser.add_argument('--up_num', help='the out dir.\n')
    parser.add_argument('--neighbor_num', help='the out dir.\n')
    parser.add_argument('--protein_size', help='the out dir.\n')
    parser.add_argument('--type', help='the out dir.\n') 
    parser.add_argument('--flag', help='known,novel\n')
    parser.add_argument('--prog_dir', help='directory of program.\n')

    args = parser.parse_args()
    if args.input_protein:
        protein_file = args.input_protein
    if args.prog_dir:
        prog_dir = args.prog_dir
    if args.input_st_info:
        st_info = args.input_st_info
    if args.input_prophage:
        prophage_file = args.input_prophage
    if args.output:
        outdir = args.output
    if args.identity:
        identity = args.identity
    else:
        identity = 0.4
    if args.coverage:
        coverage = args.coverage
    else:
        coverage = 0.7
    if args.up_num:
        up_num = args.up_num
    else:
        up_num = 3
    if args.neighbor_num:
        neighbor_num = args.neighbor_num
    else:
        neighbor_num = 5
    if args.protein_size:
        protein_size = args.protein_size
    else:
        protein_size = 400
    if args.type:
        data_format = args.type
    else:
        data_format = 'fasta'

    if args.flag:
        flag = args.flag
    else:
        flag = 'known'

    prot_key_index_value_dict = {}
    index_key_prot_info_value_dict = {}
    get_prote_dict(protein_file, prot_key_index_value_dict, index_key_prot_info_value_dict)

    #step1 annotate known Acr and Aca proteins in the whole bacteria genome
    predict_known_anti_hth(protein_file,outdir,identity,coverage,prog_dir)
    known_anti_file = os.path.join(outdir,'predict_known_anti_info.txt')
    known_aca_file = os.path.join(outdir,'predict_known_aca.txt')

    if flag=='novel':
    
        # step2 predict Acr candidate proteins by Aca proteins
        predict_novel_acr_by_aca(known_aca_file, prot_key_index_value_dict, index_key_prot_info_value_dict, outdir,prog_dir, up_num,protein_size)
        acr_candidate_by_aca_file = os.path.join(outdir, 'aca_neighbor_%s_filter_prot_sieze_%s_aa_rpsblast' % (str(up_num), str(protein_size)))

        # step3 get prophage containing self-targeting
        prophage_contain_protospacer_file = os.path.join(outdir,'prophage_containing_protospacer.txt')
        get_prophage_contain_protospacer(prophage_file, prophage_contain_protospacer_file, st_info)

        #step4:extract protein in prophage
        prophage_protein_file = os.path.join(outdir,'prophage_contain_protospacer_protein.faa')
        prophage_protein_info_file = os.path.join(outdir,'prophage_contain_protospacer_protein_info.txt')
        prophage_find_hth(prophage_contain_protospacer_file, protein_file, data_format, prophage_protein_file, prophage_protein_info_file)

        # step5 predice novel acr in prophage containg self-targeting
        hth_file = os.path.join(outdir,'prophage_contain_protospacer_protein_hth_filter_info.txt')
        predict_novel_acr_by_st_in_prophage(prophage_protein_file, hth_file, prot_key_index_value_dict,
                                            index_key_prot_info_value_dict, outdir,prog_dir, up_num, protein_size)

    # step6 integrated anti
    integrated_anti_module_file = os.path.join(outdir, 'acr_candidate_info_integrated_three_methods.txt')
    acr_homology_with_known_anti_file = known_anti_file
    acr_predict_by_aca_file = os.path.join(outdir, 'aca_neighbor_%s_filter_prot_size_%s_aa.txt' % (str(up_num),str(protein_size)))
    acr_predict_by_prophage_contain_st_file = os.path.join(outdir,'acr_predict_by_prophage_contain_self-targeting.txt')
    integrate_anti_module_result(integrated_anti_module_file, acr_homology_with_known_anti_file,
                                 acr_predict_by_aca_file, acr_predict_by_prophage_contain_st_file,
                                 prot_key_index_value_dict, index_key_prot_info_value_dict)

    # step7 whether the acr candidate in the prophage region
    integrated_anti_module_add_prophage_info_file = os.path.join(outdir, 'acr_candidate_info_integrated_three_methods_check_acr_in_prophage.txt')
    check_candodatae_protein_in_prophage(integrated_anti_module_file, integrated_anti_module_add_prophage_info_file,
                                         prophage_file)

    # step8 upstream and downstrem protein of acr candidate
    anno_candidate_neighbor_hth(integrated_anti_module_add_prophage_info_file, prot_key_index_value_dict,
                                index_key_prot_info_value_dict, outdir,prog_dir, neighbor_num, identity, coverage)