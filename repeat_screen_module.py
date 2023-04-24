#!/usr/bin
#-*-coding:utf-8-*-

import os
import argparse
import threading
import xml.etree.ElementTree as ET
import re
import multiprocessing

def short_blastn(prog_dir,query_file,subject_file,out_file,outfmt=5,evalue='10'):
    cmd_short_blastn = '%s/software/blastn -query %s -subject %s -evalue %s -outfmt %s -out %s -word_size 7 -dust no -soft_masking FALSE -gapopen 10 -penalty -1 -gapextend 2 -task blastn-short -max_target_seqs 1' \
                       % (prog_dir,query_file, subject_file, evalue, outfmt, out_file)
    os.system(cmd_short_blastn)

def short_blastn_db(prog_dir,query_file,db_file,out_file,outfmt=5,evalue='10'):
    cmd_short_blastn = '%s/software/blastn -query %s -db %s -evalue %s -outfmt %s -num_threads 30 -out %s -word_size 7 -dust no -soft_masking FALSE -gapopen 10 -penalty -1 -gapextend 2 -task blastn-short -max_target_seqs 1' \
                       % (prog_dir,query_file, db_file, evalue, outfmt, out_file)
    os.system(cmd_short_blastn)

def parse_blastn_outfmt_5(xml_file,txt_file):
    text = open(xml_file).read()
    text = re.sub(u"[\x00-\x08\x0b-\x0c\x0e-\x1f]+", u" ", text)
    root = ET.fromstring(text)
    BlastOutput_iterations = root.find("BlastOutput_iterations")
    fout = open(txt_file,'w')
    header_list = ['query_repeat_info','hit_repeat_info','hit_crispr_type','evalue','identity','query_coverage','hit_coverage']
    fout.write('\t'.join(header_list)+'\n')
    for Iteration in BlastOutput_iterations.findall("Iteration"):
        query_info = str(Iteration.find("Iteration_query-def").text)
        Iteration_hits = Iteration.find("Iteration_hits")
        Iteration_query_len = Iteration.find('Iteration_query-len').text
        for Hit in Iteration_hits.findall("Hit"):
            hit_info = str(Hit.find("Hit_def").text)
            hit_crispr_type = hit_info.split('|')[2]
            hit_len = Hit.find("Hit_len").text
            Hit_ID = Hit.find("Hit_id").text
            Hsp = Hit.find("Hit_hsps").find("Hsp")
            hit_from = Hsp.find("Hsp_hit-from").text
            hit_to = Hsp.find("Hsp_hit-to").text
            query_from = Hsp.find("Hsp_query-from").text
            query_to = Hsp.find("Hsp_query-to").text
            Hsp_evalue = Hsp.find("Hsp_evalue").text
            hit_length = abs(int(hit_to) - int(hit_from)) + 1
            hit_coverage = str(float(hit_length) / float(hit_len))
            query_length = abs(int(query_to) - int(query_from)) + 1
            query_coverage = str(float(query_length)/float(Iteration_query_len))
            identity = str(float(Hsp.find("Hsp_identity").text) / float(Hsp.find("Hsp_align-len").text))
            strwrite_list = [query_info,hit_info,hit_crispr_type,Hsp_evalue,identity,query_coverage,hit_coverage]
            fout.write('\t'.join(strwrite_list)+'\n')
    fout.close()

def blastn(prog_dir,query_file, subject_file, xml_file,txt_file,outfmt=5, evalue='10'):
    short_blastn(prog_dir,query_file, subject_file,xml_file)
    if os.path.exists(xml_file):
        parse_blastn_outfmt_5(xml_file, txt_file)

def blastn_db(prog_dir,query_file, subject_file, xml_file,txt_file,outfmt=5, evalue='10'):
    short_blastn_db(prog_dir,query_file, subject_file,xml_file)
    if os.path.exists(xml_file):
        parse_blastn_outfmt_5(xml_file, txt_file)

def create_query_file(input_file,query_file):
    fout = open(query_file, 'w')
    with open(input_file,'r')as fin:
        lines = fin.readlines()
        header_list = lines[0].strip('\n').split('\t')
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            if "assembly_id" in header_list:
                cur_assembly_id = content[header_list.index("assembly_id")]
            else:
                cur_assembly_id = content[header_list.index("strain_id")]

            cur_genome_id = content[header_list.index("genome_id")]
            ###  correct crispr type
            cur_crispr_type = content[header_list.index("correct_crispr_type")]
            #cur_crispr_type = content[header_list.index("crispr_type_by_cas_prot")]
            cur_crispr_array = content[header_list.index("crispr_array_locus_merge")]
            cur_crispr_location = content[header_list.index("crispr_array_location_merge")]
            cur_repeat_seq_list = content[header_list.index("consensus_repeat")].split(',')
            index = 1
            for cur_repeat_seq_item in cur_repeat_seq_list:
                cur_repeat_len = len(cur_repeat_seq_item)
                strwrite = '>%d|%s|%s|%s|%s|%s|%d\n%s\n'%(index,cur_assembly_id,cur_genome_id,cur_crispr_type,cur_crispr_array,cur_crispr_location,cur_repeat_len,cur_repeat_seq_item)
                fout.write(strwrite)
                index = index + 1
    fout.close()

def mkdir(dirPath):
    cmd_mkdir = 'mkdir -p %s'%dirPath
    os.system(cmd_mkdir)

def rm_file(file_path):
    if os.path.exists(file_path):
        cmd_rm = 'rm -rf %s'%file_path
        os.system(cmd_rm)

def get_blastn_result_list(blastn_file,short_blastn_id,short_blastn_cov):
    blastn_result_list = ['NA']*6
    if not os.path.exists(blastn_file):
        return blastn_result_list
    with open(blastn_file,'r')as fin:
        lines = fin.readlines()
        if len(lines) == 1:
            return blastn_result_list
        blastn_result_ori_list = lines[1].strip('\n').split('\t')
        blastn_ori_id = blastn_result_ori_list[-3]
        blastn_ori_query_cov = blastn_result_ori_list[-2]
        blastn_ori_hit_cov = blastn_result_ori_list[-1]
        if float(blastn_ori_id) >= float(short_blastn_id) and float(blastn_ori_query_cov) >= float(short_blastn_cov) and float(blastn_ori_hit_cov) >= float(short_blastn_cov):
            blastn_result_list = blastn_result_ori_list[2:]+[blastn_result_ori_list[1]]
        return blastn_result_list

def get_this_round_known_repeat_file(input_file,this_round_known_repeat_file):
    if os.path.exists(this_round_known_repeat_file):
        rm_file(this_round_known_repeat_file)
    
    with open(input_file, 'r')as fin:
        lines = fin.readlines()
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_repeat_own_strain_file = content[-1]
            cmd_cat = 'cat %s >> %s'%(cur_repeat_own_strain_file,this_round_known_repeat_file)
            os.system(cmd_cat)

def get_repeat_matching_dict(repeat_matching_result_file,repeat_matching_dict,short_blastn_id,short_blastn_cov):
    if os.path.exists(repeat_matching_result_file):
        with open(repeat_matching_result_file,'r')as fin:
            lines = fin.readlines()
            for line in lines[1:]:
                content = line.strip('\n').split('\t')
                cur_query_info_list = content[0].split('|')
                cur_key = '%s|%s|%s|%s'%(cur_query_info_list[1],cur_query_info_list[2],cur_query_info_list[-3],cur_query_info_list[-2])
                cur_id = content[-3]
                cur_query_cov = content[-2]
                cur_hit_cov = content[-1]
                if float(cur_id) >= float(short_blastn_id) and float(cur_query_cov) >= float(short_blastn_cov) and float(cur_hit_cov)>=float(short_blastn_cov):
                    if cur_key not in repeat_matching_dict:
                        repeat_matching_dict[cur_key] = content[2:]+[content[1]]

def repeat_screening(prog_dir,input_file,this_round_known_repeat_file,output_dir,prefix,evalue,short_blastn_id,short_blastn_cov,repeat_previous_strain_file):
    # step1 create repeat file
    tmp_file_dir = '%s/tmp'% (output_dir)
    mkdir(tmp_file_dir)
    query_file = '%s/all_query_repeat_seq.fa'%(tmp_file_dir)
    create_query_file(input_file, query_file)
    muthread_create_query_file = threading.Thread(target=create_query_file, args=(input_file, query_file,))
    muthread_create_query_file.start()

    muthread_create_query_file.join()

    # step2 blastn
    xml_compare_with_own_strain_file = '%s/repeat_compare_with_own_strain.xml' % (tmp_file_dir)
    txt_compare_with_own_strain_file = '%s/repeat_compare_with_own_strain.txt' % (tmp_file_dir)
    muthread_compare_with_own_strain = threading.Thread(target=blastn, args=(
    prog_dir,query_file, this_round_known_repeat_file, xml_compare_with_own_strain_file, txt_compare_with_own_strain_file, evalue,))

    xml_compare_with_previous_strain_file = '%s/repeat_compare_with_previous_strain.xml' % (tmp_file_dir)
    txt_compare_with_previous_strain_file = '%s/repeat_compare_with_previous_strain.txt' % (tmp_file_dir)
    muthread_compare_with_previous_strain = threading.Thread(target=blastn_db, args=(
        prog_dir,query_file, repeat_previous_strain_file, xml_compare_with_previous_strain_file, txt_compare_with_previous_strain_file,
        evalue,))

    muthread_compare_with_own_strain.start()
    muthread_compare_with_previous_strain.start()

    muthread_compare_with_own_strain.join()
    muthread_compare_with_previous_strain.join()

    print("write repeat matching infomation to file\n")
    # step3 add repeat matching info

    repeat_with_own_strains_mathcing_dict = {}
    muthread_repeat_with_own_strains_mathcing = threading.Thread(target=get_repeat_matching_dict,args=
    (txt_compare_with_own_strain_file,repeat_with_own_strains_mathcing_dict,short_blastn_id,short_blastn_cov,))

    repeat_with_previous_strains_mathcing_dict = {}
    muthread_repeat_with_previous_strains_mathcing = threading.Thread(target=get_repeat_matching_dict, args=(
    txt_compare_with_previous_strain_file, repeat_with_previous_strains_mathcing_dict,short_blastn_id,short_blastn_cov,))

    muthread_repeat_with_own_strains_mathcing.start()
    muthread_repeat_with_previous_strains_mathcing.start()

    muthread_repeat_with_own_strains_mathcing.join()
    muthread_repeat_with_previous_strains_mathcing.join()


    result_file = '%s/%s_add_repeat_type_flag.txt'%(output_dir,prefix)
    fout = open(result_file,'w')
    with open(input_file,'r')as fin:
        lines = fin.readlines()
        header_list = lines[0].strip('\n').split('\t')
        repeat_type_index = header_list.index("repeat_type")
        if repeat_type_index == len(header_list)-1:
            header_write_list = header_list + \
                                ['repeat_homo_with_own_strain','repeat_homo_with_own_strain_evalue','repeat_homo_with_own_strain_id','repeat_homo_with_own_strain_query_cov','repeat_homo_with_own_strain_hit_cov','repeat_homo_with_own_detail',
                                 'repeat_homo_with_previous_strain','repeat_homo_with_previous_strain_evalue','repeat_homo_with_previous_strain_id','repeat_homo_with_previous_strain_query_cov','repeat_homo_with_previous_strain_hitcov','repeat_homo_with_previous_detail']
        else:
            header_write_list = header_list[0:header_list.index("repeat_type")+1]+\
                                ['repeat_homo_with_own_strain','repeat_homo_with_own_strain_evalue','repeat_homo_with_own_strain_id','repeat_homo_with_own_strain_query_cov','repeat_homo_with_own_strain_hit_cov','repeat_homo_with_own_detail',
                                 'repeat_homo_with_previous_strain','repeat_homo_with_previous_strain_evalue','repeat_homo_with_previous_strain_id','repeat_homo_with_previous_strain_query_cov','repeat_homo_with_previous_strain_hitcov','repeat_homo_with_previous_detail']+\
                                header_list[header_list.index("repeat_type")+1:]
        fout.write('\t'.join(header_write_list)+'\n')
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            if 'strain_id' in header_list:
                cur_assembly_id = content[header_list.index("strain_id")]
            else:
                cur_assembly_id = content[header_list.index("assembly_id")]

            cur_genome_id = content[header_list.index("genome_id")]
            cur_crispr_array = content[header_list.index("crispr_array_locus_merge")]
            cur_crispr_location = content[header_list.index("crispr_array_location_merge")]
            cur_key = '%s|%s|%s|%s'%(cur_assembly_id,cur_genome_id,cur_crispr_array,cur_crispr_location)
            cur_blastn_own_strain_result = ['NA']*6
            cur_blastn_previous_strain_result = ['NA']*6

            if cur_key in repeat_with_own_strains_mathcing_dict:
                cur_blastn_own_strain_result = repeat_with_own_strains_mathcing_dict[cur_key]

            if cur_key in repeat_with_previous_strains_mathcing_dict:
                cur_blastn_previous_strain_result = repeat_with_previous_strains_mathcing_dict[cur_key]

            if repeat_type_index == len(header_list) - 1:
                strwrite_list = content+cur_blastn_own_strain_result+cur_blastn_previous_strain_result
            else:
                strwrite_list = content[0:header_list.index(
                    "repeat_type") + 1] + cur_blastn_own_strain_result + cur_blastn_previous_strain_result + content[header_list.index("repeat_type") + 1:]
            fout.write('\t'.join(strwrite_list)+'\n')
    fout.close()

    print("filter items according to repeat matching\n")
    
    # step4 screenint the input file according to repeat matching info
    result_filter_file = '%s/%s_filter_known_repeat.txt' % (output_dir, prefix)
    fout = open(result_filter_file,'w')
    with open(result_file,'r')as fin:
        lines = fin.readlines()
        header_list = lines[0].strip('\n').split('\t')
        fout.write('\t'.join(header_list)+'\n')
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_repeat_homon_with_own_strain_flag = content[header_list.index("repeat_homo_with_own_strain")]
            cur_repeat_homon_with_db_flag = content[header_list.index("repeat_homo_with_previous_strain")]
            if cur_repeat_homon_with_own_strain_flag == 'NA' and cur_repeat_homon_with_db_flag=='NA':
                fout.write(line)
    fout.close()

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',help="input file containing crispr cas info needed to screening.\n")
    parser.add_argument('--output',help="abspath of result directory.\n")
    parser.add_argument('--prefix', help="prefix of result file.\n")
    parser.add_argument('--evalue', help='short-blastn evalue.\n', default='10')
    parser.add_argument('--id',help='short-blastn identity.\n',default='0.9')
    parser.add_argument('--cov',help='short-blastn coverage.\n',default='0.9')
    parser.add_argument('--repeat_db',help="abspath of result directory.\n")
    parser.add_argument('--known_repeat_path',help="abspath of result directory.\n")
    parser.add_argument('--prog_dir', help='directory of program.\n')
    
    args = parser.parse_args()
    if args.input:
        input_file = args.input
    if args.output:
        output_dir = args.output
    if args.prog_dir:
        prog_dir = args.prog_dir
    if args.prefix:
        prefix = args.prefix
    else:
        prefix = output_dir.split('/')[-1]
    if args.evalue:
        short_blastn_evalue = args.evalue
    else:
        short_blastn_evalue = '10'
    if args.id:
        short_blastn_id = args.id
    else:
        short_blastn_id = '0.9'
    if args.cov:
        short_blastn_cov = args.cov
    else:
        short_blastn_cov = '0.9'
    
    if args.repeat_db:
        repeat_previous_strain_file = args.repeat_db
        
    if args.known_repeat_path:
        this_round_known_repeat_file = args.known_repeat_path

    mkdir(output_dir)
    repeat_screening(prog_dir,input_file,this_round_known_repeat_file, output_dir, prefix, short_blastn_evalue, short_blastn_id, short_blastn_cov, repeat_previous_strain_file)