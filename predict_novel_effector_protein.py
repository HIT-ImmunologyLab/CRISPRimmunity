import os
import numpy as np
import sys
from lxml import etree
import re,argparse
from Bio import Phylo
import json

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


def get_filter_neighbor_protein(crispr_cas_result_file,save_dir,att_spacer_num,att_min_repeat_length,att_max_repeat_length,att_protein_size,att_neigh_distance):
	with open(crispr_cas_result_file) as f:
		contents = f.readlines()
	header = contents[0].strip().split('\t')
	
	save_protein_table_file = os.path.join(save_dir,'filter_size_distance_domain_crispr_protein_info.txt')
	save_known_repeat_file = os.path.join(save_dir,'known_repeat_repeat.fa')
	
	f_save = open(save_protein_table_file,'w')
	f_save_repeat = open(save_known_repeat_file,'w')
	
	save_header = 'protein_id\tprotein_size_aa\tdistance_crispr_number\t'
	save_header = save_header+'\t'.join(header[0:(header.index('unknown_protein_around_crispr')+1)]).strip()
	f_save.write(save_header+'\n')
	f_save.flush()
	for line in contents[1:]:
		line = line.strip().split('\t')
		contig_id = line[1]
		genome_id = line[0]
		crispr_array_locus_merge = line[header.index('crispr_array_locus_merge')]
		crispr_array_location_merge = line[header.index('crispr_array_location_merge')]
		crispr_array_start = crispr_array_location_merge.split('-')[0]
		crispr_array_end = crispr_array_location_merge.split('-')[1]
		unknown_protein_around_crispr = line[header.index('unknown_protein_around_crispr')]
		cas_prot = line[header.index('prot_within_array_20000')].split(',')
		crispr_type = line[header.index('correct_crispr_type')]
		
		spacer_max_num = line[header.index('spacer_num')]
		
		repeat_length = line[header.index('repeat_length')]
		repeat_sequence = line[header.index('consensus_repeat')].split(',')
		if crispr_type in ['Unclear','Orphan']:
			if (np.max(list(map(int,repeat_length.split(','))))>=int(att_min_repeat_length)):
				if (np.max(list(map(int,repeat_length.split(','))))<=int(att_max_repeat_length)):		
					if float(spacer_max_num)>=2:				
						for protein in unknown_protein_around_crispr.split(','):
							if protein!='NA':
								protein_size = protein.split('|')[1]
								protein_location = protein.split('|')[2]
								protein_id = '|'.join(protein.split('|')[3:])
								if int(protein_size.strip('aa'))>=int(att_protein_size):
									if int(protein_location.split('_')[1])<int(att_neigh_distance):
										f_save.write(protein_id+'\t'+protein_size.strip('aa')+'\t'+protein_location+'\t'+'\t'.join(line[0:(header.index('unknown_protein_around_crispr')+1)]).strip()+'\n')
										f_save.flush()
		else:
			for index,repeat in enumerate(list(set(repeat_sequence))):
				f_save_repeat.write('>'+genome_id+'|'+contig_id+'|'+crispr_type.replace(' ','')+'|'+crispr_array_locus_merge+'|'+crispr_array_location_merge+'|'+str(index)+'|'+str(len(repeat))+'\n')
				f_save_repeat.write(repeat+'\n')
				f_save_repeat.flush()
	f_save.close()
	f_save_repeat.close()

def run_repeat_screen_module(prog_dir,intput_flie,result_dir,prefix,previous_repeat_file,this_repeat_file,repeat_identity,repeat_coverage):
	with open(intput_flie) as f:
		contents = f.readlines()
	if len(contents)<=1:
		print('No novel candidate effector was detected!')
		sys.exit(-1)
	script = "%s/repeat_screen_module.py"%(prog_dir)
	cmd_python = "python %s --input %s --output %s --prefix %s --repeat_db %s --known_repeat_path %s --id %s --cov %s --prog_dir %s"%(script,intput_flie,result_dir,prefix,previous_repeat_file,this_repeat_file,str(repeat_identity),str(repeat_coverage),prog_dir)
	print(cmd_python)
	os.system(cmd_python)

def diamond_blastp(prog_dir,file,outfile,database,format,evalue,identity=0.3,coverage=0.6):
	num_threads = 20
	script = prog_dir + "/software/diamond blastp -d "+database+" -q "+file+" -f "+str(format)+" -e "+str(evalue)+" -o "+outfile+" -p "+str(num_threads)+" --id "+str(identity)+" -k 1000000"
	os.system(script)

def ParseResult(blast_file,parse_file,filter_file,identity,coverage,filter_flag='yes'):
	root = etree.parse(blast_file)
	BlastOutput_iterations = root.find("BlastOutput_iterations")
	f = open(parse_file,'w')
	f.write('protein_id\thit_id\thit_def\tq_start\tq_end\ts_start\ts_end\tevalue\tidentity\tquery_coverage\thit_coverage\n')
	if filter_flag!='no':
		f_filter = open(filter_file,'w')
		f_filter.write('protein_id\thit_id\thit_def\tq_start\tq_end\ts_start\ts_end\tevalue\tidentity\tquery_coverage\thit_coverage\n')
	for Iteration in BlastOutput_iterations.findall("Iteration"):
		strTemp = str(Iteration.find("Iteration_query-def").text)
		query_len = str(Iteration.find("Iteration_query-len").text)
		Iteration_hits = Iteration.find("Iteration_hits")
		for Hit in Iteration_hits.findall("Hit"):
			strDef = str(Hit.find("Hit_def").text)
			Hit_len =  Hit.find("Hit_len").text
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
				query_arr.append([int(query_from),int(query_to)])
				hit_arr.append([int(Hit_from),int(Hit_to)])
				Hsp_identity = Hsp_identity + int(Hsp.find("Hsp_identity").text)
				Hsp_align_len = Hsp_align_len + int(Hsp.find("Hsp_align-len").text)
				Hsp_evalue = Hsp.find("Hsp_evalue").text
				if float(Hsp_evalue)<min_Hsp_evalue:
					min_Hsp_evalue = float(Hsp_evalue)
			hit_length,hit_start,hit_end = merge_intervals(hit_arr)
			query_length,query_start,query_end = merge_intervals(query_arr)
			hit_coverage = str(float(hit_length)/ float(Hit_len))
			query_coverage = str(float(query_length)/ float(query_len))			
			hit_identity = str(float(Hsp_identity) / Hsp_align_len)
			
			if float(query_coverage)>1:
				query_coverage = str(1)
			if float(hit_coverage)>1:
				hit_coverage = str(1)

			f.write(strTemp+'\t'+Hit_ID+'\t'+strDef+'\t'+str(query_start)+'\t'+str(query_end)+'\t'+str(hit_start)+'\t'+str(hit_end)+'\t'+str(min_Hsp_evalue)+'\t'+hit_identity+'\t'+str(query_coverage)+'\t'+str(hit_coverage)+'\n')
			f.flush()
			if filter_flag=='yes':
				if (float(query_coverage)>=float(coverage)) and (float(hit_coverage)>=float(coverage)) and ((float(hit_identity))>=float(identity)):
					f_filter.write(strTemp+'\t'+Hit_ID+'\t'+strDef+'\t'+str(query_start)+'\t'+str(query_end)+'\t'+str(hit_start)+'\t'+str(hit_end)+'\t'+str(min_Hsp_evalue)+'\t'+hit_identity+'\t'+str(query_coverage)+'\t'+str(hit_coverage)+'\n')
					f_filter.flush()
			if filter_flag=='hit':
				if (float(hit_coverage)>=float(coverage)) and ((float(hit_identity))>=float(identity)):
					f_filter.write(strTemp+'\t'+Hit_ID+'\t'+strDef+'\t'+str(query_start)+'\t'+str(query_end)+'\t'+str(hit_start)+'\t'+str(hit_end)+'\t'+str(min_Hsp_evalue)+'\t'+hit_identity+'\t'+str(query_coverage)+'\t'+str(hit_coverage)+'\n')
					f_filter.flush()

	f.close()
	if filter_flag!='no':
		f_filter.close()

def get_protein_db_crispr(protein_cluster_file,protein_table_file):
	protein_cluster_dict = {}
	with open(protein_cluster_file) as f:
		contents = f.read().split('>Cluster')
	for cluster in contents[1:]:
		cluster_lines = cluster.strip().split('\n')[1:]
		group_ids = []
		for cluster_line in cluster_lines:
			group_id = cluster_line.split(', ')[1].split('...')[0].strip('>')
			if '*' in cluster_line:
				represent_id = group_id
			group_ids.append(group_id)
		if represent_id not in protein_cluster_dict.keys():
			protein_cluster_dict.update({represent_id:group_ids})
	
	protein_crispr_table_dict = {}
	with open(protein_table_file) as f:
		contents = f.readlines()
	protein_crispr_table_header = contents[0].strip().split('\t')
	for line in contents[1:]:
		line = line.strip().split('\t')
		protein_id = line[0]
		if protein_id not in protein_crispr_table_dict.keys():
			protein_crispr_table_dict.update({protein_id:[]})
		if line not in protein_crispr_table_dict[protein_id]:
			protein_crispr_table_dict[protein_id].append(line)

	return protein_cluster_dict,protein_crispr_table_dict,protein_crispr_table_header

def get_protein_db_sequence_dict(prog_dir):
	protein_sequence_file = "%s/data/protein_genome_database/complete_incomplete_filter_crispr_protein_merge.faa"%(prog_dir)
	with open(protein_sequence_file) as f:
		contents = f.read().strip()
	protein_sequence_dict = {}
	for protein in contents.split('\n>'):
		protein_id = protein.split('\n')[0].strip()
		protein_sequence_dict.update({protein_id:protein})
	return protein_sequence_dict

def find_homo_protein_crispr(prog_dir,filter_repeat_table_file,run_dir,outdir,homo_identity,homo_coverage):
	with open(filter_repeat_table_file) as f:
		filter_table_contents = f.readlines()
	if len(filter_table_contents)<=1:
		print('No novel candidate effector was detected!')
		sys.exit(-1)

	#step1:get filter_table protein sequences
	save_pro_sequence_file = os.path.join(outdir,"filter_known_repeat_protein_sequence.faa")
	f_save_protein_first = open(save_pro_sequence_file,'w')
	header = filter_table_contents[0].strip().split('\t')
	ori_protein_header = header
	filter_protein_dict = {}
	protein_sequence_dict = {}
	for line in filter_table_contents[1:]:
		line = line.strip().split('\t')
		protein_id = line[0]
		contig_id = line[header.index('genome_id')]
		crispr_array_locus_merge = line[header.index('crispr_array_locus_merge')]
		crispr_array_location_merge = line[header.index('crispr_array_location_merge')]
		c_crispr_neighbor_pro_file = os.path.join(run_dir,'up_down',contig_id+'_'+crispr_array_locus_merge,'protein.faa')
		
		if protein_id not in filter_protein_dict.keys():
			filter_protein_dict.update({protein_id:[]})
		filter_protein_dict[protein_id].append(line)
		
		with open(c_crispr_neighbor_pro_file) as f:
			neighbor_po_sequences = f.read().strip()
		for pro in neighbor_po_sequences.split('>')[1:]:
			pro_id = pro.split('\n')[0]
			if pro_id==protein_id:
				f_save_protein_first.write('>'+pro.strip()+'\n')
				f_save_protein_first.flush()
				protein_sequence_dict.update({pro_id:pro.strip()})
	f_save_protein_first.close()	

	#step2:find homo protein
	out_blastp_file = os.path.join(outdir,'filter_repeat_protein_blastp.txt')
	out_blastp_parse_file = os.path.join(outdir,'filter_repeat_protein_blastp_info.txt')
	out_blastp_filter_file = os.path.join(outdir,'filter_repeat_protein_blastp_info_filter.txt')
	protein_database = "%s/database/crispr_neigh_protein/crispr_neigh_protein_db"%(prog_dir)
	diamond_blastp(prog_dir,save_pro_sequence_file,out_blastp_file,protein_database,5,0.1,homo_identity,homo_coverage)
	ParseResult(out_blastp_file,out_blastp_parse_file,out_blastp_filter_file,homo_identity,homo_coverage)

	#step3:filter homo protein genome crispr
	protein_cdhit_cluster_file = "%s/data/protein_genome_database/complete_incomplete_filter_crispr_protein_merge_cdhit.clstr"%(prog_dir)
	protein_crispr_table_file = "%s/data/protein_genome_database/complete_incomplete_filter_crispr_protein_merge_cdhit_add_annotation_cas_domain_table.txt"%(prog_dir)
	protein_cluster_dict,protein_crispr_table_dict,protein_crispr_table_header = get_protein_db_crispr(protein_cdhit_cluster_file,protein_crispr_table_file)
	
	new_protein_crispr_table_header = ['homo_'+item for item in protein_crispr_table_header]

	save_homo_protein_file = os.path.join(outdir,'filter_repeat_homo_protein_crispr_info.txt')
	save_homo_protein_filter_file = os.path.join(outdir,'filter_repeat_homo_protein_crispr_filter_result.txt')
	save_homo_protein_filter_sequence_file = os.path.join(outdir,'filter_repeat_homo_protein_crispr_filter_sequence.faa')
	
	save_header = '\t'.join(ori_protein_header)+'\thomo_protein_identity\thomo_protein_query_coverage\thomo_protein_hit_coverage\t'+'\t'.join(new_protein_crispr_table_header)
	f_save = open(save_homo_protein_file,'w')
	f_save.write('\t'.join(ori_protein_header)+'\thomo_protein_identity\thomo_protein_query_coverage\thomo_protein_hit_coverage\t'+'\t'.join(new_protein_crispr_table_header)+'\n')
	f_save_filter = open(save_homo_protein_filter_file,'w')
	f_save_filter.write('\t'.join(ori_protein_header)+'\thomo_protein_identity\thomo_protein_query_coverage\thomo_protein_hit_coverage\t'+'\t'.join(new_protein_crispr_table_header)+'\n')
	f_save_filter_sequence = open(save_homo_protein_filter_sequence_file,'w')
	
	save_protein_ids = []
	candidate_protein_dict = {}
	if os.path.exists(out_blastp_filter_file):
		with open(out_blastp_filter_file) as f:
			contents = f.readlines()
		if len(contents)>1:
			protein_db_sequence_dict = get_protein_db_sequence_dict(prog_dir)
			for line in contents[1:]:
				line = line.strip().split('\t')
				candidate_protein_id = line[0]
				hit_protein_id = line[1]
				if candidate_protein_id not in candidate_protein_dict.keys():
					candidate_protein_dict.update({candidate_protein_id:[]})
				candidate_protein_dict[candidate_protein_id].append(line)
	
	for candidate_protein_id in candidate_protein_dict.keys():			
		flag = 1
		candidate_protein_info_list = filter_protein_dict[candidate_protein_id]	
		for hit_info in candidate_protein_dict[candidate_protein_id]:
			hit_protein_id = hit_info[1]
			identity = hit_info[-3]
			query_coverage = hit_info[-2]
			hit_coverage = hit_info[-1]				
			for hit_protein_group_id in protein_cluster_dict[hit_protein_id]:
				hit_protein_crispr_info = protein_crispr_table_dict[hit_protein_group_id]
				
				for candidate_protein_info in candidate_protein_info_list:
					for info in hit_protein_crispr_info:
						f_save.write('\t'.join(candidate_protein_info)+'\t'+identity+'\t'+query_coverage+'\t'+hit_coverage+'\t'+'\t'.join(info)+'\n')
						f_save.flush()	
						if 'Type' in info[protein_crispr_table_header.index('correct_crispr_type')]:
							flag = 0
							break
						if 'cas' in info[protein_crispr_table_header.index('cas_name')].lower():
							flag = 0
							break
		if flag==1:
			for hit_info in candidate_protein_dict[candidate_protein_id]:
				hit_protein_id = hit_info[1]
				identity = hit_info[-3]
				query_coverage = hit_info[-2]
				hit_coverage = hit_info[-1]
				for hit_protein_group_id in protein_cluster_dict[hit_protein_id]:
					hit_protein_crispr_info = protein_crispr_table_dict[hit_protein_group_id]
					for candidate_protein_info in candidate_protein_info_list:
						for info in hit_protein_crispr_info:
							if candidate_protein_info[ori_protein_header.index('genome_id')]==info[protein_crispr_table_header.index('genome_id')]:
								if candidate_protein_info[ori_protein_header.index('crispr_array_location_merge')].split('-')[0]==info[protein_crispr_table_header.index('crispr_array_location_merge')].split('-')[0]:
									continue
							if info[protein_crispr_table_header.index('array_in_prot')]=='yes':
								continue
							f_save_filter.write('\t'.join(candidate_protein_info)+'\t'+identity+'\t'+query_coverage+'\t'+hit_coverage+'\t'+'\t'.join(info)+'\n')
							f_save_filter.flush()
							if candidate_protein_id not in save_protein_ids:
								save_protein_ids.append(candidate_protein_id)
								c_protein_sequence = protein_sequence_dict[candidate_protein_id]
								f_save_filter_sequence.write('>'+c_protein_sequence.strip().strip('>')+'\n')
								f_save_filter_sequence.flush()
							if hit_protein_group_id not in save_protein_ids:
								save_protein_ids.append(hit_protein_group_id)
								c_protein_sequence = protein_db_sequence_dict[hit_protein_group_id]
								f_save_filter_sequence.write('>'+c_protein_sequence.strip().strip('>')+'\n')
								f_save_filter_sequence.flush()										
	f_save.close()
	f_save_filter.close()
	f_save_filter_sequence.close()

	#plot crispr array protein_dict
	save_homo_protein_dict_file = os.path.join(outdir,'filter_repeat_homo_protein_crispr_filter_result_cytoscape.txt')
	get_homo_protein_dict(prog_dir,save_homo_protein_filter_file,save_homo_protein_dict_file)

def get_homo_protein_dict(prot_dir,save_homo_protein_filter_file,save_homo_protein_dict_file):
	homo_dict_file = "%s/data/protein_genome_database/save_effector_protein_db_info_dict.npy"%(prot_dir)
	homo_dict = np.load(homo_dict_file,allow_pickle=True).item()
	with open(save_homo_protein_filter_file) as f:
		contents = f.readlines()
	header = contents[0].strip().split('\t')
	save_homo_dict = {}
	for line in contents[1:]:
		line = line.strip().split('\t')
		homo_protein_id = line[header.index('homo_protein_id')]
		homo_crispr_id = line[header.index('homo_genome_id')]+'_'+line[header.index('homo_crispr_array_locus_merge')]
		homo_protein_crispr_dict = homo_dict[homo_protein_id][homo_crispr_id]
		if homo_protein_id not in save_homo_dict.keys():
			save_homo_dict.update({homo_protein_id:{}})
		if homo_crispr_id not in save_homo_dict[homo_protein_id].keys():
			save_homo_dict[homo_protein_id].update({homo_crispr_id:homo_protein_crispr_dict})

	js = json.dumps(save_homo_dict) 
	file = open(save_homo_protein_dict_file, 'w') 
	file.write(js) 
	file.close()

def mkdir(dirname):
	command = "mkdir -p %s"%(dirname)
	os.system(command)

def get_merge_protein_info(tree_file,merge_table_file):
	tree = Phylo.read(tree_file, "newick")
	all_nodes = tree.get_terminals()
	f_table = open(merge_table_file,'w')
	f_table.write('nodes\tcrispr_type\teffector\tfamily\n')
	for node_name in all_nodes:
		node_name = node_name.name
		if 'cas' in node_name.lower() and (len(node_name.split('|'))>=2):
			cas_name = node_name.split('|')[0]
			cas_type = node_name.split('|')[1]
			if cas_type=='V':
				f_table.write(node_name+'\t'+cas_type+'(Cas14d-Cas14u)'+'\t'+cas_name+'\n')
			else:
				f_table.write(node_name+'\t'+cas_type+'('+cas_name+')'+'\t'+cas_name+'\n')
		else:
			f_table.write(node_name+'\tNovel class 2 effector protein\t'+node_name.strip('|')+'\n')
		f_table.flush()
	f_table.close()

def multialign_analysis(prog_dir,prot_seq_file,save_prefix):
	prot_seq_msa_file = save_prefix+'_msa.txt'
	prot_seq_mview_txt_file = save_prefix+'_mview_txt.txt'
	prot_seq_mview_html_file = save_prefix+'_mview.html'
	cmd_mafft = '%s/software/mafft --reorder --adjustdirection --thread 20 "%s" > "%s"' % (prog_dir,prot_seq_file,prot_seq_msa_file)
	os.system(cmd_mafft)

def plot_tree(prog_dir,out_tree_file,merge_table_file,outdir):
	script = "%s/plot_tree.R"%(prog_dir)
	command = "Rscript %s %s %s %s"%(script,out_tree_file,merge_table_file,outdir)
	os.system(command)

def fasttree(prog_dir,align_fasta_file,out_tree):
	command = "%s/software/FastTree -n 20 %s > %s"%(prog_dir,align_fasta_file,out_tree)
	os.system(command)

def build_tree(prog_dir,candidate_protein_file,outdir):
	known_effector_file = "%s/data/II-V-VI_effector/cas_II_V_VI_effector_protein_sequence_cdhit"%(prog_dir)
	merge_file = os.path.join(outdir,'merge_novel_known_effector.faa')
	merge_table_file = os.path.join(outdir,'merge_novel_known_effector_info.txt')
	with open(candidate_protein_file) as f:
		candidate_contents = f.read()
	with open(known_effector_file) as f:
		known_effector_contents = f.read()
	f_sequence = open(merge_file,'w')
	f_sequence.write(candidate_contents.strip().replace('[','').replace(']','').replace('(','_').replace(')','').replace(':','_')+'\n'+known_effector_contents.strip())
	f_sequence.flush()
	f_sequence.close()
	
	mafft_save_prefix = os.path.join(outdir,'mafft')
	multialign_analysis(prog_dir,merge_file,mafft_save_prefix)

	align_fasta_file = mafft_save_prefix+'_msa.txt'
	out_tree_file = mafft_save_prefix+'_fasttree.txt'
	fasttree(prog_dir,align_fasta_file,out_tree_file)

	merge_table_file = os.path.join(outdir,'merge_novel_known_effector_table.txt')
	get_merge_protein_info(out_tree_file,merge_table_file)

	plot_tree(prog_dir,out_tree_file,merge_table_file,outdir)


def predict(prog_dir,crispr_cas_table_file,work_dir,att_spacer_num,att_min_repeat_length,att_max_repeat_length,att_protein_size,att_neigh_distance,repeat_identity,repeat_coverage,homo_identity,homo_coverage):
	#step1:filter crispr and protein
	save_dir = os.path.join(work_dir,'novel_effector_protein')
	mkdir(save_dir)
	get_filter_neighbor_protein(crispr_cas_table_file,save_dir,att_spacer_num,att_min_repeat_length,att_max_repeat_length,att_protein_size,att_neigh_distance)
	save_protein_table_file = os.path.join(save_dir,'filter_size_distance_domain_crispr_protein_info.txt')
	save_known_repeat_file = os.path.join(save_dir,'known_repeat_repeat.fa')

	#step2:filter repeat
	repeat_result_dir = os.path.join(save_dir,'repeat_analysis')
	mkdir(repeat_result_dir)
	previous_repeat_file = "%s/database/known_repeat/known_repeat_db"%(prog_dir)
	run_repeat_screen_module(prog_dir,save_protein_table_file,repeat_result_dir,"repeat",previous_repeat_file,save_known_repeat_file,repeat_identity,repeat_coverage)

	#step3:analyse hmology protein
	filter_repeat_table_file = os.path.join(repeat_result_dir,'repeat_filter_known_repeat.txt')
	find_homo_protein_crispr(prog_dir,filter_repeat_table_file,work_dir,save_dir,homo_identity,homo_coverage)

	#step4:plot tree
	candidate_protein_file = os.path.join(save_dir,'filter_repeat_homo_protein_crispr_filter_sequence.faa')
	build_tree(prog_dir,candidate_protein_file,save_dir)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',help="input file containing crispr cas info needed to screening.\n")
	parser.add_argument('--output',help="abspath of result directory.\n")
	parser.add_argument('--id',help='short-blastn identity.\n',default='0.9')
	parser.add_argument('--cov',help='short-blastn coverage.\n',default='0.9')
	parser.add_argument('--identity', help='diamond blastp identity.\n', default='0.9')
	parser.add_argument('--coverage',help='diamond blastp coverage.\n',default='0.9')
	parser.add_argument('--att_spacer_num',help='diamond blastp coverage.\n',default='0.9')
	parser.add_argument('--att_min_repeat_length',help='diamond blastp coverage.\n',default='0.9')
	parser.add_argument('--att_max_repeat_length',help='diamond blastp coverage.\n',default='0.9')
	parser.add_argument('--att_protein_size',help='diamond blastp coverage.\n',default='0.9')
	parser.add_argument('--att_neigh_distance',help='diamond blastp coverage.\n',default='0.9')
	parser.add_argument('--prog_dir', help='directory of program.\n')

	args = parser.parse_args()	
	if args.input:
		crispr_cas_table_file = args.input	
	if args.output:
		work_dir = args.output
	if args.prog_dir:
		prog_dir = args.prog_dir
	if args.id:
		repeat_identity = args.id
	else:
		repeat_identity = 0.9
	if args.cov:
		repeat_coverage = args.cov
	else:
		repeat_coverage = 0.9
	if args.identity:
		homo_identity = args.identity
	else:
		homo_identity = 0.3
	if args.coverage:
		homo_coverage = args.coverage
	else:
		homo_coverage = 0.6
	if args.att_spacer_num:
		att_spacer_num = args.att_spacer_num
	else:
		att_spacer_num = 2
	if args.att_min_repeat_length:
		att_min_repeat_length = args.att_min_repeat_length
	else:
		att_min_repeat_length = 18
	
	if args.att_max_repeat_length:
		att_max_repeat_length = args.att_max_repeat_length
	else:
		att_max_repeat_length = 45
	
	if args.att_protein_size:
		att_protein_size = args.att_protein_size
	else:
		att_protein_size = 500

	if args.att_neigh_distance:
		att_neigh_distance = args.att_neigh_distance
	else:
		att_neigh_distance = 5
	predict(prog_dir,crispr_cas_table_file,work_dir,att_spacer_num,att_min_repeat_length,att_max_repeat_length,att_protein_size,att_neigh_distance,repeat_identity,repeat_coverage,homo_identity,homo_coverage)
