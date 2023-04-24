import os,re,argparse
import numpy as np
from Bio import Entrez, SeqIO
import threading
import Levenshtein
import subprocess

def spacer_blastn_phage(prog_dir,spacer_file,outfile,format):
	phage_database = "%s/database/phage_plasmid_nucl/phage_plasmid_nucl_db"%(prog_dir)
	command = '%s/software/blastn -query %s -db %s -evalue 10 -outfmt %s -out %s -word_size 7 -dust no -soft_masking FALSE -gapopen 10 -penalty -1 -gapextend 2 -ungapped -num_threads 20 -task blastn-short' % (prog_dir,spacer_file,phage_database,str(format),outfile)
	os.system(command)

def get_num(string):
	strlist = re.findall("\d+\.?\d*",str(string))
	return strlist

def parse_spacer_blastn_phage_0(file,resufile,filter_mismatch,filter_coverage):
	with open(file) as f:
		contents = f.read().strip()
	spacers_hits = contents.split('Query=')[1:]
	hit_dict = {}
	f_result = open(resufile,'w')
	f_result.write('bac_id\tbac_def\thit_phage_id\thit_phage_def\tspacer\tspacer_length\tidentity\tcoverage\tmismatch\tevalue\tbit-score\n')
	for each_spacer_hits in spacers_hits:
		if 'No hits found' in each_spacer_hits:
			continue
		spacer_name = each_spacer_hits.split('\n')[0].strip()
		bac_id = spacer_name.split('|')[-1].split('.')[0]
		try:
			bac_def = strain_inf_dict[bac_id]
		except:
			bac_def = bac_id
		spacer_length = float(spacer_name.split('|')[2])
		spacer_hits = each_spacer_hits.split('> ')[1:]
		for spacer_hit in spacer_hits:
			hit_phage = spacer_hit.split('\n')[0].strip()
			hit_phage_id = hit_phage.split()[0].split('.')[0]
			hit_phage_def = ' '.join(hit_phage.split()[1:])
			bit_score = float(get_num(spacer_hit.split('\n\n')[1].split('\n')[0])[0])
			evalue = float(spacer_hit.split('\n\n')[1].split('\n')[0].split('=')[-1].strip())
			identity = float(get_num(spacer_hit.split('\n\n')[1].split('\n')[1])[2])
			hit_length1 = int(get_num(spacer_hit.split('\n\n')[1].split('\n')[1])[0])
			hit_length = int(get_num(spacer_hit.split('\n\n')[1].split('\n')[1])[1])
			mismatch = hit_length-hit_length1
			coverage = hit_length/spacer_length
			hit_sequence_details = spacer_hit.split('\n\n')[2]
			positions = get_num(hit_sequence_details)
			hit_details = hit_sequence_details.split('\n')
			query = hit_details[0].split()[2].strip()
			subj = hit_details[-1].split()[2].strip()
			hit_info= query+'\t'+'Query'+'\t'+str(positions[0])+':'+str(positions[1])+'\n'+hit_details[1][len(hit_details[1])-len(query):]+'\n'+subj+'\t'+'Sbjct'+'\t'+str(positions[2])+':'+str(positions[3])             
			if (int(mismatch) <= int(filter_mismatch)) and (coverage>=float(filter_coverage)):
				if bac_id not in hit_dict.keys():
					hit_dict.update({bac_id:{}})
				if hit_phage_id not in hit_dict[bac_id].keys():
					hit_dict[bac_id].update({hit_phage_id:[]})
				resuline = list(map(str,[bac_id,bac_def,hit_phage_id,hit_phage_def,spacer_name,spacer_length,hit_length,identity,coverage,mismatch,evalue,bit_score]))				
				f_result.write('\t'.join(resuline)+'\n')
				f_result.flush()	
	f_result.close()

def parse_spacer_blastn_phage_6(prog_dir,file,resufile,filter_mismatch,filter_coverage):
	phage_info_file = "%s/data/phage_plasmid/ncbi_bacteriophage_PLSDB_plasmid_info.txt"%(prog_dir)
	phage_inf_dict = {}
	with open(phage_info_file) as f:
		contents = f.readlines()
	for line in contents[1:]:
		line = line.strip().split('\t')
		phage_id = line[0].split('.')[0]
		#print(line)
		phage_def = line[1]
		if phage_id not in phage_inf_dict.keys():
			phage_inf_dict.update({phage_id:phage_def.strip('"')})

	with open(file) as f:
		contents = f.readlines()
	hit_dict = {}
	f_result = open(resufile,'w')
	f_result.write('bac_id\thit_phage_id\thit_phage_def\tspacer\tspacer_length\thit_spacer_start\thit_spacer_end\thit_phage_start\thit_phage_end\tidentity\tcoverage\tmismatch\tevalue\tbit-score\n')
	f_result.flush()
	for line in contents:
		line = line.strip().split('\t')
		#print(line)
		if '|' in line[1]:
			phage_id = line[1].split('|')[1].split('.')[0]
		else:
			phage_id = line[1].split('.')[0]
		if phage_id not in phage_inf_dict.keys():
			continue
		phage_def = phage_inf_dict[phage_id]
		spacer_name = line[0]
		spacer_length = float(spacer_name.split('|')[2])
		identity = float(line[2])
		gap = int(line[5])
		mismatch = int(line[4])
		bit_score = float(line[-1])
		query_start = int(line[-6])
		query_end = int(line[-5])
		hit_start = int(line[-4])
		hit_end = int(line[-3])
		hit_length = abs(query_start-query_end)+1
		coverage = float(hit_length)/spacer_length
		bac_id = spacer_name.split('|')[-2]
		evalue = float(line[-2])
		if (int(mismatch) <= int(filter_mismatch)) and (coverage>=float(filter_coverage)) and (gap<=1):
			if bac_id not in hit_dict.keys():
				hit_dict.update({bac_id:{}})
			if phage_id not in hit_dict[bac_id].keys():
				hit_dict[bac_id].update({phage_id:[]})
			resuline = list(map(str,[bac_id,phage_id,phage_def,spacer_name,spacer_length,str(query_start),str(query_end),str(hit_start),str(hit_end),identity,coverage,mismatch,evalue,bit_score]))			
			f_result.write('\t'.join(resuline)+'\n')
			f_result.flush()
	f_result.close()

def cal_mismatch(query_sequence,hit_sequence,direction):
	trans_dict = {'A':'T','G':'C','C':'G','T':'A','N':'N'}
	if direction=='same':
		hmm_distance = Levenshtein.hamming(query_sequence,hit_sequence)
	else:
		hit_sequence = [trans_dict[base] for base in hit_sequence]
		hit_sequence = str(''.join(hit_sequence[::-1]))
		hmm_distance = Levenshtein.hamming(query_sequence,hit_sequence)
	return hmm_distance

def multialgn(prog_dir,input_file,output_file):
	cmd_mafft = '%s/software/mafft --clustalout --reorder --adjustdirection --thread 10 "%s" > "%s"'%(prog_dir,input_file,output_file)
	os.system(cmd_mafft)

def blastdbcmd(prog_dir,db_path,query_id):
	command = "%s/software/blastdbcmd -db %s -entry %s"%(prog_dir,db_path,query_id)
	cmd_search = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output,error = cmd_search.communicate()
	output = ''.join(str(output).strip("'b").split('\\n')[1:]).strip()
	return output

def spacer_match_virus(prog_dir,spacer_file,spacer_parse_file,save_prefix):
	save_file = save_prefix+"_hit_result.txt"
	f_save = open(save_file,'w')
	f_save.write('bac_id\tspacer_id\thit_phage_id\thit_phage_def\tspacer_sequence\thit_phage_region\thit_phage_sequence\tmismatch\tcoverage\thmm_mismatch\n')	
	
	#virus information
	#virus_file = "/zrom1/ganrui/crisprimmunity/data/virus/ncbi_bacteriophage_info_100bp.txt"
	
	virus_info_file = "/zrom1/ganrui/crisprimmunity/data/phage_plasmid/ncbi_bacteriophage_PLSDB_plasmid_info.txt"
	virus_db_path = "/zrom1/ganrui/crisprimmunity/data/database/phage_plasmid_nucl/phage_plasmid_nucl_db"
	with open(virus_info_file) as f:
		contents = f.readlines()
	virus_info_dict = {}
	for line in contents[1:]:
		line = line.strip().split('\t')
		virus_id = line[0]
		virus_def = line[1]
		virus_info_dict.update({virus_id.split('.')[0]:[virus_id,virus_def]})	
	#print(virus_info_dict)
	#spacer sequence information
	spacer_dict = {}
	with open(spacer_file) as f:
		contents = f.read().strip()
	for line in contents.split('>')[1:]:
		spacer_id = line.split('\n')[0].strip()
		spacer_sequence = line.split('\n')[1].strip()
		spacer_dict.update({spacer_id:spacer_sequence})
	
	#parse information
	header = 0
	with open(spacer_parse_file) as f:
		contents = f.readlines()
	header = contents[0].strip().split('\t')
	for line in contents[1:]:
		line = line.strip().split('\t')
		spacer_id = line[header.index('spacer')]
		bac_id = line[header.index('bac_id')]		
		spacer_length = int(spacer_id.split('|')[-3])
		virus_id = line[header.index('hit_phage_id')]
		query_start = min(int(line[header.index('hit_spacer_start')]),int(line[header.index('hit_spacer_end')]))
		query_end = max(int(line[header.index('hit_spacer_start')]),int(line[header.index('hit_spacer_end')]))	
		hit_phage_start = min(int(line[header.index('hit_phage_start')]),int(line[header.index('hit_phage_end')]))
		hit_phage_end = max(int(line[header.index('hit_phage_start')]),int(line[header.index('hit_phage_end')]))
		hit_length = abs(int(query_start)-int(query_end))+1
		mismatch = line[header.index('mismatch')]
		identity = line[header.index('identity')]
		coverage = line[header.index('coverage')]
		

		if virus_id not in virus_info_dict.keys():
			continue

		virus_sequence = blastdbcmd(prog_dir,virus_db_path,virus_info_dict[virus_id][0])

		c_spacer_sequence = spacer_dict[spacer_id]
		
		if int(line[header.index('hit_phage_start')])<int(line[header.index('hit_phage_end')]):
			hit_direction = '+'
		else:
			hit_direction = '-'
		
		if int(line[header.index('hit_spacer_start')])<int(line[header.index('hit_spacer_end')]):
			query_direction = '+'
		else:
			query_direction = '-'
		
		if query_direction==hit_direction:
			direction = 'same'
			cur_subject_from = max(0,int(hit_phage_start)-(int(query_start)-1)-1)
			cur_subject_to = min(len(virus_sequence),(int(hit_phage_end)+(len(c_spacer_sequence)-int(query_end))))

		else:
			direction = 'different'	
			cur_subject_from = max(0,int(hit_phage_start)-(len(c_spacer_sequence)-int(query_end))-1)
			cur_subject_to = min(len(virus_sequence),(int(hit_phage_end)+(int(query_start)-1)))
		
		phage_def = virus_info_dict[virus_id][1]

		c_align_spacer_phage_file = save_prefix+'_'+spacer_id.replace('|','-')+'_'+virus_id+'_'+str(cur_subject_from+1)+'-'+str(cur_subject_to)+'.fa'
		with open(c_align_spacer_phage_file,'w') as f:
			f.write('>'+spacer_id+'\n'+c_spacer_sequence+'\n'+'>'+virus_id+'|'+str(cur_subject_from+1)+'-'+str(cur_subject_to)+'\n'+expand_virus_sequence)
		c_align_spacer_mutalign_file = save_prefix+'_'+spacer_id.replace('|','-')+'_'+virus_id+'_'+str(cur_subject_from+1)+'-'+str(cur_subject_to)+'_multialign'

		multialgn(prog_dir,c_align_spacer_phage_file,c_align_spacer_mutalign_file)
		
		with open(c_align_spacer_mutalign_file) as f:
			multialign_sequence = '\n'.join(f.read().split('\n')[3:]).strip()
		match_num = multialign_sequence.split('\n')[2].count('*')
		mismatch_num = int(spacer_length)-match_num
		f_save.write(bac_id+'\t'+spacer_id+'\t'+line[1]+'\t'+phage_def+'\t'+c_spacer_sequence+'\t'+str(cur_subject_from+1)+'-'+str(cur_subject_to)+'\t'+expand_virus_sequence+'\t'+str(mismatch)+'\t'+str(coverage)+'\t'+str(mismatch_num)+'\n')
		f_save.flush()
	f_save.close()

def mkdir(dirname):
	command = "mkdir -p "+dirname
	os.system(command)

def predict(prog_dir,spacer_file,outdir,mismatch,coverage):
	#step1 : spacer blastn phage database
	blastn_file = os.path.join(outdir,'spacer_blastn_phage.txt')
	spacer_blastn_phage(prog_dir,spacer_file,blastn_file,6)
	
	#step1:filter mismatch and coverage
	parse_file = os.path.join(outdir,'spacer_blastn_phage_filter.txt')
	parse_spacer_blastn_phage_6(prog_dir,blastn_file,parse_file,mismatch,coverage)

	#step3:recalculate mismatch and multiple alignment
	save_dir = os.path.join(outdir,'spacer_target_phage')
	mkdir(save_dir)
	save_prefix = os.path.join(save_dir,'spacer_target')
	spacer_match_virus(prog_dir,spacer_file,parse_file,save_prefix)

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',help='''Path of the spacer file''')
	parser.add_argument('--output',help='Path of the output directory.\n')
	parser.add_argument('--mismatch',help='Path of the output directory.\n')
	parser.add_argument('--coverage',help='Path of the output directory.\n')
	parser.add_argument('--prog_dir', help='directory of program.\n')

	args = parser.parse_args()
	if args.prog_dir:
		prog_dir = args.prog_dir
	if args.input:
		input_file = args.input
	if args.output:
		outdir = args.output
	if args.mismatch:
		mismatch = args.mismatch
	else:
		mismatch = 2
	if args.coverage:
		coverage = args.coverage
	else:
		coverage = 0.9
	predict(prog_dir,input_file,outdir,mismatch,coverage)