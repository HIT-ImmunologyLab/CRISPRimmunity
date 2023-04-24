import os
import re,argparse
from Bio import Entrez
from Bio import SeqIO
import xml.etree.ElementTree as ET
import threading
import Levenshtein
from lxml import etree
def getFaaFromGB(input_file,faa_file):   #parse protein from genbank files in phaster
	special_pros = ['capsid','head','plate','tail','coat','portal','holin','integrase','transposase','terminase','protease','lysis','bacteriocin','tRNA']
	records = SeqIO.parse(input_file, "gb")
	counter = 0
	savefile = open(faa_file, 'w')
	for record in records:
		contig_id = record.id
		c_contig_length = len(str(record.seq))
		for feature in record.features:
			if feature.type == 'CDS':                
				location = feature.location
				if str(location).find('+') != -1:
					direction = '+'
				elif str(location).find('-') != -1:
					direction = '-'
				if '<' in str(location):
					location = str(location).replace('<','')
				if '>' in str(location):
					location = str(location).replace('>','')
				locations = re.findall("\d+\.?\d*",str(location))
				min_start = locations[0]
				max_end = locations[1]
				for loc in locations:
					if int(loc)<int(min_start):
						min_start = loc
					if int(loc)>int(max_end):
						max_end = loc
				if len(locations)>=4:
					if min(int(locations[0]),int(locations[1])) < min(int(locations[2]),int(locations[3])):
						if abs(max(int(locations[0]),int(locations[1]))-min(int(locations[2]),int(locations[3])))>1000:
							continue
					if min(int(locations[0]),int(locations[1])) > min(int(locations[2]),int(locations[3])):
						if abs(min(int(locations[0]),int(locations[1]))-max(int(locations[2]),int(locations[3])))>1000:
							continue				
				location = str(int(min_start)+1)+'_'+max_end+'_'+direction
				counter = counter+1
				if 'product' in feature.qualifiers:
					product = feature.qualifiers['product'][0]	  
				else:
					product = 'unknown'
				product = product.replace(' ','-').replace('>','').replace('<','').replace('|','-')
				if 'protein_id' in feature.qualifiers:
					proteinId = feature.qualifiers['protein_id'][0]
				else:
					if 'inference' in feature.qualifiers:
						strInference = str(feature.qualifiers['inference'])
						if 'RefSeq' in strInference:
							proteinId = strInference.split('RefSeq:')[1].rstrip(']').rstrip('\'')
						elif 'SwissProt' in strInference:
							proteinId = strInference.split('SwissProt:')[1].rstrip(']').rstrip('\'')
						else:
							proteinId = 'unknown'
					else:
						proteinId = 'unknown'
				if 'translation' in feature.qualifiers:
					translation = feature.qualifiers['translation'][0]
					savefile.write('>' +contig_id+ '|' + str(proteinId)+'|' + str(location)+'|'+product+'\n')					
					if translation[-1] == '\n':
						savefile.write(translation)
					else:
						savefile.write(translation + '\n')
	savefile.close()

def pred_orf(fasta_file,faa_prefix):
	cmd_fragGeneScan = '/home/qiusuo/FragGeneScan/run_FragGeneScan.pl -genome %s -out %s -complete=1 -train=complete -thread=20' % (fasta_file, faa_prefix)
	os.system(cmd_fragGeneScan)

def diamond_blastp(prog_dir,file,outfile,database,format,evalue,identity,coverage):
	num_threads = 20
	script = prog_dir + "/software/diamond blastp -d "+database+" -q "+file+" -f "+str(format)+" -e "+str(evalue)+" -o "+outfile+" -p "+str(num_threads)+" --id "+str(identity)+" -k 1"
	os.system(script)

def find_known_anti(prog_dir,protein_file,database,outdir,identity,coverage):
	evalue = 0.01
	
	out_blastp_file = os.path.join(outdir,'blastp_known_anti.txt')
	diamond_blastp(prog_dir,protein_file,out_blastp_file,database,5,evalue,identity,coverage)
	
	predict_file = os.path.join(outdir,'predict_known_anti_info.txt')
	blastp_add_info_file = os.path.join(outdir,'blastp_known_anti_add_info.txt')
	ParseResult(out_blastp_file,blastp_add_info_file,predict_file,identity,coverage)

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

def get_inf(file,outdir):
	with open(file) as f:
		contents = f.read().strip()
	if contents[0]=='>':
		type = 'fasta'
	else:
		type = "genbank"
	return type

def find_known_aca(prot_dir,protein_file,database,outdir,identity,coverage):
	evalue = 0.01	
	out_blastp_file = os.path.join(outdir,'blastp_known_aca.txt')
	diamond_blastp(prot_dir,protein_file,out_blastp_file,database,5,evalue,identity,coverage)
	
	predict_file = os.path.join(outdir,'predict_known_aca.txt')
	blastp_add_info_file = os.path.join(outdir,'blastp_known_aca_add_info.txt')
	ParseResult(out_blastp_file,blastp_add_info_file,predict_file,identity,coverage)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', help='absPath of genome file.\n')
	parser.add_argument('--output', help='the out dir.\n')
	parser.add_argument('--protein', help='the out dir.\n')
	parser.add_argument('--identity', help='the out dir.\n')
	parser.add_argument('--coverage', help='the out dir.\n')
	parser.add_argument('--prog_dir', help='directory of program.\n')

	args = parser.parse_args()
	if args.input:
		input_file = args.input
	if args.output:
		outdir = args.output
	if args.prog_dir:
		prog_dir = args.prog_dir
	if args.identity:
		identity = args.identity
	else:
		identity = 0.4
	if args.coverage:
		coverage = args.coverage
	else:
		coverage = 0.7
	if args.protein:
		protein_file = args.protein
	else:
		type = get_inf(input_file,outdir)
		protein_prefix = os.path.join(outdir,'protein')
		protein_file = os.path.join(outdir,'protein.faa')
		if type=='genbank':
			getFaaFromGB(input_file,protein_file)
		else:
			pred_orf(input_file,protein_prefix)
	anti_database = "%s/database/anit_verified_db/anti-crispr_db_non-redundant_verified_db"%(prog_dir)
	aca_database = "%s/database/anit_verified_db/Acr_dimaond_db"%(prog_dir)
	m1 = threading.Thread(target=find_known_anti,args=(prog_dir,protein_file,anti_database,outdir,identity,coverage))
	m2 = threading.Thread(target=find_known_aca,args=(prog_dir,protein_file,aca_database,outdir,identity,coverage))
	m1.start()
	m2.start()
	m1.join()
	m2.join()