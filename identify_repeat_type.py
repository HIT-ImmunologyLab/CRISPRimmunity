import re
import subprocess,os,sys
import argparse
Repeat_families_to_types = {                
                                "F1": ["I-B", "III-A", "III-B"],
                                "F2": ["I-E"],
                                "F3": ["I-C"], 
                                "F4": ["I-C", "I-E", "II-B"],
                                "F5": ["I-F"], 
                                "F6": ["I-A"], 
                                "F7": ["I-A"], 
                                "F8": ["I-F"], 
                                "F9": ["III-B"], 
                                "F10": ["I-B", "III-B"], 
                                "F11": ["III-B"], 
                                "F12": ["II-B", "III-A"], 
                                "F13": ["I-A", "III-B"], 
                                "F14": ["I-A", "I-D", "III-A"],
                                "F15": ["I-A", "III-B"], 
                                "F16": ["III-A"], 
                                "F17": ["?"], 
                                "F18": ["I-E", "II-B"], 
                                "F19": ["?"], 
                                "F20": ["I-B"], 
                                "F21": ["I-E"], 
                                "F22": ["I-E"], 
                                "F23": ["I-D", "II-B"], 
                                "F24": ["III-A", "III-B"], 
                                "F25": ["I-A", "II-B", "III-A"],
                                "F26": ["?"], 
                                "F27": ["II-A", "II-B"], 
                                "F28": ["I-A"], 
                                "F29": ["III-A"], 
                                "F30": ["?"], 
                                "F31": ["III-A"], 
                                "F32": ["I-C"], 
                                "F33": ["I-C", "I-E", "II-B"],
                                "F34": ["II-B"], 
                                "F35": ["II-A"], 
                                "F36": ["?"], 
                                "F37": ["I-C", "III-B"], 
                                "F38": ["I-A", "III-B"], 
                                "F39": ["I-A", "I-B", "II-B"],
                                "F40": ["I-B"],
                                "Cas12a_repeats": ["V-A"],         
                                "Cas12b_repeats": ["V-B"], 
                                "Cas12c_repeats": ["V-C"], 
                                "Cas12d_repeats": ["V-D"], 
                                "Cas12e_repeats": ["V-E"], 
                                "Cas13a_repeats": ["VI-A"],    
                                "Cas13b1_repeats": ["VI-B1"],
                                "Cas13b2_repeats": ["VI-B2"],
                                "Cas13c_repeats": ["VI-C"],
                                "c2c8_V-U2_repeats": ["V-U2"], 
                                "c2c10_V-U3_repeats": ["V-U3"], 
                                "c2c9_V-U4_repeats": ["V-U4"],
                                "c2c5_V-U5_repeats": ["V-U5"],                                
                                 }
def repeat_HMM_check(prog_dir,consensus_repeat,outfile):
	repeat_HMM_file = "%s/database/repeat_hmm/REPEATS_HMMs.hmm"%(prog_dir)
	hmm_cmd = "%s/software/nhmmscan -E 1e-6 --noali %s %s"%(prog_dir,repeat_HMM_file,consensus_repeat)
	handle = subprocess.Popen(hmm_cmd.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output,error = handle.communicate()
	output = bytes.decode(output)
	with open(consensus_repeat) as f:
		contents = f.read().strip().split('>')
	repeats = {}
	for repeat in contents[1:]:
		title = repeat.strip().split('\n')[0].strip()
		sequence = repeat.strip().split('\n')[-1].strip()
		repeats.update({title:sequence})
	f_result = open(outfile,'w')
	f_result.write('repeat\tsequence\trepeat_direction\tType_repeat\tpossible types\n')
	# if error != "":
	# 	print(error)
	# 	sys.exit()
	Type_repeat = "Repeat not recognized"; repeat_direction = 0; possible_types = []
	expression1 = re.compile("\s+E-value\s+score\s+bias")   #find the table labels for the data
	for repeat in output.split('Query:')[1:]:		
		title = repeat.split('\n')[0].strip().split()[0]
		sequence = repeats[title]
		flag = 0
		for line in repeat.split('\n'):
			a = re.search(expression1, line)
			if a is not None: 
				i = 2
				e_value = 1
				while True:
					data_string = repeat.split('\n')[repeat.split('\n').index(line) + i]  #Get the data from lines below
					i += 1
					if data_string == "": #No hits found
						break
					try:
						found_e_value = float(data_string.strip().split()[0])
						if found_e_value < e_value:  #want to take the lowest e-value (best match)
							e_value = found_e_value
							repeat_group = data_string.strip().split()[3]  #The repeat group is the 4th position
							repeat_direction = int(data_string.strip().split()[5]) - int(data_string.strip().split()[4]) #positive means CRT and HMM direction match
						#Use the identified group to guess at the Type
							if repeat_group[-1] == 'R':  #Removes the R designation that was added to indicate where the original REPEATS data was backward
								repeat_group = repeat_group[:-1]
							possible_types = Repeat_families_to_types[repeat_group]
							if len(possible_types) == 1:
								Type_repeat = "Type {0}".format(possible_types[0])
							elif len(possible_types) == 2:
								Type_repeat = "Type {0} or {1}".format(possible_types[0],possible_types[1])
							elif len(possible_types) > 2:
								Type_repeat = "Type" + "".join([" {0},".format(string) for string in possible_types[:-1]]) + " or {0}".format(possible_types[-1])
							Type_repeat += " (group {0})".format(repeat_group)
							print((repeat_direction,Type_repeat,possible_types))
							flag = 1
							f_result.write(title+'\t'+sequence+'\t'+str(repeat_direction)+'\t'+Type_repeat+'\t'+','.join(possible_types)+'\n')
							f_result.flush()
					except:
						break
				break
		if flag==0:
			f_result.write(title+'\t'+sequence+'\t'+'NA'+'\t'+'NA'+'\t'+'NA'+'\n')
	return (repeat_direction,Type_repeat,possible_types)
	f_result.close()

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', help='absPath of repeat fasta file.\n')
	parser.add_argument('--output', help='out file of parsed repeat file.\n')
	parser.add_argument('--prog_dir', help='directory of program.\n')
    
	args = parser.parse_args()
	if args.input:
		repeat_file = args.input
	if args.output:
		outfile = args.output
	if args.prog_dir:
		prog_dir = args.prog_dir

	repeat_HMM_check(prog_dir,repeat_file,outfile)