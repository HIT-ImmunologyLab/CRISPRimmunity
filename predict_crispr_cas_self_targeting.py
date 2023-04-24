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
import numpy as np
import multiprocessing

CRISPR_types = {
                "Type I-A":[["Cas6"],["Csa5"],["Cas7"],["Cas5"],["Cas8a1","Cas8a2","Csa4","Csx9"],["Cas3'"],["Cas3''"],["Cas2"],["Cas4","Csa1"],["Cas1"],["Cas4"]],
                "Type I-B":[["Cas6"],["Cas8b1","Csh1"],["Cas7","Cst2"],["Cas5"],["Cas3"],["Cas4","Csa1"],["Cas1"],["Cas2"]],
                "Type I-C":[["Cas3"],["Cas5"],["Cas8c","Csd1","Csp2"],["Cas7"],["Cas4","Csa1"],["Cas1"],["Cas2"]],
                "Type I-U":[["Cas3"],["Cas8c"],["Cas7"],["Cas5","GSU0054"],["Cas6"],["Cas4","Csa1"],["Cas1"],["Cas2"]],
                "Type I-D":[["Cas3'"],["Cas3''"],["Cas10d","Csc3"],["Cas7","Csc2"],["Cas5","Csc1"],["Cas6"],["Cas4"],["Cas1"],["Cas2"]],
                "Type I-E":[["Cas3"],["Cas8e","Cse1"],["Cse2","CasB"],["Cas7","Cse4","CasC"],["Cas5","CasD"],["Cas6","Cse3","CasE"],["Cas1"],["Cas2"]],
                "Type I-F":[["Cas1"],["Cas2","Cas3"],["Cas8f","Csy1"],["Cas5","Csy2"],["Cas7","Csy3"],["Cas6","Cas6f","Csy4"]],
                "Type I-G":[["Cas6"],["Cas8a1","Cst1"],["Cas7","Cst2"],["Cas5","Cas5t"],["Cas3"],["Cas4"],["Cas1"],["Cas2"]],
                
                "Type II-A":[["Cas9","Csn1"],["Cas1"],["Cas2"],["Csn2"]],
                "Type II-B":[["Cas9","Csn1"],["Cas1"],["Cas2"],["Cas4","Csa1"],],
                "Type II-C":[["Cas9","Csn1"],["Cas1"],["Cas2"]],
                
                "Type III-A":[["Cas6"],["Cas10","Csm1"],["Csm2"],["Cas7","Csm3"],["Cas5","Csm4"],["Cas7","Csm5"],["Csm6"],["Cas1"],["Cas2"]],
                "Type III-B":[["Cas7","Cmr1"],["Cas10"],["Cas5","Cmr4"],["Cmr5"],["Cas6"],["Cas7","Cmr6"],["Cas1"],["Cas2"]],
                "Type III-C":[["Cas7","Cmr1"],["Cas7","Cmr6"],["Cas10"],["Cas7","Cmr4"],["Cmr5"],["Cas5","Cmr3"]],
                "Type III-D":[["Cas10"],["Cas7","Csm3"],["Cas5","Csx10"],["Csm2"],["Cas7","Csm3"],["Cas7","Csm3"],["all1473"],["Cas7","Csm3"]],                
                "Type IV-A":[["dinG","Csf4"],["Csf1"],["Cas7","Csf2"],["Cas5","Csf3"]],                
                "Type V-A":[["Cas12a","Cpf1"],["Cas4","Csa1"],["Cas1"],["Cas2"]],
                "Type V-B":[["Cas12b","c2c1"],["Cas1"],["Cas2"]],
                "Type V-C":[["Cas12c","c2c3"],["Cas1"]],
                "Type V-D":[["Cas12d","CasY"],["Cas1"]],
                "Type V-E":[["Cas12e","CasX"],["Cas4"],["Cas1"],["Cas2"]],
                "Type V-U1":[["c2c4"]],
                "Type V-U2":[["c2c8"]],
                "Type V-U3":[["c2c10"]],
                "Type V-U4":[["c2c9"]],
                "Type V-U5":[["c2c5"]],
                
                "Type VI-A":[["Cas13a","c2c2"],["Cas1"],["Cas2"]],
                "Type VI-B1":[["Cas13b1","c2c6"],["Csx27"]],
                "Type VI-B2":[["Cas13b2","c2c6"],["Csx28"]],
                "Type VI-C":[["Cas13c","c2c7"]],
                "Type VI-D":[["Cas13d","c2c7"]]
                
                }
                
                
Cas_proteins =  {
                "Csa1":["Type I-A","Type I-B","Type I-C","Type I-D","Type II-B","Type V-A","Type V-B","Type V-E",],
                "Csa4":["Type I-A"],
                "Csa5":["Type I-A"],
                "Csc1":["Type I-D"],
                "Csc2":["Type I-D"],
                "Csc3":["Type I-D"],
                "Csd1":["Type I-C"],
                "Cse1":["Type I-E"],
                "Cse2":["Type I-E"],
                "Cse3":["Type I-E"],
                "Cse4":["Type I-E"],
                "Csh1":["Type I-B"],
                "Csf1":["Type IV-A"],
                "Csf2":["Type IV-A"],
                "Csf3":["Type IV-A"],
                "Csf4":["Type IV-A"],
                "Csn2":["Type II-A"],
                "Csp2":["Type I-C"],
                "Csy1":["Type I-F"],
                "Csy2":["Type I-F"],
                "Csy3":["Type I-F"],
                "Csy4":["Type I-F"],
                "Csn1":["Type II-A","Type II-B","Type II-C"],
                "Csm1":["Type III-A"],
                "Csm2":["Type III-A","Type III-D"],
                "Csm3":["Type III-A","Type III-D"],
                "Csm4":["Type III-A"],
                "Csm5":["Type III-A"],
                "Csm6":["Type III-A"],
                "Cst1":["Type I-G"],
                "Cst2":["Type I-G","Type I-B"],
                "Csx9":["Type I-A"],
                "Csx10":["Type III-D"],
                "Cmr1":["Type III-B","Type III-C"],
                "Cmr3":["Type III-C"],
                "Cmr4":["Type III-B","Type III-C"],
                "Cmr5":["Type III-B","Type III-C"],
                "Cmr6":["Type III-B","Type III-C"],
                "GSU0054":["Type I-U"],
                "all1473":["Type III-D"],
                "dinG":["Type IV-A"],
                "Cas1":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type II-A","Type II-B","Type II-C","Type III-A","Type III-B","Type V-A","Type V-B","Type V-C","Type V-D","Type V-E","Type VI-A"],
                "Cas2":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type II-A","Type II-B","Type II-C","Type III-A","Type III-B","Type V-A","Type V-B","Type V-E","Type VI-A"],
                "Cas3":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G"],                    
                "Cas3'":["Type I-A","Type I-D"],
                "Cas3''":["Type I-A","Type I-D"],
                "Cas4":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-G","Type II-B","Type V-A","Type V-B"],               
                "Cas5":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type III-A","Type III-B","Type III-C","Type III-D","Type IV-A"],                              
                "Cas5t":["Type I-G"],
                "Cas6":["Type I-A","Type I-B","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type III-A","Type III-B"],                                             
                "Cas6f":["Type I-F"],
                "Cas7":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type III-A","Type III-B","Type III-C","Type III-D","Type IV-A"],
                "Cas8a1":["Type I-A","Type I-G"],   
                "Cas8a2":["Type I-A"],   
                "Cas8b1":["Type I-B"],                                                                                
                "Cas8c": ["Type I-C","Type I-U"],                                                                                                                                           
                "Cas8e": ["Type I-E"],                                                                                                                                                                                                       
                "Cas8f": ["Type I-F"],                                                                                                                                                                                                                                                                   
                "Cas9":["Type II-A","Type II-B","Type II-C"],                                                                                                                                                                                                                                                                                                                
                "Cas10":["Type III-A","Type III-B","Type III-C","Type III-D"],
                "Cas10d":["Type I-D"],
                "CasB":["Type I-E"],
                "CasD":["Type I-E"],
                "CasC":["Type I-E"],
                "CasE":["Type I-E"],
                "Cas12a":["Type V-A"],
                "Cpf1":["Type V-A"],
                "Cas12b":["Type V-B"],
                "c2c1":["Type V-B"],
                "Cas12c":["Type V-C"],
                "c2c3":["Type V-C"],
                "Cas12d":["Type V-D"],
                "CasY":["Type V-D"],
                "Cas12e":["Type V-E"],
                "CasX":["Type V-E"],
                "c2c4":["Type V-U1"],
                "c2c5":["Type V-U5"],
                "c2c8":["Type V-U2"], 
                "c2c10":["Type V-U3"], 
                "c2c9":["Type V-U4"],  
                "Cas13a":["Type VI-A"],
                "c2c2":["Type VI-A"],
                "Cas13b1":["Type VI-B1"],
                "Cas13b2":["Type VI-B2"],
                "c2c6":["Type VI-B1","Type VI-B2"],
                "Csx27":["Type VI-B1"],
                "Csx28":["Type VI-B2"],
                "Cas13c":["Type VI-C"],
                "c2c7":["Type VI-C"],            
                "csa3":["Type I-A"],
                "RT":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type III-A","Type III-B","Type III-C","Type III-D"],
                "WYL":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type III-A","Type III-B","Type III-C","Type III-D","Type VI-D"],
                "PD-DExK":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type III-A","Type III-B","Type III-C","Type III-D"],                                                                                                                                                                                                                                                                                       
                "DEDDh":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G"],                                                                                                                                                                                                                                                                    
                "PrimPol":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type III-A","Type III-B","Type III-C","Type III-D"],
              } 

def is_known_Cas_protein(product,types_list=[]):   
    protein_name = ''
    is_Cas = False
    for key,values in Cas_proteins.items():
        flag = 0
        if 'c2c10' in product.lower():
            if 'c2c10' in key.lower():
                flag = 1
        else:
            if product.lower().find(key.lower()) > -1:
                flag = 1
        if flag ==1:
            #Check to see if the protein type is annotated
            protein_name = key
            parts = [x.lower() for x in product.split(" ")]
            for part in parts:
                if "type" in part and part != parts[-1]:   #ignore things that end in type
                    type_expected = "Type " + parts[parts.index(part) + 1].upper()  #find type and look next to it for the designation
                    #If a type pops up, but the letter is not there, add all possibilities
                    for value in values:
                        if type_expected in value:
                            types_list.append(value)  #Adds weight that the proper type will be identified
                    break
            types_list += values   #the correct Type will have the most entries in this list
            is_Cas = True
            break 
    return is_Cas,protein_name,types_list    #Returns the CRISPR type if detected and whether the locus appears complete

def Type_check(types_list):   
    #Now determine what type it is and if it is complete
    #Determine the type by counting the Type that is the most prevalent
    data = Counter(types_list)   
    Type = ""
    off_count = 0
    for most_common in data.most_common():
        count = most_common[1]
        if count >= off_count:
            if Type == '':
                Type = most_common[0]
            elif Type.count(",") == 1:
                Type += ", or " + most_common[0]
            else:
                Type += ", " + most_common[0]
            off_count = count
        if count < off_count or most_common == data.most_common()[-1]:  #reached end of list or found less common entry 
            if Type.count(",") > 2:  #Too many possibilities to exactly identify
                Type = "Unclear"
            elif Type.count(",") == 1:
                Type.replace(",", " or")
                Type += "?"
            elif Type.count(",") == 2:
                Type[:Type.rfind(",")] + " or" + Type[Type.rfind(",")+1:]
                Type += "?"
            break
    if Type == "":
        Type = "Orphan"                
    return Type

def check_crispr_type_by_cas_prot(cas_prots):
    types_list = []
    type_dict = {"Cas8a1":["Type I-A","Type I-G"],   
                "Cas8a2":["Type I-A"],   
                "Cas8b1":["Type I-B"],                                                                                
                "Cas8c": ["Type I-C","Type I-U"],                                                                                                                                           
                "Cas8e": ["Type I-E"],                                                                                                                                                                                                       
                "Cas8f": ["Type I-F"],                                                                                                                                                                                                                                                                   
                "Cas9":["Type II-A","Type II-B","Type II-C"],                                                                                                                                                                                                                                                                                                                
                "Cas10":["Type III-A","Type III-B","Type III-C","Type III-D"],
                "Cas10d":["Type I-D"],
                "CasB":["Type I-E"],
                "CasD":["Type I-E"],
                "CasC":["Type I-E"],
                "CasE":["Type I-E"],
                "Cas12a":["Type V-A"],
                "Cpf1":["Type V-A"],
                "Cas12b":["Type V-B"],
                "c2c1":["Type V-B"],
                "Cas12c":["Type V-C"],
                "c2c3":["Type V-C"],
                "Cas12d":["Type V-D"],
                "CasY":["Type V-D"],
                "Cas12e":["Type V-E"],
                "CasX":["Type V-E"],
                "c2c4":["Type V-U1"],
                "c2c5":["Type V-U5"],
                "c2c8":["Type V-U2"], 
                "c2c10":["Type V-U3"], 
                "c2c9":["Type V-U4"],  
                "Cas13a":["Type VI-A"],
                "c2c2":["Type VI-A"],
                "Cas13b1":["Type VI-B1"],
                "Cas13b2":["Type VI-B2"],
                "c2c6":["Type VI-B1","Type VI-B2"],
                "Csx27":["Type VI-B1"],
                "Csx28":["Type VI-B2"],
                "Cas13c":["Type VI-C"],
                "c2c7":["Type VI-C"]}
    for cas_name in cas_prots:    
        is_Cas,protein_name,types_list = is_known_Cas_protein(cas_name,types_list)
    crispr_type = Type_check(types_list)
    for cas_prot in cas_prots:
        for c_cas_prot in type_dict.keys():
            if cas_prot!='c2c1':
                if cas_prot.lower() == c_cas_prot.lower():
                    crispr_type = crispr_type+','+','.join(type_dict[c_cas_prot])
                else:
                    if c_cas_prot in ['Cas8a1','Cas8a2','Cas8b1','Cas13b1','Cas13b2']:
                        if cas_prot.lower() == c_cas_prot.lower()[0:-1]:
                            crispr_type = crispr_type+','+','.join(type_dict[c_cas_prot])
            else:
                crispr_type = crispr_type+',Type V-B'
    
    crispr_type = ','.join(list(set(crispr_type.split(','))))
    if 'Unclear' in crispr_type and ('Type' in crispr_type):
        crispr_type = crispr_type.replace('Unclear','').strip(',')
    return crispr_type

def mkdir(dirname):
    if not os.path.exists(dirname):
        command = "mkdir -p "+dirname
        os.system(command)

def runPILERCR(prog_dir,inFilePath,outFilePath):
    cmd_pilercr = '%s/software/pilercr -minarray 1 -in %s -out %s -noinfo -quiet'%(prog_dir,inFilePath,outFilePath)
    os.system(cmd_pilercr)

def runCRISPRFinder(prog_dir,inFilePath,outDirPath):
    sudoPassword = 'dyi*AhJXV7'
    username = 'qiusuo'
    os.chdir('/'.join(outDirPath.split('/')[0:-1]))
    # cmd_casfinder = '/home/simon/CRISPRCasFinder/CRISPRCasFinder.pl -cf CasFinder-2.0.2 -def SubTyping -keep -fl 500 -i %s -out %s -so /home/simon/CRISPRCasFinder/sel392v2.so'%(inFilePath,outDirPath)
    cmd_casfinder = '%s/software/CRISPRCasFinder.pl -cf CasFinder-2.0.2 -def SubTyping -keep -fl 500 -i %s -out %s -so %s/library/sel392v2.so' % (prog_dir,
    inFilePath, outDirPath,prog_dir)
    ## if need sudo

    password_command = ""
    os.system('echo %s|sudo -S %s' % (sudoPassword, cmd_casfinder))
    os.system('wait')
    os.system('echo su %s' %(username))
    ## if not need sudo
    # os.system(cmd_casfinder)

def runCRT(prog_dir,inFilePath,outFilePath):
    repeats = 4
    min_repeat_length = 18
    max_repeat_length = 45
    min_spacer_length = 18
    max_spacer_length = 45
    CRT_params = [repeats, str(min_repeat_length), str(max_repeat_length), str(min_spacer_length), str(max_spacer_length)]
    # CRISPR_cmd = "java -cp /zrom/zfx/STSS/Self-Targeting-Spacer-Searcher/bin/CRT1.2-CLI.jar crt -minNR %s -minRL %s -maxRL %s -minSL %s -maxSL %s %s %s"%(CRT_params[0],CRT_params[1],CRT_params[2],CRT_params[3],CRT_params[4],inFilePath,outFilePath)
    CRISPR_cmd = "java -cp %s/software/CRT1.2-CLI.jar crt -minNR %s -minRL %s -maxRL %s -minSL %s -maxSL %s %s %s" % (
    prog_dir,CRT_params[0], CRT_params[1], CRT_params[2], CRT_params[3], CRT_params[4], inFilePath, outFilePath)
    crispr_search = subprocess.Popen(CRISPR_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output,error = crispr_search.communicate()

def parse_pildercr_result(pilercr_result_file,pilercr_spcfile,pilercr_cspfile):
    bfile_tmp = open(pilercr_result_file, 'r')
    lines = bfile_tmp.readlines()
    llen = len(lines)
    spc_num = 0
    spcfile = open(pilercr_spcfile, 'w')
    cspfile = open(pilercr_cspfile, 'w')
    # ------------------find if there are crisprs---------
    i = 0
    while 'putative CRISPR arrays found' not in lines[i]: i += 1
    p = re.compile(r'\d+ putative CRISPR arrays found\.')
    res = p.findall(lines[i])
    if len(res) != 0:
        if res[0][0] == '0':
            bfile_tmp.close()
            return str(0) + '\t' + str(spc_num) + '\n'
    # -------------------find all spacers----------------
    # spcfile = open(pilercr_spcfile, 'w')
    ifstart = False
    contigs = {}
    array_ids = {}
    while 'SUMMARY BY SIMILARITY' not in lines[i]:
        if 'Array' in lines[i]:
            original_array_id = lines[i].split()[-1].strip()
        if '>' in lines[i]:
            fragbact = lines[i][1:lines[i].find(' ')].split('.')[0].strip()
            if fragbact not in contigs.keys():
                array_cnt = 0
                contigs.update({fragbact:0})
        if '===' in lines[i]:
            array_cnt_spacer = 0
            ifstart = not ifstart
            if ifstart: array_cnt += 1
            i += 1
        if ifstart:
            array_cnt_spacer += 1
            tabs = lines[i].split()
            if len(tabs) == 7:
                #fragbact = tabs[1].split('.')[0]
                spc_num = spc_num + 1
                if original_array_id not in array_ids.keys():
                   array_ids.update({original_array_id:array_cnt}) 
                spcfile.write('>' + str(array_cnt) + '.' + str(array_cnt_spacer) + '|' + str(int(tabs[0]) + int(tabs[1])) + '|' + tabs[3] + '|' + fragbact + '\n')
                spcfile.write(tabs[6] + '\n')
        i += 1
    spcfile.close()
    # -------------------get the information about repeat spacers in crispr----------
    # cspfile = open(pilercr_cspfile, 'w')
    ifstart = False
    while 'SUMMARY BY POSITION' not in lines[i]:
        if '===' in lines[i] or '***' in lines[i]:
            ifstart = True
            i += 1
        # if len(lines[i].strip())==0:
        #     crispr_count = 0
        if ifstart:
            tabs = lines[i].strip().split()
            # print tabs
            if len(tabs) == 9:
                fragbact = tabs[1]
                for item in contigs.keys():
                    if fragbact.split('.')[0] in item:
                        fragbact = item
                        contigs[item] +=1
                        break
                #crispr_count = contigs[item]
                original_array_id = tabs[0].strip()
                crispr_count = array_ids[original_array_id]
                cspfile.write(
                    "%d\t%s\t%s\t%s\t%s\t%s\n" % (crispr_count, tabs[2], tabs[3], tabs[7], tabs[8].strip(), fragbact))
            if len(tabs) == 10:
                fragbact = tabs[1]
                for item in contigs.keys():
                    if fragbact.split('.')[0] in item:
                        fragbact = item
                        contigs[item] +=1
                        break
                #crispr_count = contigs[item]
                original_array_id = tabs[0].strip()
                crispr_count = array_ids[original_array_id]
                cspfile.write(
                    "%d\t%s\t%s\t%s\t%s\t%s\n" % (crispr_count, tabs[3], tabs[4], tabs[8], tabs[9].strip(), fragbact))
        i += 1
    cspfile.close()
    bfile_tmp.close()

def parse_casfinder_result(casfinder_result_file,casfinder_spcfile,casfinder_cspfile):
    jsonfile = open(casfinder_result_file)
    jsondata = json.load(jsonfile)
    spcfile = open(casfinder_spcfile, 'w')
    cspfile = open(casfinder_cspfile, 'w')
    # Parse Root
    Sequences = jsondata['Sequences']  # list
    for sequence in Sequences:
        # Parse every sequence information
        Id = sequence['Id']
        Crisprs = sequence['Crisprs']  # list
        Crisprs = sorted(Crisprs, key=lambda e: int(e.__getitem__('Start')))  # sort the crisprs
        for crispr in Crisprs:
            # Parse every crispr information in this sequence
            # output .tmp file here
            crisprStart = crispr['Start']  # integer
            crisprEnd = crispr['End']  # integer
            crisprDR_Consensus = crispr['DR_Consensus']
            number_Spacers = crispr['Spacers']  # integers
            crisprPotential_Orientation = crispr['Potential_Orientation']  # enum str : + - ND
            crisprRegions = crispr['Regions']  # list
            crisprid = crispr['Name'].split('_')[-1].strip()
            spcid = 1
            # output .csp file here
            cspfile.write("%s\t%s\t%d\t%s\t%s\t%s\t%s\n" % (
            crisprid, crisprStart, abs(int(crisprStart) - int(crisprEnd)), crisprPotential_Orientation,
            crisprDR_Consensus, Id,number_Spacers))
            for region in crisprRegions:
                # Parse every crispr region
                regionType = region['Type']  # enum str  : LeftFLANK+ (DR Spacer DR)+ RightFLANK+
                regionStart = region['Start']  # integer
                regionSequence = region['Sequence']
                # output .spc file here
                if regionType == "Spacer":
                    spcfile.write(">%s.%d|%s|%d|%s\n%s\n" % (
                    crisprid, spcid, regionStart, len(regionSequence), Id, regionSequence))
                    spcid = spcid +1
    spcfile.close()
    cspfile.close()

def get_num(string):
    strlist = re.findall("\d+\.?\d*",str(string))
    return strlist

def parse_spacer_crt(prog_dir,file,outdir):#get two files:crt.csp,crt.spc in outdir
    with open(file) as f:
        contents = f.read().strip()
    genome_id = contents.split('\n')[0].split(':')[1].strip().split()[0].strip()
    csp_file = os.path.join(outdir,'crt.csp')
    spc_file = os.path.join(outdir,'crt.spc')
    f_csp = open(csp_file,'w')
    # f_csp.write('crispr_index\tcrispr_start\tcrispr_length\trepeat_sequence\tgenome_id\n')
    f_spc = open(spc_file,'w')
    if 'REPEAT' in contents:
        spacers = contents.split('CRISPR')
        for spacer in spacers[1:]:
            spacer = spacer.strip()
            title = spacer.split('\n')[0].strip()
            title_nums = get_num(title)
            crispr_index = title_nums[0]
            crispr_length = abs(int(title_nums[-1])-int(title_nums[1]))
            #first line
            content = spacer.split('\n')[3].strip()
            crispr_start = get_num(content)[0]
            spacer_sequences = spacer.split('\n')[3:]
            flag = 1
            repeats = []
            for item in spacer_sequences:
                if ']' in item:
                    item = item.strip()
                    spacer_length = get_num(item)[-1]
                    spacer_sequence = item.split()[2]
                    repeat_sequence = item.split()[1]
                    repeats.append(repeat_sequence)
                    spacer_index = crispr_index+'.'+str(flag)
                    flag = flag+1
                    spacer_start = int(get_num(item)[0])+int(get_num(item)[1])
                    f_spc.write('>'+spacer_index+'|'+str(spacer_start)+'|'+str(spacer_length)+'|'+genome_id+'\n')
                    f_spc.write(spacer_sequence+'\n')
            repeats.append(spacer_sequences[flag-1].split()[1])
            i = 1
            fasta_string = ''
            valid_dna = "ACGTN"
            repeat_file = os.path.join(outdir,'crt_repeat')
            f_repeat = open(repeat_file,'w')
            for repeat in repeats:
                fasta_string += '>%s|%s\n%s\n'%(str(crispr_index),str(i),repeat)
                i += 1
            f_repeat.write(fasta_string)
            f_repeat.close()
            temp_file = os.path.join(outdir,'crt_repeat.fa')
            clustal_cmd = prog_dir + "/software/clustalo -i "+repeat_file+" --force --outfmt=clustal -o "+temp_file
  
            os.system(clustal_cmd)
            try:
                alignments = AlignIO.read(temp_file, "clustal")
                repeats_align = AlignInfo.SummaryInfo(alignments)
                consensus_repeat = str(repeats_align.dumb_consensus(ambiguous='N', require_multiple=1))
            except:
                consensus_repeat = repeats[0]
            f_csp.write(str(crispr_index)+'\t'+str(crispr_start)+'\t'+str(crispr_length)+'\t'+'unknown'+'\t'+consensus_repeat+'\t'+genome_id+'\n')
            f_csp.flush()
    f_csp.close()
    f_spc.close()

def run_crt_parse(prog_dir,contig_file,contig_outfile,contig_dir):
    runCRT(prog_dir,contig_file,contig_outfile)
    parse_spacer_crt(prog_dir,contig_outfile,contig_dir)

def multifasta_crt(prog_dir,strain_inf_file,contig_dir,outdir):
    with open(strain_inf_file) as f:
        contents = f.readlines()
    spacer_file = os.path.join(outdir,'crt.spc')
    crispr_file = os.path.join(outdir,'crt.csp')
    f_spc = open(spacer_file,'w')
    f_csp = open(crispr_file,'w')
    i=1
    pool = multiprocessing.Pool(processes=int(5))
    for line in contents[1:]:
        line = line.strip().split('\t')
        strain_id = line[0]
        i = i+1
        contig_file = os.path.join(contig_dir,strain_id+'.fa')
        c_contig_outdir = os.path.join(outdir,strain_id)
        mkdir(c_contig_outdir)
        c_contig_outfile = os.path.join(c_contig_outdir,'crt.out')        

        pool.apply_async(run_crt_parse,(prog_dir,contig_file,c_contig_outfile,c_contig_outdir))
        
    pool.close()
    pool.join()
    
    i=1
    for line in contents[1:]:
        line = line.strip().split('\t')
        strain_id = line[0]
        i = i+1
        contig_file = os.path.join(contig_dir,strain_id+'.fa')
        c_contig_outdir = os.path.join(outdir,strain_id)
        mkdir(c_contig_outdir)
        c_contig_outfile = os.path.join(c_contig_outdir,'crt.out')   
       
        try:
            contig_spc_file = os.path.join(c_contig_outdir,'crt.spc')
            contig_csp_file = os.path.join(c_contig_outdir,'crt.csp')
            with open(contig_spc_file) as f:
                line_resu = f.read().strip()
                if len(line_resu)>0:
                    f_spc.write(line_resu+'\n')
            with open(contig_csp_file) as f:
                line_resu = f.read().strip()
                if len(line_resu)>0:
                    f_csp.write(line_resu+'\n')
        except:
            pass
    f_spc.close()
    f_csp.close()

def inverse_complement(query_sequence):
    trans_dict = {'A':'T','G':'C','C':'G','T':'A','N':'N'}   
    hit_sequence = [trans_dict[base] for base in query_sequence]
    hit_sequence = str(''.join(hit_sequence[::-1]))
    return hit_sequence

def parse_crispridentify_result(crispridentify_dir,crispridentify_spcfile,crispridentify_cspfile):
    spcfile = open(crispridentify_spcfile, 'w')
    cspfile = open(crispridentify_cspfile, 'w')
    for contig_id in os.listdir(crispridentify_dir):
        c_dir = os.path.join(crispridentify_dir,contig_id)
        c_summary_file = os.path.join(c_dir,'Summary.csv')
        if not os.path.exists(c_summary_file):
            continue
        c_final_resu_file = os.path.join(c_dir,'Bona-Fide_Candidates.txt')
        if os.path.exists(c_final_resu_file):
            with open(c_final_resu_file) as f:
                crispr_contents = f.read()           
            contig_id = crispr_contents.split('\n')[0].strip().strip('>').split()[0]            
            with open(c_summary_file) as f:
                summary_contents = f.readlines()
            summary_header = summary_contents[0].strip().split(',')
            consensus_repeat_dict = {}
            for line in summary_contents[1:]:
                line = line.strip().split(',')
                category = line[summary_header.index('Category')]
                if category=='Bona-fide':
                    crispr_id = line[summary_header.index('ID')]
                    crispr_start = line[summary_header.index('Start')]
                    crispr_end = line[summary_header.index('End')]
                    crispr_length = line[summary_header.index('Length')]
                    repeat = line[summary_header.index('Consensus repeat')]
                    crispr_direction = line[summary_header.index('Strand')]
                    if crispr_direction=='Reversed':
                        crispr_direction = '-'
                    else:
                        crispr_direction = '+'
                    if crispr_id not in consensus_repeat_dict.keys():
                        consensus_repeat_dict.update({crispr_id:[repeat,crispr_direction]})
                    
                    
                    spacer_num = line[summary_header.index('Number of spacers')]
                    cspfile.write(crispr_id+'\t'+crispr_start+'\t'+crispr_length+'\t'+crispr_direction+'\t'+repeat+'\t'+contig_id+'\t'+spacer_num+'\n')
                    cspfile.flush()
           
            for crispr in crispr_contents.split('CRISPR:')[1:]:
                crispr_id = crispr.split('\n')[0].split(',')[0].strip()
                consensus_repeat = consensus_repeat_dict[crispr_id][0]
                crispr_direction = consensus_repeat_dict[crispr_id][1]
                i = 2
                spacer_id = 1
                while '_______' not in crispr.split('\n')[i]:
                    c_line = crispr.split('\n')[i].strip()
                    i = i +1
                    c_spacer_seq =  c_line.split()[-4].strip()
                    if '.' not in c_spacer_seq and len(c_spacer_seq)>1:
                        if crispr_direction=='+':
                            c_spacer_start = int(c_line.split()[0])+len(c_line.split()[1].replace(' ',''))
                        else:
                            c_spacer_seq = inverse_complement(c_spacer_seq)
                            c_spacer_start = int(c_line.split()[0])-len(c_spacer_seq)
                        spcfile.write('>' + str(crispr_id) + '.' +str(spacer_id)+ '|' +str(c_spacer_start)+'|'+str(len(c_spacer_seq))+'|'+contig_id.split('.')[0]+'\n')
                        spcfile.write(c_spacer_seq+'\n')
                        spacer_id = spacer_id+1
    cspfile.close()
    spcfile.close()

def run_CRISPRidentify(prog_dir,input_fasta_file,result_dir):
    mkdir(result_dir)
    # os.chdir("/zrom1/software/CRISPRidentify")
    # cmd_run_crispridentify = '/home/qiusuo/miniconda3/bin/python /zrom1/software/CRISPRidentify/CRISPRidentify.py --file %s --result_folder %s'%(input_fasta_file,result_dir)

    os.chdir("%s/software/CRISPRidentify"%(prog_dir))
    cmd_run_crispridentify = '/home/qiusuo/miniconda3/bin/python %s/software/CRISPRidentify/CRISPRidentify.py --file %s --result_folder %s'%(prog_dir,input_fasta_file,result_dir)

    print(cmd_run_crispridentify)
    os.system(cmd_run_crispridentify)

def combineCRISPRArray(pilercr_array_file,casfinder_array_file,crt_array_file,crispridentify_cspfile,merge_array_file):
    # create a dict storing the crispr array file
    crispr_info_dict = {}
    with open(pilercr_array_file,'r') as fin:
        lines = fin.readlines()
        for line in lines:
            content = line.strip('\n').split('\t')
            crispr_array_id = content[0]
            crispr_array_start = int(content[1])
            crispr_array_end = str(int(content[1])+int(content[2])-1)
            bacteria_id = content[5].strip().split()[0].split('.')[0]
            repeat_sequence = content[4].strip('-')
            if not bacteria_id in crispr_info_dict.keys():
                crispr_info_dict[bacteria_id] = []
                crispr_new = {}
                crispr_new['crispr_array_id'] = crispr_array_id
                crispr_new['crispr_array_start'] = crispr_array_start
                crispr_new['crispr_array_end'] = crispr_array_end
                crispr_new['method'] = 'PILER-CR'
                crispr_new['repeat'] = repeat_sequence
                crispr_info_dict[bacteria_id].append(crispr_new)
            else:
                crispr_new = {}
                crispr_new['crispr_array_id'] = crispr_array_id
                crispr_new['crispr_array_start'] = crispr_array_start
                crispr_new['crispr_array_end'] = crispr_array_end
                crispr_new['method'] = 'PILER-CR'
                crispr_new['repeat'] = repeat_sequence
                crispr_info_dict[bacteria_id].append(crispr_new)
    
    with open(casfinder_array_file,'r') as fin:
        lines = fin.readlines()
        for line in lines:
            content = line.strip('\n').split('\t')
            crispr_array_id = content[0]
            crispr_array_start = int(content[1])
            crispr_array_end = str(int(content[1])+int(content[2]))
            bacteria_id = content[5].strip().split()[0].split('.')[0]
            repeat_sequence = content[4]
            if not bacteria_id in crispr_info_dict.keys():
                crispr_info_dict[bacteria_id] = []
                crispr_new = {}
                crispr_new['crispr_array_id'] = crispr_array_id
                crispr_new['crispr_array_start'] = crispr_array_start
                crispr_new['crispr_array_end'] = crispr_array_end
                crispr_new['method'] = 'CRISPRCasFinder'
                crispr_new['repeat'] = repeat_sequence
                crispr_info_dict[bacteria_id].append(crispr_new)
            else:
                crispr_new = {}
                crispr_new['crispr_array_id'] = crispr_array_id
                crispr_new['crispr_array_start'] = crispr_array_start
                crispr_new['crispr_array_end'] = crispr_array_end
                crispr_new['method'] = 'CRISPRCasFinder'
                crispr_new['repeat'] = repeat_sequence
                crispr_info_dict[bacteria_id].append(crispr_new)
    
    with open(crt_array_file,'r') as fin:
        lines = fin.readlines()
    if len(lines)>0:
        for line in lines:
            content = line.strip('\n').split('\t')
            crispr_array_id = content[0]
            crispr_array_start = int(content[1])
            crispr_array_end = str(int(content[1])+int(content[2]))
            bacteria_id = content[5].strip().split()[0].split('.')[0]
            repeat_sequence = content[-2]
            if not bacteria_id in crispr_info_dict.keys():
                crispr_info_dict[bacteria_id] = []
                crispr_new = {}
                crispr_new['crispr_array_id'] = crispr_array_id
                crispr_new['crispr_array_start'] = crispr_array_start
                crispr_new['crispr_array_end'] = crispr_array_end
                crispr_new['method'] = 'CRT'
                crispr_new['repeat'] = repeat_sequence
                crispr_info_dict[bacteria_id].append(crispr_new)
            else:
                crispr_new = {}
                crispr_new['crispr_array_id'] = crispr_array_id
                crispr_new['crispr_array_start'] = crispr_array_start
                crispr_new['crispr_array_end'] = crispr_array_end
                crispr_new['method'] = 'CRT'
                crispr_new['repeat'] = repeat_sequence
                crispr_info_dict[bacteria_id].append(crispr_new)

    #crispridentify
    with open(crispridentify_cspfile,'r') as fin:
        lines = fin.readlines()
    if len(lines)>0:
        for line in lines:
            content = line.strip('\n').split('\t')
            crispr_array_id = content[0]
            crispr_array_start = int(content[1])
            crispr_array_end = str(int(content[1])+int(content[2]))
            bacteria_id = content[5].strip().split()[0].split('.')[0]
            repeat_sequence = content[-2]
            if not bacteria_id in crispr_info_dict.keys():
                crispr_info_dict[bacteria_id] = []
                crispr_new = {}
                crispr_new['crispr_array_id'] = crispr_array_id
                crispr_new['crispr_array_start'] = crispr_array_start
                crispr_new['crispr_array_end'] = crispr_array_end
                crispr_new['method'] = 'CRISPRidentify'
                crispr_new['repeat'] = repeat_sequence
                crispr_info_dict[bacteria_id].append(crispr_new)
            else:
                crispr_new = {}
                crispr_new['crispr_array_id'] = crispr_array_id
                crispr_new['crispr_array_start'] = crispr_array_start
                crispr_new['crispr_array_end'] = crispr_array_end
                crispr_new['method'] = 'CRISPRidentify'
                crispr_new['repeat'] = repeat_sequence
                crispr_info_dict[bacteria_id].append(crispr_new)
    

    for key in crispr_info_dict:
        crispr_info_dict[key].sort(key=lambda x: int(x['crispr_array_start']))

    # merge region
    # combinedList = []
    for key in crispr_info_dict:
        templist = crispr_info_dict[key][0]
        combinedList = []
        for newlist in crispr_info_dict[key][1:]:
            newListmin = min(int(newlist['crispr_array_start']),int(newlist['crispr_array_end']))
            newListmax = max(int(newlist['crispr_array_start']),int(newlist['crispr_array_end']))
            tempListmax = max(int(templist['crispr_array_start']),int(templist['crispr_array_end']))
            if int(newListmin) <= int(tempListmax):
                templist['crispr_array_end'] = str(max(int(tempListmax), int(newListmax)))
                templist['crispr_array_id'] += ',' + newlist['crispr_array_id']
                templist['method'] += ',' + newlist['method']
                templist['repeat'] += ',' + newlist['repeat']
            else:
                combinedList.append(templist)
                templist = newlist
        combinedList.append(templist)
        crispr_info_dict[key]=(combinedList)

    # save merged result
    fout = open(merge_array_file,'w')
    for key in crispr_info_dict:
        num_array = 0
        for item in crispr_info_dict[key]:
            num_array = num_array + 1
            bacteria_id = key
            crispr_array_id = str(item['crispr_array_id'])
            crispr_array_start = str(item['crispr_array_start'])
            crispr_array_end = str(item['crispr_array_end'])
            crispr_array_method = item['method']
            crispr_array_repeat = item['repeat']
            strWrite = str(num_array)+'\t'+bacteria_id+'\t'+crispr_array_id+'\t'+crispr_array_start+'\t'+crispr_array_end+'\t'+crispr_array_repeat+'\t'+crispr_array_method+'\n'
            fout.write(strWrite)
            fout.flush()
    fout.close()

def combineSpacer(pilercr_spc_file,casfinder_spc_file,crt_spc_file,crispridentify_spcfile,merge_array_file,merge_spc_file):
    dictSpacer = {}
    with open(pilercr_spc_file,'r')as fin:
        content = fin.read()
        elems = content.split('>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            bacteria_id = elem.split('\n')[0].split('|')[-1].split('.')[0]
            crispr_array_id = elem.split('\n')[0].split('|')[0].split('.')[0]
            spacer_start = elem.split('\n')[0].split('|')[1]
            spacer_defLine = elem.split('\n')[0]
            spacer_seq = elem.split('\n')[1]
            spacer_len = elem.split('\n')[0].split('|')[2]
            method = 'PILER-CR'
            if bacteria_id not in dictSpacer.keys():
                dictSpacer.update({bacteria_id:{}})
            if method not in dictSpacer[bacteria_id].keys():
                dictSpacer[bacteria_id].update({method:[]})
            spacer_info_new = {}
            spacer_info_new['crispr_array_id'] = crispr_array_id
            spacer_info_new['spacer_start'] = spacer_start
            spacer_info_new['spacer_len'] = spacer_len
            spacer_info_new['spacer_defLine'] = spacer_defLine
            spacer_info_new['spacer_seq'] = spacer_seq
            spacer_info_new['method'] = method
            dictSpacer[bacteria_id][method].append(spacer_info_new)
    
    with open(casfinder_spc_file,'r')as fin:
        content = fin.read()
        elems = content.split('>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            bacteria_id = elem.split('\n')[0].split('|')[-1].split('.')[0].strip()
            crispr_array_id = elem.split('\n')[0].split('|')[0].split('.')[0]
            spacer_start = elem.split('\n')[0].split('|')[1]
            spacer_defLine = elem.split('\n')[0]
            spacer_len = elem.split('\n')[0].split('|')[2]
            spacer_seq = elem.split('\n')[1]
            method = 'CRISPRCasFinder'
            if bacteria_id not in dictSpacer.keys():
                dictSpacer.update({bacteria_id:{}})
            if method not in dictSpacer[bacteria_id].keys():
                dictSpacer[bacteria_id].update({method:[]})
            spacer_info_new = {}
            spacer_info_new['crispr_array_id'] = crispr_array_id
            spacer_info_new['spacer_start'] = spacer_start
            spacer_info_new['spacer_len'] = spacer_len
            spacer_info_new['spacer_defLine'] = spacer_defLine
            spacer_info_new['spacer_seq'] = spacer_seq
            spacer_info_new['method'] = method
            dictSpacer[bacteria_id][method].append(spacer_info_new)
    
    with open(crt_spc_file,'r')as fin:
        content = fin.read()
        elems = content.split('>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            bacteria_id = elem.split('\n')[0].split('|')[-1].split('.')[0]
            crispr_array_id = elem.split('\n')[0].split('|')[0].split('.')[0]
            spacer_start = elem.split('\n')[0].split('|')[1]
            spacer_defLine = elem.split('\n')[0]
            spacer_len = elem.split('\n')[0].split('|')[2]
            spacer_seq = elem.split('\n')[1]
            method = 'CRT'
            if bacteria_id not in dictSpacer.keys():
                dictSpacer.update({bacteria_id:{}})
            if method not in dictSpacer[bacteria_id].keys():
                dictSpacer[bacteria_id].update({method:[]})
            spacer_info_new = {}
            spacer_info_new['crispr_array_id'] = crispr_array_id
            spacer_info_new['spacer_start'] = spacer_start
            spacer_info_new['spacer_len'] = spacer_len
            spacer_info_new['spacer_defLine'] = spacer_defLine
            spacer_info_new['spacer_seq'] = spacer_seq
            spacer_info_new['method'] = method
            dictSpacer[bacteria_id][method].append(spacer_info_new)
    
    with open(crispridentify_spcfile,'r')as fin:
        content = fin.read()
        elems = content.split('>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            bacteria_id = elem.split('\n')[0].split('|')[-1].split('.')[0].strip()
            crispr_array_id = elem.split('\n')[0].split('|')[0].split('.')[0]
            spacer_start = elem.split('\n')[0].split('|')[1]
            spacer_defLine = elem.split('\n')[0]
            spacer_len = elem.split('\n')[0].split('|')[2]
            spacer_seq = elem.split('\n')[1]
            method = 'CRISPRidentify'
            if bacteria_id not in dictSpacer.keys():
                dictSpacer.update({bacteria_id:{}})
            if method not in dictSpacer[bacteria_id].keys():
                dictSpacer[bacteria_id].update({method:[]})
            spacer_info_new = {}
            spacer_info_new['crispr_array_id'] = crispr_array_id
            spacer_info_new['spacer_start'] = spacer_start
            spacer_info_new['spacer_len'] = spacer_len
            spacer_info_new['spacer_defLine'] = spacer_defLine
            spacer_info_new['spacer_seq'] = spacer_seq
            spacer_info_new['method'] = method
            dictSpacer[bacteria_id][method].append(spacer_info_new)

    spacer_write_dict = OrderedDict()
    
    with open(merge_array_file,'r')as fin:
        lines = fin.readlines()
        for line in lines:
            content = line.strip('\n').split('\t')
            array_id = content[0]
            bacteria_id = content[1].split('.')[0]
            crispr_array_original = content[2].split(',')
            method_original = content[-1].split(',')
            if bacteria_id not in spacer_write_dict.keys():
                spacer_write_dict.update({bacteria_id:{}})
            if array_id not in spacer_write_dict[bacteria_id].keys():
                count_spacer = 0
                spacer_write_dict[bacteria_id].update({array_id:{}})
            for i in range(0,len(method_original)):
                crispr_array_original_info = crispr_array_original[i]
                method_origianl_info = method_original[i]
                for spacer_info in dictSpacer[bacteria_id][method_origianl_info]:
                    if spacer_info['crispr_array_id'] == crispr_array_original_info:
                        spacer_seq = spacer_info['spacer_seq']
                        spacer_start = spacer_info['spacer_start']
                        if spacer_start+'_'+spacer_seq not in spacer_write_dict[bacteria_id][array_id].keys():
                            count_spacer = count_spacer+1
                            spacer_write_dict[bacteria_id][array_id].update({spacer_start+'_'+spacer_seq:{'locus_id':str(array_id) + '.' + str(count_spacer),
                                'spacer_start':spacer_info['spacer_start'],
                                'spacer_len':spacer_info['spacer_len'],
                                'spacer_seq':spacer_info['spacer_seq'],
                                'method':[]}})
                        spacer_write_dict[bacteria_id][array_id][spacer_start+'_'+spacer_seq]['method'].append(method_origianl_info)    
    
    fout = open(merge_spc_file,'w')
    for bac_key in spacer_write_dict.keys():
        num_locus = len(spacer_write_dict[bac_key])
        for array_id in spacer_write_dict[bac_key].keys():
            locus_key = str(i)
            for spacer_id,spacer_item in spacer_write_dict[bac_key][array_id].items():
                spacer_defLine = '>'+str(spacer_item['locus_id'])+'|'+str(spacer_item['spacer_start'])+'|'+str(spacer_item['spacer_len'])+'|'+bac_key+'|'+','.join(spacer_item['method'])+'\n'
                spacer_seqLine = spacer_item['spacer_seq']+'\n'
                fout.write(spacer_defLine+spacer_seqLine)
                fout.flush()
    fout.close()

def combineFile(pilercr_array_file, casfinder_array_file, crt_array_file,crispridentify_cspfile, merge_array_file, pilercr_spc_file,
                casfinder_spc_file, crt_spc_file,crispridentify_spcfile, merge_spc_file):
    if not os.path.exists(pilercr_array_file):
        cmd_touch = 'touch %s' % (pilercr_array_file)
        os.system(cmd_touch)
        cmd_touch = 'touch %s' % (pilercr_spc_file)
        os.system(cmd_touch)

    if not os.path.exists(casfinder_array_file):
        cmd_touch = 'touch %s' % (casfinder_array_file)
        os.system(cmd_touch)
        cmd_touch = 'touch %s' % (casfinder_spc_file)
        os.system(cmd_touch)

    if not os.path.exists(crt_array_file):
        cmd_touch = 'touch %s' % (crt_array_file)
        os.system(cmd_touch)
        cmd_touch = 'touch %s' % (crt_spc_file)
        os.system(cmd_touch)

    if not os.path.exists(crispridentify_cspfile):
        cmd_touch = 'touch %s' % (crispridentify_cspfile)
        os.system(cmd_touch)
        cmd_touch = 'touch %s' % (crispridentify_spcfile)
        os.system(cmd_touch)

    #if os.path.exists(pilercr_array_file) and os.path.exists(casfinder_array_file) and os.path.exists(crt_array_file):
    combineCRISPRArray(pilercr_array_file, casfinder_array_file, crt_array_file,crispridentify_cspfile, merge_array_file)
    combineSpacer(pilercr_spc_file, casfinder_spc_file, crt_spc_file,crispridentify_spcfile, merge_array_file, merge_spc_file)
        

def combine_repeat(merge_cspfile,outfile):
    f_result = open(outfile,'w')
    with open(merge_cspfile) as f:
        contents = f.readlines()
    for line in contents:
        line = line.strip().split('\t')
        crispr_id = line[0]
        contig_id = line[1]
        repeat = line[5].split(',')
        method = line[-1].split(',')
        repeats = list(zip(repeat,method))
        for item in repeats:
            f_result.write('>'+contig_id+'_'+crispr_id+'_'+item[1]+'\n'+item[0].strip('-')+'\n')
    f_result.close()
    return outfile

def get_crispr_sequence(outdir,fa_file,merge_cspfile,merge_spcfile,pilercr_array_file,casfinder_array_file,crt_array_file,crispridentify_cspfile,save_file):
    contig_sequence_dict = {}
    with open(fa_file) as f:
        contents = f.read().split('>')
    for contig in contents[1:]:
        contig_id = contig.split('\n')[0].split()[0].split('.')[0]
        contig_sequence_dict.update({contig_id:''.join(contig.split('\n')[1:]).strip()})

    crispr_array_location_dict  ={}
    file_list = [['PILER-CR',pilercr_array_file],
    ['CRISPRCasFinder',casfinder_array_file],
    ['CRT',crt_array_file],
    ['CRISPRidentify',crispridentify_cspfile]]
    
    for crispr in file_list:
        crispr_file = crispr[1]
        method = crispr[0]
        if os.path.exists(crispr_file):
            with open(crispr_file) as f:
                contents = f.readlines()
            for line in contents:
                line = line.strip().split('\t')
                crispr_start = line[1]
                crispr_shift = line[2]
                if method=='PILER-CR':
                    crispr_end = int(crispr_start)+int(crispr_shift)-1
                else:
                    crispr_end = int(crispr_start)+int(crispr_shift)
                contig_id = line[5].split('.')[0]
                crispr_id = line[0]
                if contig_id not in crispr_array_location_dict.keys():
                    crispr_array_location_dict.update({contig_id:{}})
                if method not in crispr_array_location_dict[contig_id].keys():
                    crispr_array_location_dict[contig_id].update({method:{}})                   
                crispr_array_location_dict[contig_id][method].update({crispr_id:[int(crispr_start),int(crispr_end)]})
    
    merge_spc_dict = {}
    with open(merge_spcfile) as f:
        contents = f.read()
    for spacer in contents.split('>')[1:]:
        c_spacer_id = spacer.split('\n')[0].strip()
        c_spacer_start = c_spacer_id.split('|')[1]
        c_spacer_length = c_spacer_id.split('|')[2]
        contig_id = c_spacer_id.split('|')[3]
        c_spacer_method = c_spacer_id.split('|')[-1].split(',')
        c_crispr_id = c_spacer_id.split('|')[0].split('.')[0]
        c_spacer_end = int(c_spacer_start)+int(c_spacer_length)-1
        if contig_id not in merge_spc_dict.keys():
            merge_spc_dict.update({contig_id:{}})
        for index,method in enumerate(c_spacer_method):
            if method not in merge_spc_dict[contig_id].keys():
                merge_spc_dict[contig_id].update({method:{}})
            if c_crispr_id not in merge_spc_dict[contig_id][method].keys():
                merge_spc_dict[contig_id][method].update({c_crispr_id:[]})
            merge_spc_dict[contig_id][method][c_crispr_id].append([int(c_spacer_start),c_spacer_end])
    
    f_save = open(save_file,'w')
    with open(merge_cspfile) as f:
        contents = f.readlines()
    for line in contents:
        line = line.strip().split('\t')
        contig_id = line[1]
        c_contig_sequence = contig_sequence_dict[contig_id].strip()
        crispr_id = line[0]
        locus_id_list = line[2].split(',')
        method_list = line[-1].split(',')
        merge_crispr_start = line[3]
        merge_crispr_end = line[4]
        f_save.write('Array: '+contig_id+'|'+crispr_id+'|'+merge_crispr_start+'-'+merge_crispr_end+'\n')
        f_save.write('>merge|'+contig_id+'|'+crispr_id+'|'+merge_crispr_start+'-'+merge_crispr_end+'|'+','.join(method_list)+'\n')
        f_save.write(c_contig_sequence[int(merge_crispr_start)-1:int(merge_crispr_end)]+'\n\n')
        for index,method in enumerate(method_list):
            locus_id = locus_id_list[index]
            c_crispr_location = crispr_array_location_dict[contig_id][method][locus_id]
            c_spacer_locations = merge_spc_dict[contig_id][method][crispr_id]
            c_spacer_locations = sorted(c_spacer_locations,key=lambda x:x[0])
            crispr_start = c_crispr_location[0]
            crispr_end = c_crispr_location[1]
            c_repeat_start = crispr_start
            c_repeat_end = c_spacer_locations[0][0]-1
            f_save.write('>'+contig_id+'|'+crispr_id+'|'+locus_id+'|'+str(crispr_start)+'-'+str(crispr_end)+'|'+method+'\n')
            for index,c_spacer_location in enumerate(c_spacer_locations):               
                f_save.write(c_contig_sequence[c_repeat_start-1:c_repeat_end]+'\t'+c_contig_sequence[c_spacer_location[0]-1:c_spacer_location[1]]+'\n')
                c_repeat_start = c_spacer_location[1]+1
                if index!=len(c_spacer_locations)-1:
                    c_repeat_end = c_spacer_locations[index+1][0]-1
            f_save.write(c_contig_sequence[c_repeat_start-1:crispr_end]+'\n')
        f_save.write('\n\n\n')
    f_save.close()

def identify_crispr_array(prog_dir,work_dir,bac_dir,input_file,strain_inf_file):
    contig_dir = '%s/%s/protein' % (work_dir,bac_dir)
    crt_dir = '%s/%s/crt_dir' % (work_dir, bac_dir)# save the crt results
    mkdir(crt_dir)
    strain_file = input_file
    ##1.1.predict CRISPR array by CRT:
    crt_cspfile = os.path.join(crt_dir,'crt.csp')
    crt_spcfile = os.path.join(crt_dir,'crt.spc')
    ##1.2 predict CRISPR array by CRISRPCasFinder:
    casfinder_dir = '%s/%s/casfinder'%(work_dir,bac_dir)
    casfinder_result_file = '%s/%s/casfinder/result.json'%(work_dir,bac_dir)
    casfinder_spcfile = '%s/%s/%s_casfinder.spc'%(work_dir,bac_dir,bac_dir)
    casfinder_cspfile = '%s/%s/%s_casfinder.csp'%(work_dir,bac_dir,bac_dir)
    ##1.3.predict CRISPR array by PILER-CR:
    pilercr_result_file = '%s/%s/%s_pliercr.txt'%(work_dir,bac_dir,bac_dir)
    pilercr_spcfile = '%s/%s/%s_pliercr.spc'%(work_dir,bac_dir,bac_dir)
    pilercr_cspfile = '%s/%s/%s_pliercr.csp'%(work_dir,bac_dir,bac_dir)
    ##1.4 predict crispr array by CRISPRidentify
    CRISPRidentify_dir = '%s/%s/crispridentify'%(work_dir,bac_dir)
    mkdir(CRISPRidentify_dir)
    crispridentify_spcfile = '%s/%s/%s_crispridentify.spc'%(work_dir,bac_dir,bac_dir)
    crispridentify_cspfile = '%s/%s/%s_crispridentify.csp'%(work_dir,bac_dir,bac_dir)

    muthread1 = threading.Thread(target=runPILERCR, args=(prog_dir,strain_file, pilercr_result_file,))
    muthread2 = threading.Thread(target=multifasta_crt,args=(prog_dir,strain_inf_file,contig_dir, crt_dir,))
    muthread3 = threading.Thread(target=runCRISPRFinder,args=(prog_dir,strain_file, casfinder_dir,))
    muthread4 = threading.Thread(target=run_CRISPRidentify,args=(prog_dir,strain_file, CRISPRidentify_dir,))
    muthread1.start()
    muthread2.start()
    muthread3.start()
    muthread4.start()
    muthread1.join()
    muthread2.join()
    muthread3.join()
    muthread4.join()
    
    ##2.1.parse casfinder result:
    try:
        parse_casfinder_result(casfinder_result_file, casfinder_spcfile, casfinder_cspfile)
    except:
        pass
    
    ##2.2.parse pilercr result:
    parse_pildercr_result(pilercr_result_file, pilercr_spcfile, pilercr_cspfile)

    ##2.3.parse crispridentify result:
    try:
        parse_crispridentify_result(CRISPRidentify_dir, crispridentify_spcfile, crispridentify_cspfile)
    except:
        pass
    
    ##3.merge the repeat and spacer detected by casfinder,pilercr,crt
    merge_cspfile = '%s/%s/%s_merge.csp'%(work_dir,bac_dir,bac_dir)
    merge_spcfile = '%s/%s/%s_merge.spc'%(work_dir,bac_dir,bac_dir)
    combineFile(pilercr_cspfile,casfinder_cspfile,crt_cspfile,crispridentify_cspfile,merge_cspfile,pilercr_spcfile,casfinder_spcfile,crt_spcfile,crispridentify_spcfile,merge_spcfile)

    ##4. identify repeat type
    repeat_file = '%s/%s/%s_repeat.fa'%(work_dir,bac_dir,bac_dir)
    combine_repeat(merge_cspfile,repeat_file)
    out_repeatfile = '%s/%s/%s_repeat_type'%(work_dir,bac_dir,bac_dir)
    identify_repeat_type(prog_dir,repeat_file,out_repeatfile)

    ##5 write CRISPR-spacer sequence into files
    crispr_array_file = '%s/%s/merge_crispr_array.txt'%(work_dir,bac_dir)
    get_crispr_sequence(strain_inf_file,strain_file,merge_cspfile,merge_spcfile,pilercr_cspfile,casfinder_cspfile,crt_cspfile,crispridentify_cspfile,crispr_array_file)

def pred_orf(prog_dir,fasta_file,faa_prefix):
    outfile_pro = faa_prefix+'_prodigal.faa'
    gff_file = faa_prefix+'_prodigal.gff'
    command = "%s/software/prodigal -i %s -o %s -a %s -p meta"%(prog_dir,fasta_file,gff_file,outfile_pro)
    os.system(command)
    if os.path.exists(outfile_pro):
        save_pro_file = faa_prefix+'.faa'
        f_save = open(save_pro_file,'w')
        with open(outfile_pro) as f:
            pro_contents = f.read()
        for pro in pro_contents.split('\n>'):
            title = pro.split('\n')[0].strip().strip('>')
            bac_id = '_'.join(title.split('#')[0].strip().split('_')[0:-1])
            prot_start = title.split('#')[1].strip()
            prot_end = title.split('#')[2].strip()
            prot_direction = title.split('#')[3].strip()
            if prot_direction=='1':
                prot_direction = '+'
            else:
                prot_direction = '-'
            f_save.write('>'+bac_id+'_'+prot_start+'_'+prot_end+'_'+prot_direction+'\n'+''.join(pro.split('\n')[1:]).strip()+'\n')
            f_save.flush()
    f_save.close()
def getprotinfov2(prog_dir,xml_file,faa_file,cpt_file):
    cmd = '%s/software/hmmscan -E 1e-9 --tblout %s %s/database/Cas_protein/cas_transposon.hmm %s > /dev/null' % (prog_dir, xml_file,prog_dir, faa_file)
    os.system(cmd)
    f2 = open(cpt_file, 'w')
    if not os.path.exists(xml_file):
        return 0
    else:
        if os.path.getsize(xml_file)==0:
           return 0
    f1 = open(xml_file, 'r')
    genequery = {}
    for line in f1:
        if line[0] == '#': continue
        linetab = line.split()
        if linetab[2] not in genequery.keys():
            def_line = linetab[2]
            if '<' in def_line:
                def_line = def_line.replace('<', '')
            if '>' in def_line:
                def_line = def_line.replace('>', '')

            if 'join' in def_line:
                def_line = def_line.replace('join{', '')
                def_line = def_line.replace('}', '')
            genequery[linetab[2]] = []
            genequery[linetab[2]].append(linetab[0])  # 0 cas protein name
            genequery[linetab[2]].append(def_line)  # 1 query name
            genequery[linetab[2]].append(linetab[4])  # 2 e-value
            genequery[linetab[2]].append(linetab[5])  # 3 score
            genequery[linetab[2]].append(linetab[6])  # 4 bias
            genequery[linetab[2]].append(linetab[18])  # 5 descriprion
            continue
        src_evalue = float(genequery[linetab[2]][2])
        src_score = float(genequery[linetab[2]][3])
        src_bias = float(genequery[linetab[2]][4])
        dst_evalue = float(linetab[4])
        dst_score = float(linetab[5])
        dst_bias = float(linetab[6])
        if (dst_evalue < src_evalue):
            genequery[linetab[2]][0] = linetab[0]
            genequery[linetab[2]][1] = linetab[2]
            genequery[linetab[2]][2] = linetab[4]
            genequery[linetab[2]][3] = linetab[5]
            genequery[linetab[2]][4] = linetab[6]
            genequery[linetab[2]][5] = linetab[18]
        elif (dst_score > src_score and dst_evalue == src_evalue):
            genequery[linetab[2]][0] = linetab[0]
            genequery[linetab[2]][1] = linetab[2]
            genequery[linetab[2]][2] = linetab[4]
            genequery[linetab[2]][3] = linetab[5]
            genequery[linetab[2]][4] = linetab[6]
            genequery[linetab[2]][5] = linetab[18]
        elif (dst_bias < src_bias and dst_score == src_score and dst_evalue == src_evalue):
            genequery[linetab[2]][0] = linetab[0]
            genequery[linetab[2]][1] = linetab[2]
            genequery[linetab[2]][2] = linetab[4]
            genequery[linetab[2]][3] = linetab[5]
            genequery[linetab[2]][4] = linetab[6]
            genequery[linetab[2]][5] = linetab[18]
    
    # genequery = sorted(genequery.items(), key=lambda e: float(e[0].split('_')[-4]))
    for item in genequery:
        f2.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (genequery[item][1], genequery[item][0], genequery[item][2], genequery[item][3], genequery[item][4], genequery[item][5]))
    f1.close()
    f2.close()

def rpsblastp11proteins(prog_dir,fileName,outFileName):
    blast_cline = "%s/software/rpsblast -query %s -comp_based_stats 0 -max_target_seqs 1 -evalue 0.01 -seg no -outfmt 5 -num_threads 20 -db %s/database/CDD_rps/Cdd -out %s" % (prog_dir,fileName,prog_dir,outFileName)
    os.system(blast_cline)

def ParseResult(inFileName,outFileName):
    f = open(outFileName,'w')
    if not os.path.exists(inFileName):
        return 0
    else:
        if os.path.getsize(inFileName)==0:
            return 0
    text = open(inFileName).read()
    text = re.sub(u"[\x00-\x08\x0b-\x0c\x0e-\x1f]+", u" ", text)
    #print(text)
    root = ET.fromstring(text)
    BlastOutput_iterations = root.find("BlastOutput_iterations")
    
    for Iteration in BlastOutput_iterations.findall("Iteration"):
        strTemp = str(Iteration.find("Iteration_query-def").text)
        Iteration_hits = Iteration.find("Iteration_hits")
        for Hit in Iteration_hits.findall("Hit"):
            strDef = str(Hit.find("Hit_def").text)
            Hit_len =  Hit.find("Hit_len").text
            Hit_ID = Hit.find("Hit_id").text
            Hsp = Hit.find("Hit_hsps").find("Hsp")
            Hit_from = Hsp.find("Hsp_hit-from").text
            Hit_to = Hsp.find("Hsp_hit-to").text
            Hsp_evalue = Hsp.find("Hsp_evalue").text
            Hsp_score = Hsp.find("Hsp_score").text
            hit_length = int(Hit_from)-int(Hit_from)+1
            coverage = str(float(hit_length)/ float(Hit_len))
            identity = str(Hsp.find("Hsp_identity").text)
            f.write(strTemp+'\t'+Hit_ID+'\t'+strDef+'\t'+Hsp_evalue+'\t'+identity+'\t'+coverage+'\n')
            f.flush()
    f.close()

def annotation_cas_prot(prog_dir,work_dir,bac_dir,protein_file):
    xml_file = '%s/%s/%s.xml' % (work_dir, bac_dir, bac_dir)
    cpt_file = '%s/%s/%s.cpt' % (work_dir, bac_dir, bac_dir)

    getprotinfov2(prog_dir,xml_file, protein_file, cpt_file)

def blastn(prog_dir,query_file,subject_file,result_file):
    cmd_blastn = '%s/software/blastn -query %s -subject %s -out %s -evalue 1 -word_size 7 -penalty -1 -reward 1 -outfmt 6 -dust no -soft_masking FALSE  -ungapped' % (prog_dir,query_file, subject_file, result_file)
    os.system(cmd_blastn)

def get_def_dict(fasta_file_path,genome_def_dict):
    with open(fasta_file_path,'r')as fin:
        content = fin.read()
        elems = content.split('>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            defLine = elem.split('\n')[0]
            cur_genome_id = defLine.split(' ')[0].split('.')[0]
            cur_genome_def = ' '.join(defLine.split(' ')[1:])
            if cur_genome_id not in genome_def_dict.keys():
                genome_def_dict[cur_genome_id] = cur_genome_def

def add_cov(spc_blastn_fasta_file_path,blastn_add_cov_file,genome_def_dict):
    fout = open(blastn_add_cov_file, 'w')
    query_assembly_id = os.path.basename(spc_blastn_fasta_file_path).rstrip('.merge_spc_blastn_fasta')
    header = 'query_assembly_id\tquery_id\tquery_genome_def\tquery_spacer_info\thit_id\thit_genome_def\tidentity\talignment_length\tmismatch\tgapopen\tquery_start\tquery_end\thit_start\thit_end\tevalue\tbitscore\tcoverage\n'
    fout.write(header)
    with open(spc_blastn_fasta_file_path, 'r')as fin:
        spacerLines = fin.readlines()
        for spacerLine in spacerLines:
            spacerElems = spacerLine.strip('\n').split('\t')
            query_genome_id = spacerElems[0].split('|')[3]
            hit_genome_id = spacerElems[1].split('.')[0]
            query_genome_def = genome_def_dict[query_genome_id]
            hit_genome_def = genome_def_dict[hit_genome_id]
            query_genome_length = float(spacerElems[0].split('|')[2])
            align_length = float(spacerElems[3])
            coverage_value = (align_length) / query_genome_length
            strWrite = '\t'.join(spacerElems[2:])
            strWrite = query_assembly_id + '\t' + query_genome_id + '\t' + query_genome_def + '\t'+spacerElems[0]+'\t' + spacerElems[1] + '\t' + hit_genome_def + '\t' + strWrite+'\t'+ str(coverage_value)+ '\n'
            fout.write(strWrite)
            fout.flush
    fout.close()

def filter_mismatch_and_cov(blastn_add_cov_file, filter_mismatch_cov_file, mismatch='2', coverage='1'):
    fout = open(filter_mismatch_cov_file, 'w')
    with open(blastn_add_cov_file, 'r')as fin:
        lines = fin.readlines()
        fout.write(lines[0])
        for line in lines[1:]:
            cur_mismatch = line.strip('\n').split('\t')[8]
            cur_cov = line.strip('\n').split('\t')[-1]
            if int(cur_mismatch) <= int(mismatch) and float(cur_cov) >= float(coverage):
                fout.write(line)
    fout.close()

def get_array_info_dict(merge_csp_file,array_info_dict):
    with open(merge_csp_file,'r')as fin:
        lines = fin.readlines()
        for line in lines:
            content = line.strip('\n').split('\t')
            cur_genome_id = content[1]
            cur_array_start = min(int(content[3]),int(content[4]))
            cur_array_end = max(int(content[3]),int(content[4]))
            if cur_genome_id not in array_info_dict.keys():
                array_info_dict[cur_genome_id] = []
            temp_list = []
            temp_list.append(str(cur_array_start))
            temp_list.append(str(cur_array_end))
            array_info_dict[cur_genome_id].append(temp_list)

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

def filter_in_array(prog_dir,filter_mismatch_cov_file,filter_hit_in_array,array_info_dict,st_dir,spacer_file,contig_dir):
    fout = open(filter_hit_in_array,'w')
    with open(filter_mismatch_cov_file,'r')as fin:
        lines = fin.readlines()
        header_list = lines[0].strip('\n').split('\t')
        add_head = ['crispr_array_of_hit_genome']
        fout.write('\t'.join(header_list[0:14]+add_head+header_list[14:])+'\n')
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_hit_id = content[4].split('.')[0]
            cur_hit_from = min(int(content[-4]),int(content[-5]))
            cur_hit_to = max(int(content[-4]),int(content[-5]))
            hit_in_array_flag = 0
            if cur_hit_id in array_info_dict:
                array_info_list = array_info_dict[cur_hit_id]
                for array_info_item in array_info_list:
                    if cur_hit_from>=int(array_info_item[0]) and cur_hit_to<=int(array_info_item[1]):
                        hit_in_array_flag = 1
                        break
                if hit_in_array_flag == 0:
                    str_array_info = []
                    for array_info_item in array_info_list:
                        str_array_info.append(array_info_item[0]+'-'+array_info_item[1])
                    content_list = line.strip('\n').split('\t')
                    fout.write('\t'.join(content_list[0:14])+'\t'+'|'.join(str_array_info)+'\t'+'\t'.join(content_list[14:]) + '\n')
            else:
                content_list = line.strip('\n').split('\t')
                str_array_info = ['no crispr array in hit genome!']
                fout.write('\t'.join(content_list[0:14]+str_array_info+content_list[14:])+'\n')
    fout.close()

    #multi-align file
    spacer_sequence_dict = {}
    with open(spacer_file) as f:
        spacer_contents = f.read().split('>')
    for spacer in spacer_contents[1:]:
        spacer_id = spacer.split('\n')[0]
        spacer_sequence_dict.update({spacer_id:spacer})
    
    spacer_target_file = os.path.join(st_dir,'spacer_target_hit_result.txt')
    f_save = open(spacer_target_file,'w')
    f_save.write('bac_id\tbac_def\tspacer_id\thit_contig_id\thit_contig_def\tidentity\talignment_length\tmismatch\tgapopen\tquery_start\tquery_end\thit_start\thit_end\tcrispr_array_of_hit_genome\tevalue\tbitscore\tcoverage\tspacer_sequence\thit_contig_region\thit_contig_sequence\thmm_mismatch\n') 
    f_save.flush()

    with open(filter_hit_in_array) as f:
        contents = f.readlines()
    header = contents[0].strip().split('\t')
    for line in contents[1:]:
        line = line.strip().split('\t')
        spacer_id = line[header.index('query_spacer_info')]
        hit_id = line[header.index('hit_id')]
        hit_contig_file = os.path.join(contig_dir,hit_id+'.fa')
        with open(hit_contig_file) as f:
            hit_contig_sequence = f.read()
        hit_contig_sequence = ''.join(hit_contig_sequence.split('\n')[1:]).strip().replace(' ','')
        spacer_start = line[header.index('query_start')]
        spacer_end = line[header.index('query_end')]
        hit_start = line[header.index('hit_start')]
        hit_end = line[header.index('hit_end')]
        spacer_start = min(int(line[header.index('query_start')]),int(line[header.index('query_end')]))
        spacer_end = max(int(line[header.index('query_start')]),int(line[header.index('query_end')]))
        hit_start = min(int(line[header.index('hit_start')]),int(line[header.index('hit_end')]))
        hit_end =  max(int(line[header.index('hit_start')]),int(line[header.index('hit_end')]))
        
        if int(line[header.index('query_start')])<int(line[header.index('query_end')]):
            query_direction = '+'
        else:
            query_direction = '-'       
        if int(line[header.index('hit_start')])<int(line[header.index('hit_end')]):
            hit_direction = '+'
        else:
            hit_direction = '-'               
        c_spacer_length = int(spacer_id.split('|')[2])       
        if query_direction==hit_direction:
            direction = 'same'
            cur_subject_from = max(0,int(hit_start)-(int(spacer_start)-1)-1)
            cur_subject_to = min(len(hit_contig_sequence),(int(hit_end)+(c_spacer_length-int(spacer_end))))
        else:
            direction = 'different' 
            cur_subject_from = max(0,int(hit_start)-(c_spacer_length-int(spacer_end))-1)
            cur_subject_to = min(len(hit_contig_sequence),(int(hit_end)+(int(spacer_start)-1)))
        
        expand_sequence = hit_contig_sequence[cur_subject_from:cur_subject_to]                 
        c_spacer_sequence = spacer_sequence_dict[spacer_id]
        
        spacer_length = len(c_spacer_sequence.split('\n')[1].strip())
        c_align_st_file = os.path.join(st_dir,spacer_id.replace('|','-')+'_'+hit_id+'_'+str(cur_subject_from+1)+'-'+str(cur_subject_to)+'.fa')
        with open(c_align_st_file,'w') as f:
            f.write('>'+c_spacer_sequence.strip('>').strip()+'\n'+'>'+hit_id+'|'+str(cur_subject_from+1)+'-'+str(cur_subject_to)+'\n'+expand_sequence)
       
        c_align_st_mutalign_file = os.path.join(st_dir,spacer_id.replace('|','-')+'_'+hit_id+'_'+str(cur_subject_from+1)+'-'+str(cur_subject_to)+'_multialign')
        
        multialgn(prog_dir,c_align_st_file,c_align_st_mutalign_file)
        with open(c_align_st_mutalign_file) as f:
            multialign_sequence = '\n'.join(f.read().split('\n')[3:]).strip()
        
        match_num = multialign_sequence.split('\n')[2].count('*')
        mismatch_num = int(spacer_length)-match_num        
        f_save.write('\t'.join(line[1:]).strip()+'\t'+c_spacer_sequence.split('\n')[1].strip()+'\t'+str(cur_subject_from+1)+'-'+str(cur_subject_to)+'\t'+expand_sequence+'\t'+str(mismatch_num)+'\n')
        f_save.flush()
    f_save.close()

def identify_st(prog_dir,work_dir,bac_dir,input_file,mismatch='2',coverage='1'):
    st_dir = '%s/%s/self-targeting' % (work_dir, bac_dir)
    mkdir(st_dir)
    spacer_file = '%s/%s/%s_merge.spc' % (work_dir,bac_dir,bac_dir)
    contig_dir = '%s/%s/protein' % (work_dir,bac_dir)
    
    fasta_file = input_file
    spc_file =  '%s/%s/%s_merge.spc' % (work_dir, bac_dir, bac_dir)
    space_blastn_fasta = '%s/%s/%s.merge_spc_blastn_fasta' % (work_dir, bac_dir, bac_dir)
    merge_csp_file = '%s/%s/%s_merge.csp' % (work_dir, bac_dir, bac_dir)
    
    # step1. sequence alignment
    blastn(prog_dir,spc_file, fasta_file, space_blastn_fasta)
    # step2: add genome description and coverage
    genome_def_dict = {}
    get_def_dict(fasta_file, genome_def_dict)
    blastn_add_cov_file = '%s/%s/%s.merge_spc_blastn_fasta_add_cov' % (work_dir, bac_dir, bac_dir)
    add_cov(space_blastn_fasta, blastn_add_cov_file, genome_def_dict)
    # step3: filter mismatch > 2 and coverage < 1
    filter_mismatch_cov_file = '%s/%s/%s.merge_spc_blastn_fasta_add_cov_filter_mismatch_%s_cov_%s' % (work_dir, bac_dir, bac_dir,mismatch,coverage)
    filter_mismatch_and_cov(blastn_add_cov_file, filter_mismatch_cov_file)
    # step4: filter hit in array
    array_info_dict = {}
    get_array_info_dict(merge_csp_file, array_info_dict)
    filter_hit_in_array = filter_mismatch_cov_file + '_hit_not_in_array'
    filter_in_array(prog_dir,filter_mismatch_cov_file, filter_hit_in_array, array_info_dict,st_dir,spacer_file,contig_dir)

def get_faa_prot_dict(faa_file_path,faa_prot_location_dict):
    with open(faa_file_path) as f:
        contents = f.read()
        elems = contents.split('>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            defLine = elem.split('\n')[0]
            cur_genome = defLine.split('_')[0].split('.')[0]
            prot_start = min(int(defLine.split('_')[-3]),int(defLine.split('_')[-2]))
            prot_end = max(int(defLine.split('_')[-3]),int(defLine.split('_')[-2]))
            if cur_genome not in faa_prot_location_dict:
                faa_prot_location_dict[cur_genome] = []
            faa_prot_location_dict[cur_genome].append([prot_start,prot_end])

def check_array_in_prot(genome_id,array_start,array_end,faa_prot_location_dict):
    array_in_prot = "no"
    if genome_id in faa_prot_location_dict.keys():
        prot_locatio_list = faa_prot_location_dict[genome_id]
        for prot_elem in prot_locatio_list:
            prot_start = min(int(prot_elem[0]),int(prot_elem[1]))
            prot_end = max(int(prot_elem[0]),int(prot_elem[1]))
            if int(array_start)>=int(prot_start) and int(array_end)<= int(prot_end):
                array_in_prot = "yes"
                break
    return array_in_prot

def get_spacer_dict(spacer_file,spc_info_dict):
    with open(spacer_file, 'r') as fin:
        content = fin.read()
        elems = content.split('>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            defLine = elem.split('\n')[0]
            genome_id = defLine.split('|')[-1].split('.')[0]
            crispr_locus_id = defLine.split('|')[0].split('.')[0]
            if genome_id not in spc_info_dict:
                spc_info_dict[genome_id] = {}
                spc_info_dict[genome_id][crispr_locus_id] = 1
            else:
                if crispr_locus_id not in spc_info_dict[genome_id]:
                    spc_info_dict[genome_id][crispr_locus_id] = 1
                else:
                    spc_info_dict[genome_id][crispr_locus_id] = spc_info_dict[genome_id][crispr_locus_id]+1

def get_spacer_num(contig_id,crispr_id,spacer_file):
    spacer_num = 0
    ###get the spacer numbers:
    if os.path.exists(spacer_file):
        with open(spacer_file) as f:
            spacers = f.read().strip().split('>')
        if len(spacers)>1:
            for spa in spacers[1:]:
                title = spa.strip().split('\n')[0].strip()
                crispr_cur_id = title.split('|')[0].split('.')[0]
                contig_cur_id = title.split('|')[-1].split('.')[0]
                if (crispr_cur_id == crispr_id) and (contig_cur_id == contig_id):
                    spacer_num = spacer_num+1
    return spacer_num

def get_st_dict(st_result_file,st_info_dict):
    spacer_list = []
    if not os.path.exists(st_result_file):
        return 0
    with open(st_result_file,'r')as fin:
        lines = fin.readlines()
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_genome_id = content[1]
            spacer_info = content[3]
            
            spacer_start = content[3].split('|')[1]
            spacer_end = int(content[3].split('|')[1])+int(content[3].split('|')[2])-1
            
            cur_array_locus_id = content[3].split('|')[0].split('.')[0]
            spacer_location = spacer_start+'-'+str(spacer_end)
            hit_contig_id = content[4]
            protospacer_start = content[-6]
            protospacer_end = content[-5]
            protospacer_location = hit_contig_id+'_'+protospacer_start+'-'+protospacer_end
            if cur_genome_id not in st_info_dict:
                st_info_dict[cur_genome_id] = {}
                st_info_dict[cur_genome_id][cur_array_locus_id] = []
                st_info_dict[cur_genome_id][cur_array_locus_id].append('1')
                st_info_dict[cur_genome_id][cur_array_locus_id].append('1')
                st_info_dict[cur_genome_id][cur_array_locus_id].append(spacer_location)
                st_info_dict[cur_genome_id][cur_array_locus_id].append(protospacer_location)
            else:
                if cur_array_locus_id not in st_info_dict[cur_genome_id].keys():
                    st_info_dict[cur_genome_id][cur_array_locus_id] = []
                    st_info_dict[cur_genome_id][cur_array_locus_id].append('1')
                    st_info_dict[cur_genome_id][cur_array_locus_id].append('1')
                    st_info_dict[cur_genome_id][cur_array_locus_id].append(spacer_location)
                    st_info_dict[cur_genome_id][cur_array_locus_id].append(protospacer_location)
                else:
                    if not spacer_info in spacer_list:
                        st_info_dict[cur_genome_id][cur_array_locus_id][0] = str(int(st_info_dict[cur_genome_id][cur_array_locus_id][0])+1)
                    st_info_dict[cur_genome_id][cur_array_locus_id][1] = str(int(st_info_dict[cur_genome_id][cur_array_locus_id][1]) + 1)
                    st_info_dict[cur_genome_id][cur_array_locus_id][2] = st_info_dict[cur_genome_id][cur_array_locus_id][2]+ '|'+spacer_location
                    st_info_dict[cur_genome_id][cur_array_locus_id][3] = st_info_dict[cur_genome_id][cur_array_locus_id][3]+'|'+protospacer_location
            if not spacer_info in spacer_list:
                spacer_list.append(spacer_info)

def get_cas_prots(cpt_file,genome_id,crispr_array_start,crispr_array_end,array_neigh_range):
    within_20kb_prots = []
    genome_prots = []
    with open(cpt_file,'r')as fin:
        lines = fin.readlines()
        for line in lines:
            content = line.strip('\n').split('\t')
            prot_info = content[0].split('_')
            if '.' in content[0]:
                cur_genome_id = content[0].split('.')[0]
            else:
                cur_genome_id = '_'.join(prot_info[0:-3])
            if (cur_genome_id in genome_id) or (genome_id in cur_genome_id):
                cas_prot = content[1]
                if cas_prot not in genome_prots:
                    genome_prots.append(cas_prot)
                cur_prot_start = min(int(prot_info[-2]),int(prot_info[-3]))
                cur_prot_end = max(int(prot_info[-2]),int(prot_info[-3]))
                array_neighbour_range_start = int(crispr_array_start)-int(array_neigh_range)
                array_neighbour_range_end = int(crispr_array_end)+int(array_neigh_range)
                if cur_prot_start>=array_neighbour_range_start and cur_prot_end<=array_neighbour_range_end:
                    if cas_prot not in within_20kb_prots:
                        within_20kb_prots.append(cas_prot)
    return within_20kb_prots,genome_prots

def getUpStreamProt(faa_file,crispr_start,num_prot,genome_id):
    up_stream_prot_temp = []
    hitFlag = 0
    titles = []
    index = 0
    with open(faa_file, 'r')as fin:
        content = fin.read()
        elems = content.split('>')
        if '' in elems:
            elems.remove('')
        prot_count = 0
        for elem in elems[::-1]:
            cur_genome_id = '_'.join(elem.split('\n')[0].strip().split('_')[0:-3])
            if (genome_id in cur_genome_id) or (cur_genome_id in genome_id):
                cur_prot_start = elem.split('\n')[0].split('_')[-3]
                cur_prot_end = elem.split('\n')[0].split('_')[-2]
                if int(cur_prot_end) < int(crispr_start):
                    if prot_count == num_prot:
                        break
                    else:
                        prot_count = prot_count + 1
                        prot_size = str(int((int(cur_prot_end) - int(cur_prot_start) + 1) / 3))
                        prot_info = '%s-%s|%saa' % (cur_prot_start, cur_prot_end, prot_size)
                        titles.append(prot_size+'aa|up_'+str(index)+'|'+elem.strip().split('\n')[0])
                        index = index+1
                        up_stream_prot_temp.append('>'+elem.strip())
        titles = titles[::-1]
        if len(up_stream_prot_temp)>0:
            num_prior_prot = len(up_stream_prot_temp)
            if len(up_stream_prot_temp)==num_prot:
                up_stream_prot_list = up_stream_prot_temp
            else:
                titles = ['NA']*(num_prot-len(up_stream_prot_temp))+titles
                up_stream_prot_list = up_stream_prot_temp
        else:
            titles = ['NA']*num_prot
            up_stream_prot_list = ['NA']*num_prot
    return [up_stream_prot_list,titles]

def getDownStreamProt(faa_file,crispr_end,num_prot,genome_id):
    down_stream_prot_list = []
    titles = []
    index = 0
    with open(faa_file,'r')as fin:
        content = fin.read()
        elems = content.split('>')
        if '' in elems:
            elems.remove('')
        prot_count = 0
        for elem in elems:
            cur_genome_id = '_'.join(elem.split('\n')[0].strip().split('_')[0:-3])
            if (genome_id in cur_genome_id) or (cur_genome_id in genome_id):
                cur_prot_start = elem.split('\n')[0].split('_')[-3]
                cur_prot_end = elem.split('\n')[0].split('_')[-2]
                if int(cur_prot_start) > int(crispr_end):
                    if prot_count == num_prot:
                        break
                    else:
                        prot_count = prot_count + 1
                        prot_size = str(int((int(cur_prot_end) - int(cur_prot_start) + 1) / 3))
                        titles.append(prot_size+'aa|down_'+str(index)+'|'+elem.strip().split('\n')[0])
                        index = index+1
                        prot_info = '%s-%s|%saa'%(cur_prot_start,cur_prot_end,prot_size)
                        down_stream_prot_list.append('>'+elem.strip())
        if len(down_stream_prot_list)>0:
            num_prior_prot = len(down_stream_prot_list)
            if len(down_stream_prot_list)<num_prot:
                titles = titles+['NA']*(num_prot-len(down_stream_prot_list))
                down_stream_prot_list = down_stream_prot_list
        else:
            down_stream_prot_list = ['NA']*num_prot
            titles = ['NA']*num_prot
    return [down_stream_prot_list,titles]

def pro_annotation(anno_file,annotaion_cas_file,proteins):
    if os.path.exists(anno_file):
        with open(anno_file) as f:
            contents = f.read()
    else:
        contents = ''
    proteins_new = []
    unknown_proteins = []
    for protein_1 in proteins[1]:
        pro = protein_1.strip()
        new_pro = 'NA'+'|'+protein_1
        if os.path.exists(annotaion_cas_file):
            with open(annotaion_cas_file) as f:
                annotations = f.readlines()
            if len(annotations)>0:
                flag = 0
                for an in annotations:
                    an = an.strip().split('\t')
                    prot = '_'.join(an[0].strip().split('_')[0:-1])
                    cas = an[1]
                    if prot in protein_1:
                        new_pro = cas+'|'+protein_1
                        flag = 1
                        break
        proteins_new.append(new_pro)
        if pro != 'NA':
            pro_title = pro.split('|')[-1].strip()
            anno_start = contents.find(pro_title)
            if anno_start!=-1:
                anno_end = contents.find('\n',anno_start)
                annotation = contents[anno_start:anno_end].strip().split('\t')[2].split('.')[0]
                proteins_new.append(annotation)
            else:
                proteins_new.append('NA')
                unknown_proteins.append(new_pro)
        else:
            proteins_new.append('NA')
    if len(unknown_proteins)==0:
        unknown_proteins = ['NA']
    return proteins_new,unknown_proteins

def getSpecifiNumProt(prog_dir,work_dir,bac_dir,genome_id,crispr_array_locus_merge,faa_file_path,crispr_array_start,crispr_array_end,num_prot_near_array,cpt_file):
    protein_cas_dict = {}
    with open(cpt_file,'r') as fin:
        lines = fin.readlines()
        for line in lines:
            content = line.strip('\n').split('\t')
            protein_cas_dict.update({content[0]:line})

    up_down_dir = '%s/%s/up_down' % (work_dir, bac_dir)
    mkdir(up_down_dir)
    array_dir = '%s/%s' % (up_down_dir, genome_id + '_' + crispr_array_locus_merge)
    mkdir(array_dir)
    up_proteins = getUpStreamProt(faa_file_path, crispr_array_start, num_prot_near_array, genome_id)
    down_proteins = getDownStreamProt(faa_file_path, crispr_array_end, num_prot_near_array, genome_id)
    cas_result_file = '%s/up_down_cas' % (up_down_dir)
    faa_file = '%s/protein.faa' % (array_dir)
    with open(faa_file, 'w') as f:
        if up_proteins[0].count('NA')==0:
            f.write('\n'.join(up_proteins[0]).strip()+'\n')
        if down_proteins[0].count('NA')==0:
            f.write('\n'.join(down_proteins[0]).strip())
    
    annotaion_cas_file = '%s/protein_cas' % (array_dir)
    blast_file = os.path.join(array_dir, 'blast.xml')
    anno_file = os.path.join(array_dir, 'annotation')
    xml_file_path = '%s/protein.xml' % (array_dir)
    
    f_annotation_cas = open(annotaion_cas_file,'w')
    for protein in up_proteins[0]+down_proteins[0]:
        protein_id = protein.strip().split('\n')[0].strip('>')
        if protein_id in protein_cas_dict.keys():
            f_annotation_cas.write(protein_cas_dict[protein_id].strip()+'\n')
    f_annotation_cas.close()
    
    rpsblastp11proteins(prog_dir,faa_file,blast_file)
    
    ParseResult(blast_file,anno_file)
    up_proteins_new,up_unknown_proteins = pro_annotation(anno_file, annotaion_cas_file, up_proteins)
    down_proteins_new,down_unknown_proteins = pro_annotation(anno_file,annotaion_cas_file, down_proteins)
    return up_proteins_new,down_proteins_new,up_unknown_proteins+down_unknown_proteins

def check_crispr_type_by_cas_prot_old(cas_prot):
    if len(cas_prot) == 0:
        crispr_type = 'Orphan'
    else:
        crispr_type_list = []
        cas_prot_str = ''.join(cas_prot)
        cas_prot_str = cas_prot_str.lower()
        type_dict = {'I':['cas8a','cas8b','cas8c','cas10d','cas8e','cas8f','cas8u'],
                     'II':['cas9'],'III':['cas10'],
                     'V':['cpf1','cas12a','cas12b','cas12c','cas12d','cas12e','c2c9','c2c10','tnpb','cas14a','cas14b','cas14c','cas14u'],
                     'VI':['cas13a','cas13b','cas13c','cas13d']}
        for key in type_dict:
            for type_item in type_dict[key]:
                if type_item in cas_prot_str:
                    if key not in crispr_type_list:
                        crispr_type_list.append(key)
        if crispr_type_list:
            crispr_type = ','.join(crispr_type_list)
        else:
            crispr_type = 'Unclear'
    return crispr_type

def statis_result(prog_dir,work_dir, bac_dir,input_file,mismatch='2',coverage='1',array_neigh_range='20000',num_prot_near_array='10'):
    result_file_original = '%s/%s/statis_result_original.txt' % (work_dir,bac_dir)
    fout = open(result_file_original,'w')
    header = 'assembly_id\tgenome_id\tgenome_def\tcrispr_array_locus_merge\tcrispr_array_location_merge\tcrispr_locus_id\tcrispr_pred_method\tarray_in_prot\tprot_within_array_'+array_neigh_range+'\tprot_in_genome\tcrispr_type_by_cas_prot\tconsensus_repeat\trepeat_length\tself-targeting_spacer_number\tself-targeting_target_number\tspacer_location\tprotospacer_location\trepeat_type\tunknown_protein_around_crispr\t'
    for i in range(1,int(num_prot_near_array)+1):
        header = header + 'L%s\tL%s_domain\t'%(str(int(num_prot_near_array)+1-i),str(int(num_prot_near_array)+1-i))
    for j in range(1,int(num_prot_near_array)+1):
        header = header + 'R%s\tR%s_domain\t' % (str(j),str(j))
    header = header.rstrip('\t')+'\n'
    fout.write(header)
    spacer_file = '%s/%s/%s_merge.spc' % (work_dir, bac_dir, bac_dir)
    merge_csp_file = '%s/%s/%s_merge.csp' % (work_dir, bac_dir, bac_dir)
    fasta_file_path = input_file
    faa_file_path = '%s/%s/%s.faa' % (work_dir, bac_dir, bac_dir)
    cpt_file_path = '%s/%s/%s.cpt' % (work_dir, bac_dir, bac_dir)
    st_result_file = '%s/%s/%s.merge_spc_blastn_fasta_add_cov_filter_mismatch_%s_cov_%s_hit_not_in_array' % (work_dir, bac_dir, bac_dir,mismatch,coverage)
    assembly_id = os.path.basename(merge_csp_file).rstrip('_merge.csp')
    genome_def_dict= {}
    get_def_dict(fasta_file_path,genome_def_dict)
    faa_prot_location_dict = {}
    get_faa_prot_dict(faa_file_path,faa_prot_location_dict)
    spc_info_dict = {}
    get_spacer_dict(spacer_file,spc_info_dict)
    
    st_info_dict = {}
    get_st_dict(st_result_file,st_info_dict)

    #get the repeat_type_dict
    out_repeat_file = '%s/%s/%s_repeat_type' % (work_dir, bac_dir, bac_dir)
    with open(out_repeat_file) as f:
        repeat_types = f.readlines()[1:]
    
    crispr_repeat_types = {}
    for line in repeat_types:
        line = line.strip().split('\t')
        crispr_id = '_'.join(line[0].split('_')[0:-1])
        type = line[-1].strip()
        if crispr_id not in crispr_repeat_types.keys():
            crispr_repeat_types.update({crispr_id:[]})
            crispr_repeat_types[crispr_id].append(type)
        else:
            crispr_repeat_types[crispr_id].append(type)

    with open(merge_csp_file,'r')as fin:
        lines = fin.readlines()
        for line in lines:
            content = line.strip('\n').split('\t')
            genome_id = content[1]
            genome_def = genome_def_dict[genome_id]
            crispr_array_locus_merge = content[0]
            crispr_array_location_merge = content[3]+'-'+content[4]
            crispr_locus_id = content[2]
            
            crispr_pred_method = content[-1]
            crispr_array_start = min(int(content[3]),int(content[4]))
            crispr_array_end = max(int(content[3]),int(content[4]))
            
            # get protein within 20kb of crispr array
            array_in_prot = check_array_in_prot(genome_id,crispr_array_start,crispr_array_end,faa_prot_location_dict)
            if array_in_prot=='yes':
                continue           
            within_20kb_prots, genome_prots = get_cas_prots(cpt_file_path, genome_id, crispr_array_start, crispr_array_end, array_neigh_range)
            crispr_type_by_cas_prot = check_crispr_type_by_cas_prot(within_20kb_prots)
            
            # get specific number protein upstream and downstream of crispr array
            up_proteins_new, down_proteins_new,unknown_proteins = getSpecifiNumProt(prog_dir,work_dir, bac_dir,genome_id,crispr_array_locus_merge,faa_file_path,crispr_array_start,crispr_array_end,int(num_prot_near_array),cpt_file_path)
            consensus_repeat = content[5]
            repeat_length = ''
            repeat_list = consensus_repeat.split(',')
            for repeat_item in repeat_list:
                if repeat_length == '':
                    repeat_length = str(len(repeat_item))
                else:
                    repeat_length = repeat_length + ',' + str(len(repeat_item))
                
            if genome_id in st_info_dict.keys():
                if crispr_array_locus_merge in st_info_dict[genome_id].keys():
                    st_spacer_num = st_info_dict[genome_id][crispr_array_locus_merge][0]
                    st_target_num = st_info_dict[genome_id][crispr_array_locus_merge][1]
                    spacer_location = st_info_dict[genome_id][crispr_array_locus_merge][2]
                    protospacer_location = st_info_dict[genome_id][crispr_array_locus_merge][3]
                else:
                    st_spacer_num = '0'
                    st_target_num = '0'
                    spacer_location = 'NA'
                    protospacer_location = 'NA'
            else:
                st_spacer_num = '0'
                st_target_num = '0'
                spacer_location = 'NA'
                protospacer_location = 'NA'
            
            #get the repeat type
            if genome_id.split('.')[0]+'_'+crispr_array_locus_merge in crispr_repeat_types.keys():
                repeat_type = crispr_repeat_types[genome_id.split('.')[0]+'_'+crispr_array_locus_merge]
            else:
                repeat_type = ['NA']
            strWrite = assembly_id+'\t'+genome_id+'\t'+genome_def+'\t'+crispr_array_locus_merge+'\t'+crispr_array_location_merge+'\t'+crispr_locus_id+'\t'+crispr_pred_method+'\t'+array_in_prot+'\t'+','.join(within_20kb_prots)+'\t'+','.join(genome_prots)+'\t'+crispr_type_by_cas_prot+'\t'+consensus_repeat+'\t'+str(repeat_length)+'\t'+str(st_spacer_num)+'\t'+str(st_target_num)+'\t'+spacer_location+'\t'+protospacer_location+'\t'+':'.join(repeat_type)+'\t'+','.join(unknown_proteins)+'\t'+'\t'.join(up_proteins_new)+'\t'+'\t'.join(down_proteins_new)+'\n'
            fout.write(strWrite)
    fout.close()

    #add maximum spacer ,correct crispr type and genome cas prots information
    save_add_info_prefix = '%s/%s/crispr_cas_self-targeting_statis_result.txt' % (work_dir,bac_dir)
    st_add_info(result_file_original,save_add_info_prefix,cpt_file_path,spacer_file)

def identify_repeat_type(prog_dir,file,outfile):
    py = "/home/qiusuo/miniconda3/bin/python %s/identify_repeat_type.py"%(prog_dir)
    command = py+" --input "+file+" --output "+outfile + " --prog_dir " + prog_dir
    os.system(command)
    return outfile

def identify_interacting_phage(prog_dir,input_file,outdir,mismatch,coverage):
    cmd_python = '/home/qiusuo/miniconda3/bin/python %s/predict_spacer_interact_phage.py --input %s --output %s --mismatch %s --coverage %s --prog_dir %s'%(prog_dir,input_file,outdir,str(mismatch),str(coverage),prog_dir)
    #print(cmd_python)
    os.system(cmd_python)
    with open(os.path.join(outdir,'script'),'w') as f:
        f.write(cmd_python)

def add_spacer_num(st_result_file,save_file,run_dir,prefix):
    with open(st_result_file) as f:
        contents = f.readlines()
    f_save = open(save_file,'w')
    header = contents[0].strip().split('\t')
    f_save.write('\t'.join(header[0:header.index('repeat_length')]+['spacer_locus_num','spacer_max_num']+header[header.index('repeat_length'):])+'\n')
    f_save.flush()
    header = contents[0].strip().split('\t')
    method_list = ['casfinder','pliercr','crt','crispridentify']
    spacer_dict = {}
    for method in method_list:
        if method!='crt':
            c_spacer_file = os.path.join(run_dir,prefix+'_'+method+'.spc')
        else:
            c_spacer_file = os.path.join(run_dir,'crt_dir','crt.spc')
        if os.path.exists(c_spacer_file):
            with open(c_spacer_file) as f:
                c_contents = f.read()
            for spacer in c_contents.split('>')[1:]:
                c_locus_id = spacer.split('\n')[0].split('|')[0].split('.')[0]
                c_contig_id = spacer.split('\n')[0].split('|')[-1]
                if method not in spacer_dict.keys():
                    spacer_dict.update({method:{}})
                if c_contig_id not in spacer_dict[method].keys():
                    spacer_dict[method].update({c_contig_id:{}})
                if c_locus_id not in spacer_dict[method][c_contig_id].keys():
                    spacer_dict[method][c_contig_id].update({c_locus_id:[]})
                spacer_dict[method][c_contig_id][c_locus_id].append(spacer.split('\n')[0])
    
    for line in contents[1:]:
        line = line.strip().split('\t')
        method = line[header.index('crispr_pred_method')]
        crispr_locus_id = line[header.index('crispr_locus_id')]
        method_list = method.split(',')
        crispr_locus_id_list = crispr_locus_id.split(',')
        contig_id = line[header.index('genome_id')]
        spacer_num_list = []
        for index,method_name in enumerate(method_list):
            if method_name=='CRISPRCasFinder':
                c_name = 'casfinder'
            if method_name=='PILER-CR':
                c_name = 'pliercr'
            if method_name=='CRT':
                c_name = 'crt'
            if method_name=='CRISPRidentify':
                c_name = 'crispridentify'
            crispr_id = crispr_locus_id_list[index]
            c_spacer_num = len(spacer_dict[c_name][contig_id][crispr_id])
            spacer_num_list.append(c_spacer_num)
        max_spacer = round(np.max(spacer_num_list),0)
        f_save.write('\t'.join(line[0:header.index('repeat_length')]+[','.join(list(map(str,spacer_num_list))),str(max_spacer)]+line[header.index('repeat_length'):])+'\n')
        f_save.flush()
    f_save.close()

def getFaaFromGB(input_file,faa_file):   #parse protein from genbank files in phaster
    special_pros = ['capsid','head','plate','tail','coat','portal','holin','integrase','transposase','terminase','protease','lysis','bacteriocin','tRNA']
    records = SeqIO.parse(input_file, "gb")
    counter = 0
    savefile = open(faa_file, 'w')
    for record in records:
        contig_id = record.id
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
                location = min_start+'_'+max_end+'_'+direction
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

def getFnaFromGB(gb_file_path,fa_file_path):
    handle = open(gb_file_path)
    if os.path.exists(fa_file_path):
        os.remove(fa_file_path)
    SeqIO.convert(handle, 'genbank', fa_file_path, 'fasta')

def split_genbank(gb_file,outdir,prefix):
    records = SeqIO.parse(gb_file, "gb")
    gb_faa_file = os.path.join(outdir,prefix+'_gb.faa')
    faa_file = os.path.join(outdir,prefix+'.faa')            
    fa_file = os.path.join(outdir,prefix+'.fa')    
    sequence_dir = os.path.join(outdir,'gb')
    protein_dir = os.path.join(outdir,'protein')
    protein_name_file = os.path.join(outdir,'protein_name_list.txt')
    strain_inf_file = os.path.join(outdir,'strain_info.txt')
    
    mkdir(sequence_dir)
    mkdir(protein_dir)
    f_gb_protein = open(gb_faa_file,'w')
    f_nucl = open(fa_file,'w')
    f_protein = open(faa_file,'w')
    f_protein_name_list = open(protein_name_file,'w')
    f_strain_inf = open(strain_inf_file,'w')
    f_strain_inf.write('strain_id\tdescription\tflag\tlength\n')
    strain_proteins = ''
    for record in records:
        contig_id = record.id
        contig_inf = record.description
        
        contig_file = os.path.join(sequence_dir,contig_id+'.gb')
        SeqIO.write(record,contig_file, "gb")
        contig_faa_file = os.path.join(protein_dir,contig_id+'.faa')
        contig_fa_file = os.path.join(protein_dir,contig_id+'.fa')
        if not os.path.exists(contig_faa_file):
            getFaaFromGB(contig_file, contig_faa_file)
        if not os.path.exists(contig_fa_file):
            getFnaFromGB(contig_file, contig_fa_file)
        with open(contig_faa_file) as fp:
            protein = fp.read().strip()
        proteins = protein.split('\n>')
        for pro in proteins:
            if pro.strip()=='':
                continue
            title = pro.strip().split('\n')[0].strip('>')
            sequence = ''.join(pro.strip().split('\n')[1:]).strip()
            new_title = '>'+title
            f_gb_protein.write(new_title.strip()+'\n'+sequence+'\n')           
            strain_id = title.strip('|').split('|')[0]
            protein_location = title.strip('|').split('|')[2]
            #print(protein_location)
            if str(protein_location).find('+') != -1:
                direction = '+'
            else:
                direction = '-'           
            new_prot_id = strain_id+'_'+protein_location
            ori_prot_id = new_title          
            f_protein.write('>'+new_prot_id+'\n'+sequence.strip()+'\n')
            f_protein_name_list.write(new_prot_id+'\t'+ori_prot_id.strip('>')+'\n')       
        with open(contig_fa_file) as f1:
            nucl = f1.read()
        f_strain_inf.write(contig_id+'\t'+contig_inf+'\tgenbank\t'+str(len(''.join(nucl.split('\n')[1:]).strip()))+'\n')
        f_nucl.write(nucl.strip()+'\n')    
    f_gb_protein.close()
    f_protein.close()
    f_nucl.close()
    f_strain_inf.close()

class MyThread(threading.Thread):
    def __init__(self,work_dir,bac_dir,result_file):
        threading.Thread.__init__(self)
        self.work_dir = work_dir
        self.bac_dir = bac_dir
        self.result_file = result_file
    def run(self):
        identify_selftargeting(self.work_dir,self.bac_dir,self.result_file)

def get_inf(file,outdir):
    with open(file) as f:
        contents = f.read().strip()
    if contents[0]=='>':
        type ='fasta'
    else:
        type="genbank"
    return type


def correct_type(cas_prots,original_type):
    type_dict = {"Cas8a1":["Type I-A","Type I-G"],
        "Cas8a2":["Type I-A"],   
        "Cas8b1":["Type I-B"],                                                                                
        "Cas8c": ["Type I-C","Type I-U"],                                                                                                                                           
        "Cas8e": ["Type I-E"],                                                                                                                                                                                                       
        "Cas8f": ["Type I-F"],                                                                                                                                                                                                                                                                   
        "Cas9":["Type II-A","Type II-B","Type II-C"],                                                                                                                                                                                                                                                                                                                
        "Cas10":["Type III-A","Type III-B","Type III-C","Type III-D"],
        "Cas10d":["Type I-D"],
        "CasB":["Type I-E"],
        "CasD":["Type I-E"],
        "CasC":["Type I-E"],
        "CasE":["Type I-E"],
        "Cas12a":["Type V-A"],
        "Cpf1":["Type V-A"],
        "Cas12b":["Type V-B"],
        "c2c1":["Type V-B"],
        "Cas12c":["Type V-C"],
        "c2c3":["Type V-C"],
        "Cas12d":["Type V-D"],
        "CasY":["Type V-D"],
        "Cas12e":["Type V-E"],
        "Cas12k":["Type V"],
        "CasX":["Type V-E"],
        "c2c4":["Type V-U1"],
        "c2c5":["Type V"],
        "c2c8":["Type V-U2"], 
        "c2c10":["Type V-U3"], 
        "c2c9":["Type V-U4"],  
        "Cas13a":["Type VI-A"],
        "c2c2":["Type VI-A"],
        "Cas13b1":["Type VI-B1"],
        "Cas13b2":["Type VI-B2"],
        "c2c6":["Type VI-B1","Type VI-B2"],
        "Csx27":["Type VI-B1"],
        "Csx28":["Type VI-B2"],
        "Cas13c":["Type VI-C"],
        "c2c7":["Type VI-C"],
        "cas14a":['Type V'],
        "cas14b":['Type V'],
        "cas14c":['Type V'],
        "cas14d":['Type V'],
        "cas14e":['Type V'],
        "cas14f":['Type V'],
        "cas14g":['Type V'],
        "cas14h":['Type V'],
        "cas14u":['Type V'],
        "cas14i":['Type V'],
        "cas14j":['Type V'],
        "cas14k":['Type V'],
        "Cas14c_CAS-V-F":['Type V'],
        "Cas13d":['Type VI-D'],
        "Cas12g":['Type V']
        }
    crispr_type = ""
    for cas_prot in cas_prots:
        for c_cas_prot in type_dict.keys():
            if cas_prot!='c2c1':
                if cas_prot.lower() == c_cas_prot.lower():
                    crispr_type = crispr_type+','+','.join(type_dict[c_cas_prot])
                else:
                    if c_cas_prot in ['Cas8a1','Cas8a2','Cas8b1','Cas13b1','Cas13b2']:
                        if cas_prot.lower() == c_cas_prot.lower()[0:-1]:
                            crispr_type = crispr_type+','+','.join(type_dict[c_cas_prot])
            else:
                crispr_type = crispr_type+',Type V-B'
    crispr_type = (crispr_type+','+original_type).strip(',')
    crispr_type = ','.join(list(set(crispr_type.split(','))))
    if 'Unclear' in crispr_type and ('Type' in crispr_type):
        crispr_type = crispr_type.replace('Unclear','').strip(',')
    return crispr_type

def st_add_info(st_file,save_file,cpt_file,spacer_file):
    with open(spacer_file) as f:
        contents = f.read()
    spacer_info_dict = {}
    for spacer in contents.split('>')[1:]:
        spacer_id = spacer.split('\n')[0].strip()
        spacer_method = spacer_id.split('|')[-1]
        crispr_locus = spacer_id.split('|')[0].split('.')[0]
        contig_id = spacer_id.split('|')[3]
        if contig_id not in spacer_info_dict.keys():
            spacer_info_dict.update({contig_id:{}})
        if crispr_locus not in spacer_info_dict[contig_id].keys():
            spacer_info_dict[contig_id].update({crispr_locus:{}})
        for method in spacer_method.split(','):
            if method not in spacer_info_dict[contig_id][crispr_locus].keys():
                spacer_info_dict[contig_id][crispr_locus].update({method:[]})
            spacer_info_dict[contig_id][crispr_locus][method].append(spacer_id)
    
    genome_cas_prots = []
    with open(cpt_file) as f:
        contents = f.readlines()
    for line in contents:
        line = line.strip().split('\t')
        cas_prot = line[1]
        if cas_prot not in genome_cas_prots:
            genome_cas_prots.append(cas_prot)
    
    fout1 = open(save_file,'w')
    with open(st_file,'r')as fin:
        lines = fin.readlines()
        header_list = lines[0].strip('\n').split('\t')
        header = header_list
        header_write_list = header_list[0:(header_list.index("repeat_type")+1)]+['spacer_locus_num','spacer_num','correct_crispr_type','genome_cas_prots']+header_list[(header_list.index("repeat_type")+1):]
        fout1.write('\t'.join(header_write_list)+'\n')      
        for line in lines[1:]:
            line = line.strip('\n').split('\t')               
            crispr_type = line[header.index('crispr_type_by_cas_prot')]
            if 'prot_within_array_20000' in header:
                cas_prot = line[header.index('prot_within_array_20000')].split(',')
            else:
                cas_prot = line[header.index('prot_within_array_20kb')].split(',')
            att_name_list = ['DEDDh','WYL','RT',
            'DinG','csa3','TnsE_C','RT','PrimPol','PD-DExK','cas1','cas2','cas4']
            pred_method = line[header.index('crispr_pred_method')]  
            repeat_sequence = line[header.index('consensus_repeat')].split(',')         
            contig_id = line[1]
            genome_id = line[0]
            crispr_array_locus_merge = line[header.index('crispr_array_locus_merge')]
            crispr_array_location_merge = line[header.index('crispr_array_location_merge')]
            method_list = line[header.index('crispr_pred_method')]
            array_in_prot = line[header.index('array_in_prot')]
            crispr_array_start = crispr_array_location_merge.split('-')[0]
            crispr_array_end = crispr_array_location_merge.split('-')[1]        
            diff_cas = list(set(cas_prot).difference(set(att_name_list)))
            
            c_spacer_list = []
            for method in method_list.split(','):
                c_spacer_num = len(spacer_info_dict[contig_id][crispr_array_locus_merge][method])
                c_spacer_list.append(c_spacer_num)
            spacer_num = np.max(c_spacer_list)
            spacer_locus_num = ','.join(list(map(str,c_spacer_list)))
            
            up_down_dir = os.path.join('/'.join(st_file.split('/')[0:-1]),'up_down')
            up_down_save_dir = '%s/%s' % (up_down_dir, contig_id.split('.')[0] + '_' + crispr_array_locus_merge)
            correct_crispr_type = crispr_type           
            if len(cas_prot)>0:
                if len(diff_cas)==0:
                    if len(list(set(cas_prot).intersection(set(['cas1','cas2','cas4']))))>0:
                        correct_crispr_type = 'Unclear'
                    else:
                        correct_crispr_type = 'Orphan'
                else:
                    correct_crispr_type = correct_type(cas_prot,crispr_type)                    
            correct_crispr_type = correct_crispr_type.strip().replace(' ','')
            strwrite_list = line[0:(header_list.index("repeat_type")+1)]+[str(spacer_locus_num),str(spacer_num),correct_crispr_type,','.join(genome_cas_prots)]+line[(header_list.index("repeat_type")+1):]
            fout1.write('\t'.join(strwrite_list)+'\n')  
    fout1.close()

def mkdir(dirPath):
    cmd_mkdir = 'mkdir -p %s'%dirPath
    os.system(cmd_mkdir)


def PredProphageByDBScan(prog_dir,bacSequence,workDir):
    cmd_python = '/home/qiusuo/miniconda3/bin/python %s/dbscan-swa.py --input %s --output %s --add_annotation none --prog_dir %s' % (prog_dir,bacSequence, workDir,prog_dir)
    with open(os.path.join(workDir,'script'),'w') as f:
        f.write(cmd_python)
    os.system(cmd_python)

def get_strain_sequence(prog_dir,input_file,outdir,prefix):
    faa_file = os.path.join(outdir,prefix+'.faa')   
    faaFilePrefix = os.path.join(outdir,prefix)    
    with open(input_file) as f:
        contents = f.read().strip()
    if '>' in contents[0]:
        pred_orf(prog_dir,input_file,faaFilePrefix)
        sequence_dir = os.path.join(outdir,'protein')
        mkdir(sequence_dir)
        strain_inf_file = os.path.join(outdir,'strain_info.txt')
        f_result = open(strain_inf_file,'w')
        f_result.write('strain_id\tdescription\tflag\tlength\n')     
        for strain in contents.split('\n>'):
            strain_title = strain.split('\n')[0].strip()
            strain_id = strain_title.split()[0].strip('>')               
            if '>' in strain_id or ('<') in strain_id:
                print('sequence title error! please remove special characters')
                sys.exit(-1)
            
            strain_def = strain_title.strip('>')
            f_result.write(strain_id+'\t'+strain_def+'\tfasta\t'+str(len(''.join(strain.split('\n')[1:]).strip()))+'\n')
            f_result.flush()
            c_contig_fa_file = os.path.join(sequence_dir,strain_id+'.fa')
            with open(c_contig_fa_file,'w') as f:
                f.write('>'+strain.strip('>').strip()+'\n')
    else:
        
        split_genbank(input_file,outdir,prefix)
        fa_file = os.path.join(outdir,prefix+'.fa')
        faaFilePrefix = os.path.join(outdir,prefix)
        if os.path.getsize(faa_file)==0:
            pred_orf(prog_dir,fa_file,faaFilePrefix)

class MyThread_prophage(threading.Thread):
    def __init__(self,prog_dir,file,outdir):
        threading.Thread.__init__(self)
        self.prog_dir = prog_dir
        self.file = file
        self.outdir = outdir
    def run(self):
        PredProphageByDBScan(self.prog_dir,self.file,self.outdir)

def predict_bac_prophage(prog_dir,strain_inf_file,work_dir,outdir,thread_num=3):
    bac_name = work_dir.split('/')[-1]
    tsk = []
    with open(strain_inf_file) as f:
        contents = f.readlines()
    for line in contents[1:]:
        line = line.strip().split('\t')
        strain_id = line[0]
        flag = line[2]
        protein_file = os.path.join(work_dir,'protein',strain_id+'.faa')
        if flag=='genbank':
            if os.path.getsize(protein_file)>0:
                contig_file = os.path.join(work_dir,'gb',strain_id+'.gb')
            else:
                contig_file = os.path.join(work_dir,'protein',strain_id+'.fa')
        else:
            contig_file = os.path.join(work_dir,'protein',strain_id+'.fa')
        c_outdir = os.path.join(outdir,strain_id)
        mkdir(c_outdir)
        try:
            tsk.append(MyThread_prophage(prog_dir,contig_file,c_outdir))
        except:
            print('Error: unable to start thread')
    
    for t in tsk:
        t.start()
        while True:
            if (len(threading.enumerate()) <= int(thread_num)):
                break
    for t in tsk:
        t.join()
    
    #collect prophage protein and nucl
    prophage_pro_file = os.path.join(outdir,'complete_prophage_protein.faa')
    prophage_file = os.path.join(outdir,'complete_prophage_nucl.fa')
    prophage_summary_file = os.path.join(outdir,'complete_prophage_summary.txt')
    prophage_detail_file = os.path.join(outdir,'complete_prophage_detail.txt')
    
    f_pro = open(prophage_pro_file,'w')
    f_nucl = open(prophage_file,'w')
    f_summary = open(prophage_summary_file,'w')    
    f_detail = open(prophage_detail_file,'w')    
    header_flag = 0
    f_detail.write('''The following contents displays predicted prophage regions
first line of each prophage describes the prophage information and the following lines describe the proteins and homology proteins in uniprot database
bacteria_id bac_def genome_size prophage_start  prophage_end    key_proteins    best_hit_species    CDS_number  attl_region attr_region
''')
    for line in contents[1:]:
        line = line.strip().split('\t')
        strain_id = line[0]
        c_outdir = os.path.join(outdir,strain_id)
        pro_file = os.path.join(c_outdir,'bac_DBSCAN-SWA_prophage.faa')
        nucl_file = os.path.join(c_outdir,'bac_DBSCAN-SWA_prophage.fna')
        summary_file = os.path.join(c_outdir,'bac_DBSCAN-SWA_prophage_summary.txt')
        detail_file = os.path.join(c_outdir,'bac_DBSCAN-SWA_prophage.txt')
        if os.path.exists(pro_file):
            with open(pro_file) as f:
                pro_sequence = f.read().strip()
            if os.path.exists(nucl_file):
                with open(nucl_file) as f:
                    nucl_sequence = f.read().strip()
                f_pro.write(pro_sequence+'\n')
                f_nucl.write(nucl_sequence+'\n')
                
                with open(summary_file) as f:
                    summary_contents = f.readlines()
                summary_header = summary_contents[0].strip().split('\t')
                if header_flag == 0:
                    header_flag = 1
                    f_summary.write('prefix\t'+'\t'.join(summary_header)+'\n')
                for prophage_info in summary_contents[1:]:
                    prophage_info = prophage_info.strip().split('\t')
                    f_summary.write(bac_name+'\t'+'\t'.join(prophage_info).strip()+'\n')
                    f_summary.flush() 

                with open(detail_file) as f:
                    contents = f.read().strip()
                f_detail.write('\n'.join(contents.split('\n')[3:]).strip()+'\n')
                f_detail.flush()
    
    f_pro.close()
    f_nucl.close()
    f_summary.close()
    f_detail.close()

def predict_anti(prog_dir,protein_file,prophage_summary_file,st_info_file,outdir,identity,coverage,up_num,anti_flag):
    script = "%s/anti.py"%(prog_dir)
    command = "/home/qiusuo/miniconda3/bin/python %s --input_protein %s --input_prophage %s --input_st_info %s --output %s --identity %s --coverage %s --up_num %s --flag %s --prog_dir %s"%(script,protein_file,prophage_summary_file,st_info_file,outdir,str(identity),str(coverage),str(up_num),anti_flag,prog_dir)
    os.system(command)
    with open(os.path.join(outdir,'script'),'w') as f:
        f.write(command) 

def statis_anti(anti_info_file,crispr_st_info_file,detail_st_file,save_file):   
    spacer_info_dict = {}
    with open(detail_st_file) as f:
        detail_st_contents = f.readlines()
    header = detail_st_contents[0].strip().split('\t')    
    
    for line in detail_st_contents[1:]:
        line = line.strip().split('\t')
        spacer_id = line[header.index('spacer_id')]
        spacer_start = spacer_id.split('|')[1]
        spacer_length = spacer_id.split('|')[2]
        spacer_end = int(spacer_start)+int(spacer_length)-1
        spacer_location = str(spacer_start)+'-'+str(spacer_end)
        contig_id = line[header.index('bac_id')].split('.')[0]
        if contig_id not in spacer_info_dict.keys():
            spacer_info_dict.update({contig_id:{}})
        spacer_info_dict[contig_id].update({spacer_location:spacer_id})      

    with open(crispr_st_info_file) as f:
        crispr_contents = f.readlines()
    crispr_header = crispr_contents[0].strip().split('\t')   
    crispr_spacer_target_info_dict = {}
    protospacer_location_dict = {}
    genome_crispr_type_list = []
    genome_st_flag = 0
    header = crispr_contents[0].strip().split('\t')
 
    for line in crispr_contents[1:]:
        line = line.strip().split('\t')
        correct_type = line[header.index('correct_crispr_type')]
        crispr_locus_id = line[header.index('crispr_array_locus_merge')]
        crispr_neigh_cas = line[header.index('prot_within_array_20000')]
        contig_id = line[header.index('genome_id')]
        st_number = line[header.index('self-targeting_target_number')]
        spacer_location_list = line[header.index('spacer_location')].split('|')
        protospacer_location_list = line[header.index('protospacer_location')].split('|')
        
        spacer_number = line[header.index('spacer_num')]
        
        if correct_type not in genome_crispr_type_list:
            genome_crispr_type_list.append(correct_type)
        if contig_id not in crispr_spacer_target_info_dict.keys():
            crispr_spacer_target_info_dict.update({contig_id:[]})
        
        crispr_spacer_target_info_dict[contig_id].append([crispr_locus_id,crispr_neigh_cas,correct_type,spacer_number])
        genome_st_flag = genome_st_flag+int(st_number)
        
        if line[header.index('spacer_location')]!='NA':
            for index,spacer_location in enumerate(spacer_location_list):
                spacer_id = spacer_info_dict[contig_id][spacer_location]
                protospacer_location = protospacer_location_list[index]
                protospacer_contig_id = '_'.join(protospacer_location.split('_')[0:-1]).split('.')[0]
                protospacer_location = protospacer_location.split('_')[-1]
                if protospacer_contig_id not in protospacer_location_dict.keys():
                    protospacer_location_dict.update({protospacer_contig_id:{}})            
                if protospacer_location not in protospacer_location_dict[protospacer_contig_id].keys():
                    protospacer_location_dict[protospacer_contig_id].update({protospacer_location:[]})           
                protospacer_location_dict[protospacer_contig_id][protospacer_location].append([spacer_id,spacer_location,contig_id,crispr_locus_id,crispr_neigh_cas,correct_type,spacer_number])
                   
    genome_crispr_type = ','.join(genome_crispr_type_list)
    if not os.path.exists(anti_info_file):
        print('No anti-CRISPR protein detected!')
        return 0
    
    with open(anti_info_file) as f:
        contents = f.readlines()
    anti_header = contents[0].strip().split('\t')
    header_flag = 0
    f_save = open(save_file,'w')
    for line in contents[1:]:
        line = line.strip().split('\t')
        anti_id = line[0]
        anti_contig_id = '_'.join(anti_id.split('_')[0:-3]).split('.')[0]
        anti_prophage_region = line[anti_header.index('Acr_candidate_in_prophage')]       
        contig_crispr_type_list = 'NA'
        if anti_contig_id in crispr_spacer_target_info_dict.keys():
            contig_crispr_type_list = list(set([item[2] for item in crispr_spacer_target_info_dict[anti_contig_id]]))
            contig_crispr_type_list = ','.join(contig_crispr_type_list)                       
        anti_contig_protospacer_number = 0 
        st_target_prophage_spacer_list = []
        st_target_prophage_protospacer_list = []      
        st_target_prophage_crispr_type_list = []
        if anti_contig_id in protospacer_location_dict.keys():
            anti_contig_protospacer_number = len(protospacer_location_dict[anti_contig_id])
            if anti_prophage_region!='No':
                anti_prophage_start = anti_prophage_region.split('-')[0]
                anti_prophage_end = anti_prophage_region.split('-')[1]
                for protospacer_location in protospacer_location_dict[anti_contig_id].keys():
                    protospacer_start = protospacer_location.split('-')[0]
                    protospacer_end = protospacer_location.split('-')[1]
                    flag = 'no'
                    if int(protospacer_start) >= int(anti_prophage_start) and (int(protospacer_end) <= int(anti_prophage_end)):
                        flag = 'yes'
                    elif int(protospacer_start) <= int(anti_prophage_start) and (int(protospacer_end) >= int(anti_prophage_start)):
                        flag = 'yes'
                    elif int(protospacer_start) <= int(anti_prophage_end) and (int(protospacer_end) >= int(anti_prophage_end)):
                        flag = 'yes'
                    if flag == 'yes':    
                        spacer_info = protospacer_location_dict[anti_contig_id][protospacer_location]
                        spacer_list = [item[0] for item in spacer_info]
                        st_target_prophage_spacer_list = st_target_prophage_spacer_list+spacer_list
                        crispr_type_list = [item[4] for item in spacer_info]
                        st_target_prophage_crispr_type_list = st_target_prophage_crispr_type_list+crispr_type_list
                        #st_target_prophage_protospacer_list.append(anti_contig_id+'_'+protospacer_start+'-'+protospacer_end)
                        st_target_prophage_protospacer_list = st_target_prophage_protospacer_list+[anti_contig_id+'_'+protospacer_start+'-'+protospacer_end]*len(spacer_list)
        if len(st_target_prophage_spacer_list)>0:
            st_target_prophage_flag = 'yes'
        else:
            st_target_prophage_flag = 'no'
        
        if header_flag == 0:
            header_flag = 1
            f_save.write('\t'.join(anti_header).strip()+'\tgenome_crispr_type\tgenome_self-targeting_number\tanti_contig_protospacer_number\tprotospacer_in_prophage\tself-targeting_spacer\tself-targeting_protospacer\tself-targeting_crispr_type\n')
            f_save.flush()        
        
        if len(st_target_prophage_spacer_list)>0:
            st_target_prophage_spacer = '&&'.join(st_target_prophage_spacer_list)
        else:
            st_target_prophage_spacer = 'NA'
        
        if len(st_target_prophage_protospacer_list)>0:
            st_target_prophage_protospacer = '&&'.join(st_target_prophage_protospacer_list)
        else:
            st_target_prophage_protospacer = 'NA'
        
        if len(st_target_prophage_crispr_type_list)>0:
            st_target_prophage_crispr_type = ','.join(list(set(st_target_prophage_crispr_type_list))).strip(',')
        else:
            st_target_prophage_crispr_type = 'NA'        
        f_save.write('\t'.join(line).strip()+'\t'+genome_crispr_type+'\t'+str(genome_st_flag)+'\t'+str(anti_contig_protospacer_number)+'\t'+st_target_prophage_flag+'\t'+st_target_prophage_spacer+'\t'+st_target_prophage_protospacer+'\t'+st_target_prophage_crispr_type+'\n')
        f_save.flush()
    f_save.close()

def predict_effector_protein(prog_dir,crispr_st_info_file,outdir,att_spacer_num,att_min_repeat_length,att_max_repeat_length,att_protein_size,att_neigh_distance,repeat_identity,repeat_coverage,homo_identity,homo_coverage):
    script = "%s/predict_novel_effector_protein.py"%(prog_dir)
    command = "/home/qiusuo/miniconda3/bin/python %s --input %s --output %s \
    --att_spacer_num %s --att_min_repeat_length %s --att_max_repeat_length %s \
    --att_protein_size %s --att_neigh_distance %s \
    --id %s --cov %s --identity %s --coverage %s --prog_dir %s"%(script,crispr_st_info_file,outdir,att_spacer_num,att_min_repeat_length,att_max_repeat_length,att_protein_size,att_neigh_distance,repeat_identity,repeat_coverage,homo_identity,homo_coverage,prog_dir)
    os.system(command)
    with open(os.path.join(outdir,'effector_protein_script'),'w') as f:
        f.write(command)

def predict(outdir,input_file,prog_dir,protein_file,strain_inf_file,mismatch='2',coverage='1',array_neigh_range='20000',num_prot_near_array='10',annotation_modules='anti,self-targeting,MGE-targeting,prophage',anti_flag='known',prophage_thread_num=3,anti_identity=0.4,anti_coverage=0.7,anti_up_num=3,anti_prot_size=400,phage_mismatch=2,phage_coverage=0.9,att_spacer_num=2,att_min_repeat_length=18,att_max_repeat_length=45,att_protein_size=500,att_neigh_distance=5,repeat_identity=0.9,repeat_coverage=0.9,homo_identity=0.3,homo_coverage=0.6):
    print('att_spacer_num:',att_spacer_num)
    work_dir = '/'.join(outdir.split('/')[0:-1])
    bac_dir = outdir.split('/')[-1]

    #step1:predict crispr, and annotation cas prot and annotate known anti
    muthread1 = threading.Thread(target=identify_crispr_array, args=(prog_dir,work_dir, bac_dir,input_file,strain_inf_file))
    muthread2 = threading.Thread(target=annotation_cas_prot, args=(prog_dir,work_dir, bac_dir,protein_file))
    muthread1.start()
    muthread2.start()
    muthread1.join()
    muthread2.join()

    
    #step2: identify self-targeting and interacting phage
    int_phage_dir = '%s/%s/interacting_phage' % (work_dir, bac_dir)
    mkdir(int_phage_dir)
    merge_spcfile = '%s/%s/%s_merge.spc' % (work_dir, bac_dir,bac_dir)

    detail_st_file = '%s/%s/self-targeting/spacer_target_hit_result.txt' % (work_dir, bac_dir)
    if 'self-targeting' in annotation_modules:
        muthread1 = threading.Thread(target=identify_st,args=(prog_dir,work_dir, bac_dir,input_file))
    if 'MGE-targeting' in annotation_modules:   
        muthread2 = threading.Thread(target=identify_interacting_phage,args=(prog_dir,merge_spcfile, int_phage_dir,phage_mismatch,phage_coverage))
    if 'self-targeting' in annotation_modules:
        muthread1.start()
    if 'MGE-targeting' in annotation_modules:
        muthread2.start()
    if 'self-targeting' in annotation_modules:
        muthread1.join()
    if 'MGE-targeting' in annotation_modules:
        muthread2.join()
    
    #step3: statis crispr result and predict prophage
    prophage_dir = '%s/%s/prophage' % (work_dir, bac_dir)
    crispr_st_info_file = '%s/%s/crispr_cas_self-targeting_statis_result.txt' % (work_dir,bac_dir)
    mkdir(prophage_dir)
    if 'prophage' in annotation_modules:
        muthread1 = threading.Thread(target=predict_bac_prophage,args=(prog_dir,strain_inf_file,outdir,prophage_dir,prophage_thread_num))
    
    muthread2 = threading.Thread(target=statis_result,args=(prog_dir,work_dir,bac_dir,input_file,mismatch,coverage,array_neigh_range,num_prot_near_array))
    if 'prophage' in annotation_modules:
        muthread1.start()
    muthread2.start()
    if 'prophage' in annotation_modules:
        muthread1.join()
    muthread2.join()

    #step4:predcit anti and novel effector protein
    anti_outdir = '%s/%s/anti' % (work_dir, bac_dir)
    mkdir(anti_outdir)
    prophage_summary_file = '%s/%s/prophage/complete_prophage_summary.txt' % (work_dir, bac_dir)
    if 'anti' in annotation_modules:
        muthread1 = threading.Thread(target=predict_anti,args=(prog_dir,protein_file,prophage_summary_file,crispr_st_info_file,anti_outdir,anti_identity,anti_coverage,anti_up_num,anti_flag))
    if 'effector' in annotation_modules:
        muthread2 = threading.Thread(target=predict_effector_protein,args=(prog_dir,crispr_st_info_file,outdir,att_spacer_num,att_min_repeat_length,att_max_repeat_length,att_protein_size,att_neigh_distance,repeat_identity,repeat_coverage,homo_identity,homo_coverage))
    if 'anti' in annotation_modules:
        muthread1.start()
    if 'effector' in annotation_modules:
        muthread2.start()
    if 'anti' in annotation_modules:
        muthread1.join()
    if 'effector' in annotation_modules:
        muthread2.join()
    
    #step5:statis anti summary info
    anti_info_file = os.path.join(anti_outdir,'acr_candidate_info_integrated_three_methods_check_acr_in_prophage.txt')
    save_anti_file = '%s/%s/anti_statis_result.txt' % (work_dir,bac_dir)
    if 'anti' in annotation_modules:
        statis_anti(anti_info_file,crispr_st_info_file,detail_st_file,save_anti_file)


def delete_temp_files(outdir,mismatch,coverage):
    #sequence
    protein_dir = os.path.join(outdir,'protein')
    gb_dir = os.path.join(outdir,'gb')
    
    #self-targeting
    temp_fil1 = os.path.join(outdir,'%s.merge_spc_blastn_fasta'%outdir.split('/')[-1])
    temp_file2 = os.path.join(outdir,'%s.merge_spc_blastn_fasta_add_cov'%outdir.split('/')[-1])
    temp_fil3 = os.path.join(outdir,'%s.merge_spc_blastn_fasta_add_cov_filter_mismatch_%s_cov_%s'%(outdir.split('/')[-1],str(mismatch),str(coverage)))
    
    #cas
    xml_file = os.path.join(outdir,'%s.xml'%outdir.split('/')[-1])
    files = [protein_dir,gb_dir,temp_fil1,temp_file2,temp_fil3,xml_file]
    for file in files:
        command = "rm -rf %s"%file
        os.system(command)
   
    #crt
    crt_dir = os.path.join(outdir,'crt_dir')
    for file in os.listdir(crt_dir):
        c_file = os.path.join(crt_dir,file)
        if os.path.isdir(c_file):
            command = "rm -rf %s"%file
            os.system(command)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='absPath of genome file.\n')
    parser.add_argument('--output', help='the out dir.\n')
    parser.add_argument('--mismatch', help='value of self-targeting mismatch.\n',default='2')
    parser.add_argument('--cov', help='value of self-targeting coverage.\n', default='1')
    parser.add_argument('--range', help='array neighbour range.\n', default='20000')
    parser.add_argument('--neighbor', help='array neighbour number.\n', default='10')
    parser.add_argument('--anti_identity', help='array neighbour number.\n', default='0.4')
    parser.add_argument('--anti_coverage', help='array neighbour number.\n', default='0.7')
    parser.add_argument('--anti_up_num', help='array neighbour number.\n', default='3')
    parser.add_argument('--phage_mismatch', help='array neighbour number.\n', default='2')
    parser.add_argument('--phage_coverage', help='array neighbour number.\n', default='0.7')
    parser.add_argument('--prophage_thread_num', help='array neighbour number.\n', default='3')
    parser.add_argument('--anti_prot_size', help='array neighbour number.\n', default='3')
    parser.add_argument('--repeat_identity',help='short-blastn identity.\n',default='0.9')
    parser.add_argument('--repeat_coverage',help='short-blastn coverage.\n',default='0.9')
    parser.add_argument('--homo_identity', help='diamond blastp identity.\n', default='0.3')
    parser.add_argument('--homo_coverage',help='diamond blastp coverage.\n',default='0.6')
    parser.add_argument('--att_spacer_num',help='diamond blastp coverage.\n',default='2')
    parser.add_argument('--att_min_repeat_length',help='diamond blastp coverage.\n',default='18')
    parser.add_argument('--att_max_repeat_length',help='diamond blastp coverage.\n',default='45')
    parser.add_argument('--att_protein_size',help='diamond blastp coverage.\n',default='500')
    parser.add_argument('--att_neigh_distance',help='diamond blastp coverage.\n',default='5')
    parser.add_argument('--annotation',help='diamond blastp coverage.\n')
    parser.add_argument('--anti_flag',help='diamond blastp coverage.\n')
    parser.add_argument('--prog_dir', help='directory of program.\n')
    args = parser.parse_args()
    if args.input:
        input_file = args.input
    if args.output:
        outdir = args.output
    if args.prog_dir:
        prog_dir = args.prog_dir
    if args.mismatch:
        mismatch = args.mismatch
    else:
        mismatch = 2
    if args.cov:
        coverage = args.cov
    else:
        coverage = 1
    if args.range:
        range_prot_array_neighbour = args.range    
    else:
        range_prot_array_neighbour = 20000
    if args.neighbor:
        num_prot_array_neighbour = args.neighbor
    else:
        num_prot_array_neighbour = 5
    if args.anti_identity:
        anti_identity = args.anti_identity
    else:
        anti_identity = 0.4
    if args.anti_coverage:
        anti_coverage = args.anti_coverage
    else:
        anti_coverage = 0.7
    if args.anti_up_num:
        anti_up_num = args.anti_up_num
    else:
        anti_up_num = 3
    if args.anti_prot_size:
        anti_prot_size = args.anti_prot_size
    else:
        anti_prot_size = 400
    if args.phage_mismatch:
        phage_mismatch = args.phage_mismatch
    else:
        phage_mismatch = 2    
    if args.phage_coverage:
        phage_coverage = args.phage_coverage
    else:
        phage_coverage = 0.9 
    if args.prophage_thread_num:
        prophage_thread_num = args.prophage_thread_num
    else:
        prophage_thread_num = 3   
    if args.repeat_identity:
        repeat_identity = args.repeat_identity
    else:
        repeat_identity = 0.9
    if args.repeat_coverage:
        repeat_coverage = args.repeat_coverage
    else:
        repeat_coverage = 0.9
    if args.homo_identity:
        homo_identity = args.homo_identity
    else:
        homo_identity = 0.3
    if args.homo_coverage:
        homo_coverage = args.homo_coverage
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
    
    if args.repeat_identity:
        repeat_identity = args.repeat_identity
    else:
        repeat_identity = 0.9

    if args.repeat_coverage:
        repeat_coverage = args.repeat_coverage
    else:
        repeat_coverage = 0.9

    if args.homo_identity:
        homo_identity = args.homo_identity
    else:
        homo_identity = 0.3

    if args.homo_coverage:
        homo_coverage = args.homo_coverage
    else:
        homo_coverage = 0.6

    if args.annotation:
        annotation_modules = args.annotation
    else:
        annotation_modules = 'anti,self-targeting,MGE-targeting,prophage'
    if args.anti_flag:
        anti_flag = args.anti_flag
    else:
        anti_flag = 'known'
    
    prefix = outdir.split('/')[-1]
    
    get_strain_sequence(prog_dir,input_file,outdir,prefix)
    strain_inf_file = os.path.join(outdir,'strain_info.txt')
    protein_file = os.path.join(outdir,prefix+'.faa')
    fa_file = os.path.join(outdir,prefix+'.fa')
    if os.path.exists(fa_file):
        input_file = fa_file       
    print('att_spacer_num',att_spacer_num)
    predict(outdir,input_file,prog_dir,protein_file,strain_inf_file,mismatch,coverage,range_prot_array_neighbour,num_prot_array_neighbour,annotation_modules,anti_flag,prophage_thread_num,anti_identity,anti_coverage,anti_up_num,anti_prot_size,phage_mismatch,phage_coverage,att_spacer_num,att_min_repeat_length,att_max_repeat_length,att_protein_size,att_neigh_distance,repeat_identity,repeat_coverage,homo_identity,homo_coverage)
    delete_temp_files(outdir,mismatch,coverage)