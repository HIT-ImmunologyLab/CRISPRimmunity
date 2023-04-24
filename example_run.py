#!/usr/bin
#-*-coding:utf-8-*-

import os

import predict_novel_effector_protein

def mkdir(dirPath):
    command = "mkdir -p %s"%(dirPath)
    os.system(command)

def run_novel_effector(input_file,res_dir,prog_dir):
    command = "python predict_crispr_cas_self_targeting.py --input %s --output %s --annotation 'effector' --prog_dir %s"%(input_file,res_dir,prog_dir)
    print(command)
    os.system(command)

def run_Acr_pred(input_file,res_dir,prog_dir):
    command = "python predict_crispr_cas_self_targeting.py --input %s --output %s --annotation 'anti,self-targeting,prophage' --anti_flag novel --prog_dir %s"%(input_file,res_dir,prog_dir)
    print(command)
    os.system(command)

def run_IME(input_file,res_dir,prog_dir):
    command = "python predict_crispr_cas_self_targeting.py --input %s --output %s --annotation 'self-targeting,MGE-targeting,prophage,anti' --prog_dir %s"%(input_file,res_dir,prog_dir)
    print(command)
    os.system(command)

if __name__=="__main__":
    prog_dir = "/zrom1/zfx/CRISPRimmunity/github"

    ## run example of Acr prediction
    input_file = "%s/example/sequence/Staphylococcus-schleiferi-strain-5909-02.gb" % (os.getcwd())
    res_dir = "%s/example/result/Acr" % (os.getcwd())
    mkdir(res_dir)
    run_Acr_pred(input_file,res_dir,prog_dir)

    ## run example of identification of novel class 2 effector protein
    input_file = "%s/example/sequence/Armatimonadetes-bacterium-isolate-ATM2-J3.gb"%(os.getcwd())
    res_dir = "%s/example/result/Effector"%(os.getcwd())
    mkdir(res_dir)
    run_novel_effector(input_file, res_dir,prog_dir)

    ## run example of annotation of important molecular events
    input_file = "%s/example/sequence/Staphylococcus-schleiferi-strain-5909-02.gb" % (os.getcwd())
    res_dir = "%s/example/result/IME" % (os.getcwd())
    mkdir(res_dir)
    run_IME(input_file, res_dir, prog_dir)
