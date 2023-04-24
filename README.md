# CRISPRimmunity: the tool for CRISPR-associated Important Molecular events and Modulators Used in geNome edIting Tool identifYing



## Table of Contents

- [Background](#Background)

- [Requirements](#Requirements)
- [Install](#Install)
- [Usage](#Usage)
- [Example](#Example)
- [Contributors](#Contributors)
- [License](#License)

<p id="Background"></p>

## Background

The CRISPR-Cas system is a highly adaptive and RNA-guided immune system found in bacteria and archaea that has applications as a genome editing tool and is a valuable system for studying co-evolutionary dynamics of bacteriophage interactions. This repository presents CRISPRimmunity, a new tool designed for Acr prediction, identification of novel class 2 CRISPR-Cas loci, and dissection of CRISPR-associated important molecular events. CRISPRimmunity is built on a series of CRISPR-oriented databases that offer a comprehensive co-evolutionary perspective of the CRISPR-Cas and anti-CRISPR protein systems. CRISPRimmunity also provide a web server, offering a well-designed graphical interface, a detailed tutorial, and multi-faceted information, making it easy to use and facilitating future experimental design and further data mining. The web server is accessible without registration at http://www.microbiome-bigdata.com/CRISPRimmunity/index/home.

<p id="Requirements"></p>

## Requirements

The source code is written by python3. In addition, several tools have been applied in CRISPRimmunity.  <br>

First, please check for the following packages and install them as needed:

1. Python packages：numpy，Biopython，sklearn，Levenshtein，re，lxml，requests，shutil，pandas，numpy，scipy，matplotlib，codecs，multiprocessing
2. R packages：treeio，ggtree，ggrepel

Second, please check for the following tools and install them as needed:

1. Prokka in https://github.com/tseemann/prokka<br>
2. diamond in https://github.com/bbuchfink/diamond<br>
3. BLAST+ in https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/<br>
4. HMMER in http://hmmer.org/download.html<br>
5. cd-hit in https://www.bioinformatics.org/cd-hit<br>
6. pilercr in https://www.drive5.com/pilercr/<br>
7. CRISPRCasFinder in https://github.com/dcouvin/CRISPRCasFinder<br>
8. CRISPRidentify in https://github.com/BackofenLab/CRISPRidentify<br>
9. AcRanker-master in https://github.com/amina01/AcRanker<br>
10. mafft in https://mafft.cbrc.jp/alignment/software/<br>
19. clustalo in http://www.clustal.org/omega/<br>
19. mview in https://github.com/desmid/mview<br>
19. FastTree in http://www.microbesonline.org/fasttree/<br>

<p id="Install"></p>

## Install

### Linux

- step1: Clone the repository of  CRISPRimmunity from Github

```shell
git clone https://github.com/HIT-ImmunologyLab/CRISPRimmunity.git
```

- step2: Download CRISPRimmunity database and data for standalone from CRISPRimmunity webserver and check MD5

Due to the large size of the database file (~XXXG), we recommend using -c (continue) and -b (background) parameters of wget to avoid the data loss caused by network outages, and checking MD5 to verify the integrity of the downloaded files.

```shell
## md5 checksum: XXX
### Download Database
wget -c -b XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
### Check md5 (e2ac0981ab07b4529d857c5c778f30e6 )
md5sum database.tar.gz
### Download Data 
wget -c -b XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
### Check md5 (6f577bbc661276e597b237f80870f2b0)
md5sum data.tar.gz
```

- step3: Unzip the database and data file to specified subdirectory under CRISPRimmunity installation directory

```shell
### Unzip the database file
tar -zxvf database.tar.gz
cp -r <download_database_path>/database <CRISPRimmunity_path>/CRISPRimmunity
### Unzip the database file
tar -zxvf data.tar.gz
cp -r <download_data_path>/data <CRISPRimmunity_path>/CRISPRimmunity
```

When you organize the whole files well,the corresponding directory structure are displayed as shown below. <br>

![image](https://github.com/LYC2015000421/CRISPRimmunity/raw/main/picture/tree2.png)

- step4:  Create soft links to dependent tools in the software directory

```shell
### Create soft links
ln -s <path>/prokka <CRISPRimmunity_path>/software/prokka
ln -s <path>/diamond <CRISPRimmunity_path>/software/diamond
ln -s <path>/blastn <CRISPRimmunity_path>/software/blastn
ln -s <path>/blastdbcmd <CRISPRimmunity_path>/software/blastdbcmd
ln -s <path>/rpsblast <CRISPRimmunity_path>/software/rpsblast
ln -s <path>/prodigal <CRISPRimmunity_path>/software/prodigal
ln -s <path>/hmmscan <CRISPRimmunity_path>/software/hmmscan
ln -s <path>/nhmmscan <CRISPRimmunity_path>/software/nhmmscan
ln -s <path>/CRISPRCasFinder.pl <CRISPRimmunity_path>/software/CRISPRCasFinder.pl
ln -s <path>/CRT1.2-CLI.jar <CRISPRimmunity_path>/software/CRT1.2-CLI.jar
ln -s <path>/pilercr <CRISPRimmunity_path>/software/pilercr
ln -s <path>/clustalo <CRISPRimmunity_path>/software/clustalo
ln -s <path>/mafft <CRISPRimmunity_path>/software/mafft
ln -s <path>/mview <CRISPRimmunity_path>/software/mview
ln -s <path>/FastTree <CRISPRimmunity_path>/software/FastTree
ln -s <path>/acranker.py <CRISPRimmunity_path>/software/acranker.py
### Copy CRISPRidentify to software directory
cp -r <path>/CRISPRidentify <CRISPRimmunity_path>/software
```

<p id="Usage"></p>

## Usage

CRISPRimmunity is an comprehensive tool for Acr prediction, identification of novel class 2 CRISPR-Cas loci, and dissection of CRISPR-associated important molecular events.

### Command line options

```shell
General:
    --input <file name>        : Genome file path: FASTA or GenBank file
    --output <folder name>     : Output folder in which results will be stored
    --annotation <x>           : function selection: effector,anti,self-targeting,prophage,MGE-targeting
    --prog_dir <x>             : directory of CRISPRimmunity

Acr prediction module:
    --anti_flag <x>            : flag for predicting homologs for known or novel Acrs, values: novel, known
    --anti_identity <x>        : minimum identity to report an alignment with konwn Acrs, default:0.4 [0-1]
    --anti_coverage <x>        : minimum coverage to report an alignment with konwn Acrs, default:0.7 [0-1]
    --anti_up_num <x>          : maximum protein number separated from HTH domain-containing proteins, default:3
    --anti_prot_size <x>       : maximum protein size for Acr (aa), default:400

Important molecular events annotation module:
    Detection and classification of CRISPR-Cas event:
        --att_spacer_num <x>       : minimum spacer of the CRISPR array, default:2
        --att_min_repeat_length <x>: minimum repeat length  of the CRISPR array, default:18
        --att_max_repeat_length <x>: maximum repeat length  of the CRISPR array, default:45
        --range <x>                : length of the sequence of interest in the vicinity of the CRISPR array (bp), default:20000
        --neighbor <x>             : number of the protein of interest in the vicinity of the CRISPR array, default:10
        --repeat_identity <x>      : minimum identity to report an alignment with konwn type repeats, default:0.9 [0-1]
        --repeat_coverage <x>      : minimum coverage to report an alignment with konwn type repeats, default:0.9 [0-1]
	
    Detection of self-targeting event:
    --mismatch <x>             : maximum mismatch to report an alignment with self-genome for self-targeting , default:2
    --cov <x>                  : minimum coverage to report an alignment with self-genome for self-targeting, default:1
    
    Detection of prophage event:
    --prophage_thread_num <x>  : number of threads (CPUs) to use in the prophage detection, default:3
    
    Detection of self-targeting event:
    --mismatch <x>             : maximum mismatch to report an alignment with self-genome for self-targeting , default:2
    --cov <x>                  : minimum coverage to report an alignment with self-genome for self-targeting, default:1
    
    Detection of MGE-targeting event:
    --phage_mismatch <x>       : maximum mismatch to report an alignment with MGE database for MGE-targeting, default:2
	--phage_coverage <x>       : minimum coverage to report an alignment with MGE database for MGE-targeting, default:0.7

Novel class 2 CRISPR-Cas loci identification module:
    --att_protein_size <x>     : minimum protein size for novel class 2 effector candidate, default:500
    --att_neighbor_distance <x>: maximum protein number of interest in the vicinity of the CRISPR array for mining novel class 2 effector candidate, default:5
    --homo_identity <x>        : minimum identity to report an alignment with novel class 2 effector candidate, , default:0.3
    --homo_coverage <x>        : minimum coverage to report an alignment with novel class 2 effector candidate, default:0.6
```

### Start 

1. Identification of novel class 2 effector protein:

```shell
python <path>/predict_crispr_cas_self_targeting.py --input <genome file path> --output <outdir> --annotation 'effector' --prog_dir <program directory>
```

2. Acr prediction:

```shell
python <path>/predict_crispr_cas_self_targeting.py --input <genome file path> --output <outdir> --annotation 'anti,self-targeting,prophage' --anti_flag novel --prog_dir <program directory>
```

3. Annotation of important molecular events:

```shell
python <path>/predict_crispr_cas_self_targeting.py --input <genome file path> --output <outdir> --annotation 'self-targeting,MGE-targeting,prophage,anti' --prog_dir <program directory>
```

### Outputs

| File Name                                                    | Description                                                  | Module        |
| :----------------------------------------------------------- | ------------------------------------------------------------ | ------------- |
| crispr_cas_self-targeting_statis_result.txt                  | Overview of the result  of CRISPR-associated event, including CRISPR array detection, Cas protein annotation, repeat type annotation, self-targeting detection. | IME, ACR, EFF |
| anti_statis_result.txt                                       | Overview of the result of predicted Acrs.                    | IME, ACR      |
| <prefix>.merge_spc_blastn_fasta_add_<br>cov_filter_mismatch_2_cov_1_hit_not_in_array | Detail of the result of sequence alignment between spacer and self-genome. | IME, ACR, EFF |
| merge_crispr_array.txt<prefix>\_merge.spc                    | Overview  of the result of CRISPR array detection.           | IME, ACR, EFF |
| merge_crispr_array.txt                                       | Detail of the result of CRISPR array detection.              | IME, ACR, EFF |
| <prefix>\_merge.csp                                          | Sequence of the spacers.                                     | IME, ACR, EFF |
| self-targeting/spacer_target_hit_result.txt                  | Detail of the result of self-targeting detection.            | IME, ACR      |
| self-targeting/\*.\*\_multialign                             | Detail of the result of spacer-protospacer alignment.        | IME, ACR      |
| self-targeting/\*.\*\_fa                                     | Detail of the result of spacer-protospacer sequence.         | IME, ACR      |
| prophage/complete_prophage_summary.txt                       | Overview  of the result of prophage detection.               | IME, ACR      |
| prophage/complete_prophage_detail.txt                        | Detail of the result of prophage region.                     | IME, ACR      |
| prophage/complete_prophage_protein.faa                       | Protein sequence of prophage region.                         | IME, ACR      |
| prophage/complete_prophage_nucl.fa                           | Nucleotide sequence of prophage region.                      | IME, ACR      |
| interacting_phage/spacer_blastn_phage_filter.txt             | Overview  of the result of MGE-targeting detection.          | IME, ACR      |
| anti/acr_candidate_info_integrated_three_methods.txt         | Overview of the predicted ant.                               | IME, ACR      |
| novel_effector_protein/filter_repeat_homo_<br>protein_crispr_filter_result.txt | Overview of the identified novel class 2 effector proteins   | EFF           |
| novel_effector_protein/effector_tree.png                     | Phylogenetic tree of novel class 2 effector candidates and known class 2 effector proteins | EFF           |

IME：important molecular events detection module
ACR：Acr prediction module
EFF：novel class 2 effector proteins identification module

<p id="Example"></p>

## Example

In the CRISPRimmunity package, there is a directory named "example" which includes input sequence and result for examples of three modules.  "Armatimonadetes-bacterium-isolate-ATM2-J3.gb"  is used for identifying novel class 2 effector proteins.  "Staphylococcus-schleiferi-strain-5909-02.gb" is used for Acr prediction and annotation of important molecular events, respectively.

```python
## run the example
python example_run.py
```

(1) Result for identifying novel class 2 effector proteins

![image](https://github.com/LYC2015000421/CRISPRimmunity/raw/main/picture/eff.png)

(2) Result for Acr prediction and annotation of important molecular events

![image](https://github.com/LYC2015000421/CRISPRimmunity/raw/main/picture/anti.png)

<p id="Contributors"></p>

## Contributors

This project exists thanks to all the people who contribute.

<p id="License"></p>

## License

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.



