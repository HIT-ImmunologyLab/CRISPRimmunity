# CRISPRimmunity: CRISPR and Immunity



## Table of Contents

- [Background](#Background)

- [Requirements](#Requirements)
- [Install](#Install)
- [Usage](#Usage)
- [Visualizations](#Visualizations)
- [Contributors](#Contributors)
- [License](#License)

<p id="Background"></p>

## Background

The CRISPR-Cas system is a highly adaptive and RNA-guided immune system found in bacteria and archaea that has applications as a genome editing tool and is a valuable system for studying co-evolutionary dynamics of bacteriophage interactions. This repository presents CRISPRimmunity, a new tool designed for Acr prediction, identification of novel class 2 CRISPR-Cas loci, and dissection of CRISPR-associated important molecular events. CRISPRimmunity is built on a series of CRISPR-oriented databases that offer a comprehensive co-evolutionary perspective of the CRISPR-Cas and anti-CRISPR protein systems. CRISPRimmunity also provide a web server, offering a well-designed graphical interface, a detailed tutorial, and multi-faceted information, making it easy to use and facilitating future experimental design and further data mining. The web server is accessible without registration at http://www.microbiome-bigdata.com/CRISPRimmunity/index/home.

<p id="Requirements"></p>

## Requirements

The source code is written by python3. In addition, several tools have been applied in DBSCAN-SWA. Among these, Prokka requires installtion by users. <br>

First, please install the following python packages:

1. numpy

2. Biopython

3. sklearn

Second, please install the following tools:

1. Prokka in https://github.com/tseemann/prokka<br>

```shell
git clone https://github.com/tseemann/prokka.git
# install the dependencies:
sudo apt-get -y install bioperl libdatetime-perl libxml-simple-perl libdigest-md5-perl
# install perl package
sudo bash
export PERL_MM_USE_DEFAULT=1
export PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"
perl -MCPAN -e 'install "XML::Simple"'
# install the prokka databases
prokka --setupdb
# test the installed prokka databases
prokka --listdb
```

**warning**: Prokka needs blast+ 2.8 or higher, so we provide the version of blast+ in bin directory, the users can install a latest blast+ and add it to the environment or use the blast+ provided by DBSCAN-SWA. Please ensure the usage of blast+ in your environment by eg: 

```shell
which makeblastdb
```

2. FragGeneScan in https://github.com/gaberoo/FragGeneScan<br>
3. diamond in https://github.com/bbuchfink/diamond<br>
4. BLAST+ in https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/<br>
6. AcRanker-master in https://github.com/amina01/AcRanker<br>
7. mafft in https://mafft.cbrc.jp/alignment/software/
8. HMMER in http://hmmer.org/download.html<br>
9. cd-hit in https://www.bioinformatics.org/cd-hit<br>
11. SeqKit in https://bioinf.shenwei.me/seqkit<br>
12. DBSCAN-SWA in https://github.com/HIT-ImmunologyLab/DBSCAN-SWA<br>
13. pilercr in https://www.drive5.com/pilercr/<br>
14. CRISPRCasFinder in https://github.com/dcouvin/CRISPRCasFinder<br>
15. Self-Targeting-Spacer-Searcher in https://github.com/kew222/Self-Targeting-Spacer-Searcher<br>
17. PHIS in http://www.phis.inra.fr/<br>
18. Phage_finder in http://phage-finder.sourceforge.net/<br>
19. CRISPRLeader in https://github.com/BackofenLab/CRISPRleader<br>

<p id="Install"></p>

## Install

### Linux

- step1: Clone the repository of  CRISPRimmunity from Github

```shell
git clone https://github.com/HIT-ImmunologyLab/CRISPRimmunity.git
```

- step2: Download CRISPRimmunity database for standalone from CRISPRimmunity webserver and check MD5

Due to the large size of the database file (~XXXG), we recommend using -c (continue) and -b (background) parameters of wget to avoid the data loss caused by network outages, and checking MD5 to verify the integrity of the downloaded files.

```shell
## md5 checksum: XXX
### Download Database
wget -c -b XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
### Check md5
md5sum XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
```

- step3: Unzip the database file to specified subdirectory under CRISPRimmunity installation directory

```shell
### Unzip the database file
tar -zxvf XXXXXXXXXXXXXXXXXXXXXXXXXX
cp XXXXXXXXXXXXXXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXX
```

When you organize the whole files well,the corresponding directory structure are displayed as shown below. 

<p id="Usage"></p>

## Usage

CRISPRimmunity is an integrated tool for Acr prediction, identification of novel class 2 CRISPR-Cas loci, and dissection of CRISPR-associated important molecular events.

### Command line options

```shell
General:
--input <file name>        : Genome file path: FASTA or GenBank file
--output <folder name>     : Output folder in which results will be stored
--mismatch <x>             : value of self-targeting mismatch, default:2
--cov <x>                  : value of self-targeting coverage, default:1
--range <x>                : array neighbour range, default:20000
--neighbor <x>             : array neighbour number, default:10
--anti_identity <x>        : array neighbour number, default:0.4
--anti_coverage <x>        : array neighbour number, default:0.7
--anti_up_num <x>          : array neighbour number, default:3
--phage_mismatch <x>       : array neighbour number, default:2
--phage_coverage <x>       : array neighbour number, default:0.7
--prophage_thread_num <x>  : array neighbour number, default:3
--anti_prot_size <x>       : array neighbour number, default:3
--repeat_identity <x>      : short-blastn identity, default:0.9
--repeat_coverage <x>      : short-blastn coverage, default:0.9
--homo_identity <x>        : diamond blastp identity, default:0.3
--homo_coverage <x>        : diamond blastp coverage, default:0.6
--att_spacer_num <x>       : diamond blastp coverage, default:2
--att_min_repeat_length <x>: diamond blastp coverage, default:18
--att_max_repeat_length <x>: diamond blastp coverage, default:45
--att_protein_size <x>     : diamond blastp coverage, default:500
--att_neighbor_distance <x>: diamond blastp coverage, default:5
--annotation <x>           : diamond blastp coverage
--anti_flag <x>            : diamond blastp coverage
--prog_dir <x>             : directory of program
```

### Start 

The python script is also provided for expert users<br>

1. identification of novel class 2 effector protein:

```shell
python <path>/predict_crispr_cas_self_targeting.py --input <genome file path> --output <outdir> --annotation 'effector' --prog_dir <program directory>
```

2. Acr prediction:

```shell
python <path>/predict_crispr_cas_self_targeting.py --input <genome file path> --output <outdir> --annotation 'anti,self-targeting,prophage' -- prog_dir <program directory>
```

3. annotation of important molecular events:

```shell
python <path>/predict_crispr_cas_self_targeting.py --input <genome file path> --output <outdir> --annotation 'self-targeting,MGE-targeting,prophage,anti' -- prog_dir <program directory>
```

### Outputs

| File Name                                     | Description                                                  |
| --------------------------------------------- | ------------------------------------------------------------ |
| \<prefix\>\_DBSCAN-SWA\_prophage\_summary.txt | the tab-delimited table contains predicted prophage informations including prophage location, specific phage-related key words, CDS number, infecting virus species by a majority vote and att sites |
| \<prefix\>\_DBSCAN-SWA\_prophage.txt          | this table not only contains the information in <prefix>\_DBSCAN-SWA\_prophage\_summary.txt but also contains the detailed information of prophage proteins and hit parameters between the prophage protein and hit uniprot virus protein |
| <prefix>\_DBSCAN-SWA\_prophage.fna            | all predicted prophage Nucleotide sequences in FASTA format  |
| <prefix>\_DBSCAN-SWA\_prophage.faa            | all predicted prophage protein sequences in FASTA format     |
| **Phage Annotation**                          | if add\_annotation!=none, the following files are in "prophage\_annotation" |
| <prefix>\_prophage\_annotate\_phage.txt       | the tab-delimited table contains the information of prophage-phage pairs with prophage\_homolog\_percent, prophage\_alignment\_identity and prophage\_alignment\_coverage |
| <prefix>\_prophage\_annotate\_phage.txt       | the table contains the detailed information of bacterium-phage interactions including blastp and blastn results |

<p id="Visualizations"></p>

## Visualizations

You can find a directory named "test" in the DBSCAN-SWA package. The following visualzations come from the predicted results of Xylella fastidiosa Temecula1(NC\_004556)<br>
(1) the genome viewer to display all predicted prophages and att sites

![image](https://github.com/LYC2015000421/work-progress/raw/master/Pictures/AFDB1.png)

(2) the detailed information of predicted prophages

![image](https://github.com/LYC2015000421/work-progress/raw/master/Pictures/AFDB2.png)

(3) If the users set add_annotation as PGPD or the phage file path, the detailed information of bacterium-phage interaction will be illustrated as follows:

![image](https://github.com/LYC2015000421/work-progress/raw/master/Pictures/Pythonpack.png)

<p id="Contributors"></p>

## Contributors

This project exists thanks to all the people who contribute.

<p id="License"></p>

## License

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.



