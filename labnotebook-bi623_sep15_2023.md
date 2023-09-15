---
title: "QAA Final Report"
author: "Ike Sanderson"
date: "2023-09-14"
output:
  pdf_document: default
  word_document: default
  html_document:
    df_print: paged
toc: yes
extra_dependencies: float
---
# labnotebook-bi623.md

## In Class Activity #1 Aug 28, 2023 in Bi623

Go to: 
https://www.ncbi.nlm.nih.gov/search/
Consult query manual:
https://www.ncbi.nlm.nih.gov/books/NBK3837/

Type myeloperoxidase in the search box
1. How many hits are there for the Nucleotide, Gene, and GEO DataSets databases?
```
Nucleotide: 3591 hits

Gene: 1101 hits

Geo Datasets: 654

```
2. What types of data are in each of these 3 databases?
```
Nucleotide: Full sequences. Apart from sequence data in the EST (Expressed Sequence Tag) and GSS (Genome Survey Sequence divisions of GenBank, the Nucleotide database contains all the sequence data from GenBank, EMBL, and DDBJ, the members of the International Nucleotide Sequence Databases Collaboration (INSDC). Nucleotide also includes NCBI-curated Reference Sequences (RefSeqs), submitted assemblies and annotations from the Third Party Annotation (TPA) database, and nucleotide sequences extracted from structure records from the Protein Databank (PDB).

Gene: Gene description. Gene is a searchable database of genes, focusing on genomes that have been completely sequenced and that have an active research community to contribute gene-specific data. Information in Gene records includes nomenclature, chromosomal localization, gene products and their attributes (e.g., protein interactions), associated markers, phenotypes, interactions, and links to citations, sequences, variation details, maps, expression reports, homologs, protein domain content, and external databases.

Geo Datasets: Gene expression omnibus. Fully analyzed databases. GEO Datasets stores curated gene expression and molecular abundance data sets assembled by NCBI from the Gene Expression Omnibus (GEO) repository of microarray data.

```

3. Using information from NCBI (not Wikipedia, or another source), can you find out what the gene does? (hint: think about NCBI’s literature database)


Followed the Gene database link and the first link
Summary tells me it is an enzyme with hypohalous acids central to the microbicidal activity of neutrophils. Kills microbes.

4. By looking at a few entries, what is the gene’s “symbol”?
MPO

Do a second search using just the symbol
5. Did the number of hits in the databases above change? If so, why might this be?
Fewer hits- perhaps we are narrowing our database search. It's a narrowed field- down to one; gene symbols are not consistent across organisms

Narrow your original search using a Boolean operator with an indexed field of your choice
6. How did this affect your results?
Different results. Example: MPO AND NOT human
or
**************** for quiz!!!! Indexed fields: [organism] NOT human
[organisms zebrafish]

Use an indexed field search to find the gene in zebrafish (use the scientific name)
7. Is the gene symbol consistent across organisms?
No

8. Find the entry for the protein encoded by the gene.
******************** for quiz!!!!! play in manual for things abvout proteins

9. Are there redundant entries for the protein, and if so how many in total?
yes! 
[whatever search is] AND protein


Next: SRA at NIH
https://www.ncbi.nlm.nih.gov/sra/
See slide 51 in Aug 28 Lecture Bi 621 (Lecture 1) 

I WILL USE THIS LATER   

ICA 1
log on to talapas

My fastq-tools version was 2.9.3 and gave me errors when I was 
create new environment
$ conda create --name bgmp-sra

2. Using your Unix skills, count the number of sequences per file. Also, take a few longer sequences (10 or so) and BLAST them (you can use the NCBI Blast web tool this time; GOOGLE TIME!). Can you tell what type of data these are?

cat and head to read 
wc -l to count lines and divide by 4 to count reads

$ cat SRR849828.fastq | wc -l
520136 / 4 =
130034 reads

$ cat SRR849968.fastq | wc -l
228852 / 4 =
57213 reads

Looked for long reads with awk:
cat SRR849828.fastq | awk 'length($0) > 600'

result still included quality scores and headers, but I could select the reads. 

Made a VSCode file and pasted in reads, created a header. Example below:
>SRA1SRR849828
CGGCCGACATGTTTTGTTTTTTTTTCTTTTTTTTTCAGTATGATACACACAGTTCAGTTGTATTCAACTTTTAAGTTATTTTAAGTTATTTAATCTTTTAAAACGGTTTCTTTTCAAACTTTTATTAGGGAACAATTTTGATATAACTTAATTCAATAATTTTAATTGAATCATTAAGAAAGTTTGGGGGGGGGGTTCATTTTCCTATCTTTAAGTACAGTTAACACATGGGCCTTCCGTGCCGTTATATTTTGATGTTGGTCATTAAAATGGTATTGGTTAGGCAAGTATTTTTGGAGGTCGTACTTGGTGTAAAAACCTTGAGAATCAATGAACTAGTCACAAAATCAGTTAGGTAAATTAGAATTGGGTTCTTGGACTTGCAGTGTAAATCCACCCTTGGAATGGAAGACGTGGAAATGCTTAAACCTTTATTGGGATTTTCACCCAAATTGAAAAGTTTCTTACTTGACACGTTCCGCCCTGCCACAAAGTTTCATGTTTTGTAACAGTGATCAAAACAATCCGCCTAGTCTCTACACGTTCACGCTTGGAGATGCTCTGTACAGAATGTGTTGTGTGTATGTAGATAATGTGTAGCG

SAved file with .fa and no other suffix

Uploaded to 
https://blast.ncbi.nlm.nih.gov/Blast.cgi

Predicted as
Syngnathus scovelli NADH:ubiquinone oxidoreductase complex assembly factor 8 (ndufaf8), mRNA

Accesion ID:
Syngnathus floridae non-pregnant male brood pouch (SRR849968)

I checked the accession ID:
https://www.ncbi.nlm.nih.gov/sra/SRR849968


https://www.ncbi.nlm.nih.gov/sra/?term=SRR849828
Syngnathus scovelli pregnant male brood pouch (SRR849828)

But when I checked the reads in the BLAST tool, the second accession reads mapped out as:
Syngnathus acus

Because these reads have long stretches of AAAAAAAAAAs, these are mRNA reads. RNA seq.

Part 2: Ensembl

1. Explore basic information on Ensembl for the following genomes: yeast, C. elegans, green spotted puffer (aka Tetraodon), and human.

2. Download the protein sequence .fasta files for each.
3. Do different isoforms appear to be included in each file?

4. Using your Unix skills (and wandering through Ensembl assembly statistics), figure out what the density of annotated coding genes (number of genes per genome assembly megabase) is for each genome.

5. Write a simple Unix one-liner to print a list of protein lengths (including all isoforms when present), and apply it to each of the 4 .fasta.gz protein files.




##############################################################################

## Bi623 Assignment 1: RBH

I decided to use the files uploaded by Leslie just because I know those have solid data. Copied them
 to my own directory.

Ran a bash sort command first, because it took away a required check in the python code:
```
sort -g -k1 -k11 H_to_zfishdb.blastp > H_to_zfishSORTE.blastp
sort -g -k1 -k11 Z_to_homodb.blastp > Z_to_homodbSORTE.blastp
```
I wouldn't need to check that the evalue was the lowest. First hit will be lowest.

I know I need to look at each line of the Z to H protein blast & H to Z protein blast,
and store both proteinID. I know I didn't need to store evalue but I did anyway

I had some coding decisions to make about how to handle some things:


Had an error resulting from proteins in zfish BH not existing in human BH.

Had to add line to write statement to ignore (i.e. continue) any protein 
  that didn't exist
  
Still getting too many RBHs compared to everyone else:

(base) [ikes@login2 Bi623]$ ./syn.py 
length of protdictH is 20849
length of protdictZ is 27927
length of biomart H dictionary is 121767
length of biomart Z dictionary is 52090
bannedH="5803
bannedZ="6082
(base) [ikes@login2 Bi623]$ wc -l Human_Zebrafish_RBH.tsv 
20826 Human_Zebrafish_RBH.tsv


Realized I wasn't deleting the banned protein hits. THey still sat in my dictionaries. 
SO I tried writing code that would delete them as I identified them within the loop.
 That, predictably, did not go well. Python does not like it when I change the 
 length of a dictionary while I am iterating through it.
Thus I had to resort to making a separate list to keep track of the banned
 proteins. THen go through the list and delete the proteins from the dictionary.
 
```
if not protdictZ.get(proteinidZ) and proteinidZ not in banlist:
  protdictZ[proteinidZ] = [proteinidH,evalueZ] 
elif proteinidZ in protdictZ and evalueZ == protdictZ[proteinidZ][1]:
  banlist.append(proteinidZ) 
  del protdictZ[proteinidZ] 
else:
  continue
```
 
Made a final dictionary from which I will write. This was a change. 
First version of code was a complex set of variables but somehow I kept getting
 too many lines in final file. 
 
Eventually got code to work and got correct number of genes:

(base) [ikes@login4 Bi623]$ ./syn.py 
length of protdictH is 15046
length of protdictZ is 21845
length of biomart H dictionary is 121767
length of biomart Z dictionary is 52090
Length of result dict: 7975
(base) [ikes@login4 Bi623]$ wc -l Human_Zebrafish_RBH.tsv 
7975 Human_Zebrafish_RBH.tsv

Submitted assignment on time.


##############################################################################

## Assignment 2: Dot Plot
01. Using Ensembl’s Biomart, extract Ensembl GeneIDs, and chromosomal positions for all genes of human and zebrafish (separately).
Supposed to use 109, which means I need to use ensembl archives. Google search turns up:
https://useast.ensembl.org/info/website/archives/index.html

On ensembl: 
select the attributes: gene stable ID, gene name, protein stable ID, Chromosome/scaffold name, 
   gene start (bp), gene end (bp)
and press Return


notes to self about graphing:
search about graphing gg plot point (that's the format I'll use)
to put axis labels in order, look into factors and levels to give me an ordered vector that I'll put into the 
  gg plot argument as the X and Y axes
ensemble will put tegther the file and begin download

02. read in files

```
RBHtable = read.table("/Users/ikesanderson/bioinfo/Bi623/dot-plot/Human_Zebrafish_RBH.tsv", header = FALSE, sep = "\t") 
zgenes = read.csv("/Users/ikesanderson/bioinfo/Bi623/dot-plot/zfish_genes.txt", header = TRUE, sep = "\t")
hgenes = read.csv("/Users/ikesanderson/bioinfo/Bi623/dot-plot/human_genes.txt", header = TRUE, sep = "\t")
```

reason for header FALSE in RBH is that I didn't include headers in that file. So I need to add column headers. 
I think I need to use rename() for that

Also, when I pipe to head and take a look at the other two files, I see that their headers have . in them

Gene.stable.ID
Protein.stable.ID
Gene.start..bp.
Gene.end..bp.
Gene.name
Chromosome.scaffold.name

Gene.stable.ID
Gene.start..bp.
Gene.end..bp.
Gene.name
Protein.stable.ID
Chromosome.scaffold.name

Changed the RBH columns to:

Human.Gene.stable.ID
Human.Protein.stable.ID
Human.Gene.name
ZFish.Gene.Stable.ID
ZFish.Protein.stable.ID
ZFish.Gene.name

using:
```
RBHtable %>% 
  rename(., Human.Gene.stable.ID = V1, Human.Protein.stable.ID = V2, Human.Gene.name = V3, ZFish.Gene.Stable.ID = V4, ZFish.Protein.stable.ID = V5,ZFish.Gene.name = V6)
```
Also renamed human and zebra fish columns using:
```
hgenes %>% 
  rename(., Human.Gene.stable.ID = Gene.stable.ID, Human.Gene.start.bp = Gene.start..bp., Human.Gene.end.bp = Gene.end..bp., Human.Gene.name = Gene.name, Human.Protein.stable.ID = Protein.stable.ID, Human.Chromosome = Chromosome.scaffold.name)
```

Fri Sep 8, 2023
Was able to use left join but I needed to un-rename columns to prevent duplication of columns.
The only column that can have the same name is the one I join by
Also, my first left_join() command attempt used keep = TRUE in the arguments,
but that's not what I want.keep = FALSE

#need to figure out how to filter out the scaffolds and patches



#part 2

Sat Sep 9, 2023

Pick a RBH pair:
H bp start  H bp start  HGene Name  H Protein       H Chr
31488688	  31618588	  SFI1	      ENSP00000416931	22

ZFish Protein       Z bp st   Z bp end  Zgn name  ZFish Chr
ENSDARP00000156723	40898877	40922971	sfi1	    6

Here is the Ensembl page about how to make the tree:
http://useast.ensembl.org/Help/View?id=137


##############################################################################


## Assignment: QAA (Leslie)
#Date: Wed 9/6/2023#
001. followed instructions and created new git repo "QAA"
note: git makes new dir nested inside whatever dir I pick for repo. 
in future: do not premake a repo dir with same name as repo

002. login to talapas using iTerm:
❯ ssh ikes@login.talapas.uoregon.edu

003.followed link in assignment description to:
/projects/bgmp/shared/Bi623/QAA_data_assignments.txt

and see that I have 
Ike	24_4A_control_S18_L008	1_2A_control_S1_L008
23_4A_control_S17_L008_R1_001.fastq.gz
23_4A_control_S17_L008_R2_001.fastq.gz
1_2A_control_S1_L008_R1_001.fastq.gz
1_2A_control_S1_L008_R2_001.fastq.gz

So I will need to look at four files: R1 and R2 for reach of those two reads

004. while in login node on talapas I used:
$ module spider fastqc

that returned:
This module can be loaded directly: module load fastqc/0.11.5

So I entered:
$ module load fastqc/0.11.5
$ ml

Currently Loaded Modules:
  1) gcc/7.3   2) fastqc/0.11.5
  
Got help/documentation:
$ fastqc -h


SYNOPSIS

	fastqc seqfile1 seqfile2 .. seqfileN

    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
           [-c contaminant file] seqfile1 .. seqfileN
           
wonder if I will need to use the lists of adapters in the index list

might need this:
-o --outdir     Create all output files in the specified output directory.
                    Please note that this directory must exist as the program
                    will not create it.  If this option is not set then the
                    output file for each sequence file is created in the same
                    directory as the sequence file which was processed.
                    
                    
-a              Specifies a non-default file which contains the list of
--adapters      adapter sequences which will be explicity searched against
                    the library. The file must contain sets of named adapters
                    in the form name[tab]sequence.  Lines prefixed with a hash
                    will be ignored.


#Date: Thu Sep 7, 2023#
fastqc --outdir /projects/bgmp/ikes/bioinfo/Bi623/QAA /projects/bgmp/shared/2017_sequencing/demultiplexed/23_4A_control_S17_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/23_4A_control_S17_L008_R2_001.fastq.gz

and then 

fastqc --outdir /projects/bgmp/ikes/bioinfo/Bi623/QAA /projects/bgmp/shared/2017_sequencing/demultiplexed/1_2A_control_S1_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/1_2A_control_S1_L008_R2_001.fastq.gz

005. Made clone of git repo on my computer by navigating (using iTerm) to the 
  directory I have for Bi623 and 

```
   >mkdir QAA
```   
and then I cloned using

```
❯ git clone https://github.com/Ike-Sanderson/QAA.git
```
note: the cloned repo makes a new folder inside my QAA. so it's a nested QAA
future: don't premake a folder for the repo

006.I copied (via SCP) the files from talapas down to the cloned repo 

scp ikes@login.talapas.uoregon.edu:<source path> <destination path, which= . if iTerm pwd>

scp ikes@login.talapas.uoregon.edu:/projects/bgmp/ikes/bioinfo/Bi623/QAA/*.* .

unzipped the four zipped archives (via gui)

007. go to interactive mode on talapas:

could not remember how to do srun command, so used 
$srun --account=bgmp --partition=compute --time=2:00:00 --pty bash
$./qualscoretohisto.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/23_4A_control_S17_L008_R1_001.fastq.gz -o 23_4A_control_S17_L008_R1_001.png -l 101

opened second terminal, 
$srun --account=bgmp --partition=compute --time=2:00:00 --pty bash
$./qualscoretohisto.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/23_4A_control_S17_L008_R2_001.fastq.gz -o 23_4A_control_S17_L008_R2_001.png -l 101

note: argparse future: make base outfile from user, but script the extenstion (eg .png)
note: also write code to figure out length so the third argparse not needed
Next time: write an sbatch wrapper for the other two jobs
tried the wrapper for the other two files but ran into errors in my script. Went back to running it in interactive mode.

Run time seemed really long, but I don't have comparison data.
Copied all four .png files to the repo on my computer.

scp ikes@login.talapas.uoregon.edu:/projects/bgmp/ikes/bioinfo/Bi623/QAA/*.png .

Note for FastQC analysis: 
GC content for mouse avg 51.3% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2909565/)


008. Made new environment

$ conda create --name bgmp_qaa
$ conda activate bgmp_qaa

installed cutadapt and trimmomatic

$conda install <app name>

verified version

$trimmomatic -version
cutadapt version 4.4
trimmomatic version 0.39

###next Part 2 #5
Date Fri Sep 8, 2023
-read documentation for cutadapt and for trimmomatic

cutadapt -a -A

cutadapt removes adapter sequences from high-throughput sequencing reads.
Usage:
    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq
For paired-end reads:
    cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq

Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
characters are supported. All reads from input.fastq will be written to
output.fastq with the adapter sequence removed. Adapter matching is
error-tolerant. Multiple adapter sequences can be given (use further -a
options), but only the best-matching adapter will be removed.

Input may also be in FASTA format. Compressed input and output is supported and
auto-detected from the file name (.gz, .xz, .bz2). Use the file name '-' for
standard input/output. Without the -o option, output is sent to standard output.

How to figure out universal Illumina adapter:
https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

Confirmed adapter sequence:
$ zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/1_2A_control_S1_L008_R1_001.fastq.gz | grep 'AGATCGGAAGAG'

For 1:
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o 1_2A_control_S1_L008_R1_001cut.fastq -p 1_2A_control_S1_L008_R2_001cut.fastq /projects/bgmp/shared/2017_sequencing/demultiplexed/1_2A_control_S1_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/1_2A_control_S1_L008_R2_001.fastq.gz

For 23:
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o 23_4A_control_S17_L008_R1_001cut.fastq -p 23_4A_control_S17_L008_R2_001cut.fastq /projects/bgmp/shared/2017_sequencing/demultiplexed/23_4A_control_S17_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/23_4A_control_S17_L008_R2_001.fastq.gz

for 1_2A
=== Summary ===

Total read pairs processed:         44,303,262
  Read 1 with adapter:               1,359,563 (3.1%)
  Read 2 with adapter:               1,657,134 (3.7%)
Pairs written (passing filters):    44,303,262 (100.0%)

Total basepairs processed: 8,949,258,924 bp
  Read 1: 4,474,629,462 bp
  Read 2: 4,474,629,462 bp
Total written (filtered):  8,925,098,135 bp (99.7%)
  Read 1: 4,463,208,431 bp
  Read 2: 4,461,889,704 bp

=== First read: Adapter 1 ===

Sequence: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA; Type: regular 3'; Length: 33; Trimmed: 1359563 times

Minimum overlap: 3
No. of allowed errors:
1-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-33 bp: 3

Bases preceding removed adapters:
  A: 23.4%
  C: 29.6%
  G: 30.7%
  T: 15.0%
  none/other: 1.3%
=== Second read: Adapter 2 ===

Sequence: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT; Type: regular 3'; Length: 33; Trimmed: 1657134 times

Minimum overlap: 3
No. of allowed errors:
1-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-33 bp: 3

Bases preceding removed adapters:
  A: 26.3%
  C: 28.5%
  G: 33.0%
  T: 11.1%
  none/other: 1.1%


It looks like I only saved the output from 1_2A trimmomatic?? how can I have been this forgetful?
I looked up on the web to see if there's a way to recover a standard out log and apparantly no.
In the future I should log in with this script:
```
~/terminal_logs/$(date +%Y%m%d-%H%M%S)-$(tty)-$$.log
```
which is supposed to save the session in a file snamed wit the date and time. I need to check this with 
the profs first.

Date: 9/11/2023
Trimmomatic manual:  
  http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
  and config the trim
  
  Instructions: 
  LEADING: quality of 3; 
  TRAILING: quality of 3; 
  SLIDING WINDOW: window size of 5 and required quality of 15; 
  MINLENGTH: 35 bases
  
SLIDINGWINDOW
Perform a sliding window trimming, cutting once the average quality within the window falls
below a threshold. By considering multiple bases, a single poor quality base will not cause the
removal of high quality data later in the read.
SLIDINGWINDOW:<windowSize>:<requiredQuality>
windowSize: specifies the number of bases to average across
requiredQuality: specifies the average quality required.

LEADING
Remove low quality bases from the beginning. As long as a base has a value below this
threshold the base is removed and the next base will be investigated.
LEADING:<quality>
quality: Specifies the minimum quality required to keep a base.

TRAILING
Remove low quality bases from the end. As long as a base has a value below this threshold
the base is removed and the next base (which as trimmomatic is starting from the 3‟ prime end
would be base preceding the just removed base) will be investigated. This approach can be
used removing the special illumina „low quality segment‟ regions (which are marked with
quality score of 2), but we recommend Sliding Window or MaxInfo instead
TRAILING:<quality>
quality: Specifies the minimum quality required to keep a base.

MINLEN
This module removes reads that fall below the specified minimal length. If required, it should
normally be after all other processing steps. Reads removed by this step will be counted and
included in the „dropped reads‟ count presented in the trimmomatic summary.
MINLEN:<length>
length: Specifies the minimum length of reads to be kept.

init interactive mode:
$srun --account=bgmp --partition=compute --time=2:00:00 --pty bash

Command:
trimmomatic PE <2 input filenames go here> <4 output files go here> LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35

trimmomatic PE 23_4A_control_S17_L008_R1_001cut.fastq 23_4A_control_S17_L008_R2_001cut.fastq trimcontrol/23_4A_forward_paired.fq.gz trimcontrol/23_4A_forward_unpaired.fq.gz trimcontrol/23_4A_reverse_paired.fq.gz trimcontrol/23_4A_reverse_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35
Quality encoding detected as phred33
Input Read Pairs: 44303262 Both Surviving: 42056576 (94.93%) Forward Only Surviving: 2176754 (4.91%) Reverse Only Surviving: 31679 (0.07%) Dropped: 38253 (0.09%)
TrimmomaticPE: Completed successfully
 
trimmomatic PE 1_2A_control_S1_L008_R1_001cut.fastq 1_2A_control_S1_L008_R2_001cut.fastq trimcontrol/1_2A_forward_paired.fq.gz trimcontrol/1_2A_forward_unpaired.fq.gz trimcontrol/1_2A_reverse_paired.fq.gz trimcontrol/1_2A_reverse_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35
Message:
Quality encoding detected as phred33
Input Read Pairs: 8477859 Both Surviving: 7966343 (93.97%) Forward Only Surviving: 500642 (5.91%) Reverse Only Surviving: 5482 (0.06%) Dropped: 5392 (0.06%)
TrimmomaticPE: Completed successfully

For #7 we are referred to ICA4 to plot read lengths:
"Plot the trimmed read length distributions for both R1 and R2 reads (on the same plot - yes, 
you will have to use Python or R to plot this, or Excel if you must. See ICA4 from Bi621). 
You can produce 2 different plots for your 2 different RNA-seq samples. There are a number 
of ways you could possibly do this. One useful thing your plot should show, for example, is 
whether R1s are trimmed more extensively than R2s, or vice versa. Comment on whether you 
expect R1s and R2s to be adapter-trimmed at different rates and why""

Here is the code I will use: (note: on Mac, the command is gzcat, but on Unix (Talapas it is zcat))
```
zcat /projects/bgmp/ikes/bioinfo/Bi623/QAA/trimcontrol/1_2A_forward_paired.fq.gz | grep -E "^[ATCGN]+$" | awk '{print length($0)}' | sort | uniq -c | sort -n > 1_2A_forward_pairedREADS.txt

zcat /projects/bgmp/ikes/bioinfo/Bi623/QAA/trimcontrol/1_2A_reverse_paired.fq.gz | grep -E "^[ATCGN]+$" | awk '{print length($0)}' | sort | uniq -c | sort -n > 1_2A_reverse_pairedREADS.txt

zcat /projects/bgmp/ikes/bioinfo/Bi623/QAA/trimcontrol/23_4A_forward_paired.fq.gz | grep -E "^[ATCGN]+$" | awk '{print length($0)}' | sort | uniq -c | sort -n > 23_4A_forward_pairedREADS.txt

zcat /projects/bgmp/ikes/bioinfo/Bi623/QAA/trimcontrol/23_4A_reverse_paired.fq.gz | grep -E "^[ATCGN]+$" | awk '{print length($0)}' | sort | uniq -c | sort -n > 23_4A_reverse_pairedREADS.txt
```
I plan on using R and ggplot2 to make the graphs. Need to scp the read files over to the dir on the computer:

scp ikes@login.talapas.uoregon.edu:/projects/bgmp/ikes/bioinfo/Bi623/QAA/trimcontrol/*READ.txt .

gzcat trimmed_reads.fq.gz | grep -E "^[ATCGN]+$" | awk '{print length($0)}' | sort | uniq -c | sort -n > <insert filename here?>

#script to make dictionary with the file output, Key = length
#graph that


I'll try to write code in Python to do the same sort of thing for the trimmed files
I think these are the files I want to count:
/projects/bgmp/ikes/bioinfo/Bi623/QAA/trimcontrol/1_2A_forward_paired.fq.gz 
reads 7966343
bases 796748295
/projects/bgmp/ikes/bioinfo/Bi623/QAA/trimcontrol/1_2A_reverse_paired.fq.gz
and
/projects/bgmp/ikes/bioinfo/Bi623/QAA/trimcontrol/23_4A_forward_paired.fq.gz
/projects/bgmp/ikes/bioinfo/Bi623/QAA/trimcontrol/23_4A_reverse_paired.fq.gz


???????????? that's not the kind of counting I need. Need help getting Python code counting reads


Date: Mon Sep 11, 2023
##Part 3
8. Install software (record details in lab notebook!!!). In your QAA environment, use conda to install:
star (Spliced Transcripts Alignment to a Reference)
numpy
matplotlib
htseq 

Done on Talapas in my bgmp_qaa environment.

9. Find publicly available mouse genome fasta files (Ensemble release 110) and 
generate an alignment database from them.

When I navigate to:
https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/

I see what appears to be three complete sets of the mouse genome.
Ex:
Mus_musculus.GRCm39.dna.chromosome.1.fa.gz
Mus_musculus.GRCm39.dna_rm.chromosome.1.fa.gz
Mus_musculus.GRCm39.dna_sm.chromosome.1.fa.gz

here are the differences:
dna_sm - Repeats soft-masked (converts repeat nucleotides to lowercase)
dna_rm - Repeats masked (converts repeats to to N's)
dna - No masking

also, there is 
.toplevel - Includes haplotype information (not sure how aligners deal with this)
and
.primary_assembly - Single reference base per position <----------this is the most useful
 
I used 
Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
```
$ wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
```
PS8 told us to get this version:
"Go to the Ensembl website. Navigate to find the zebrafish reference genome by 
chromosome (FASTA) and gene set (GTF) respectively:"

the mouse gene set was at:
https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/

Mus_musculus.GRCm39.110.gtf.gz

I downloaded it to the QAA dir on talapas
```
$ wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
```
My instructions are to use a splice aware aligner. STAR and HISAT2 are splice aware.
"Splice aware" means:
RNA-seq reads are derived from mature mRNA, so there are typically no introns in the sequence. 
But aligners use a reference genome to aid in the process, so a read typically spans two or more 
exons while the reference would have one exon followed by an intron. So the reference genome would 
find a matching sequence in only one of the exons, while the rest of the read would not match 
the intron in the reference, so the read can't be properly aligned. A splice-aware aligner would 
know not to try to align RNA-seq reads to introns, and would somehow identify possible downstream 
exons and try to align to those instead, ignoring introns altogether.

(https://www.biostars.org/p/175454/#:~:text=A%20splice%2Daware%20aligner%20would,those%20instead%2C%20ignoring%20introns%20altogether.)
  
to make STAR database next. Use script adapted from PS8
$ sbatch starreference.sh

Make sure to run sbatch job for alignment:
$ sbatch staralign.sh 


9. Using your script from PS8 in Bi621, report the number of mapped and unmapped reads 
from each of your 2 sam files. Make sure that your script is looking at the bitwise flag 
to determine if reads are primary or secondary mapping (update/fix your script if necessary).

I ran samparser version 1.1 on ikes_1_2A_Aligned.out.sam
15627427 reads were mapped
305259 reads were unmapped
That's 98% mapped

I ran samparser version 1.1 on ikes_24_4A_Aligned.out.sam
79472970 reads were mapped
4640182 reads were unmapped
that's 94.2% mapped

Run htseq to count how many reads map to each feature or gene:

?<1_2A_htseqout_Aligned.tsv>
?<24_4A_htseqout_Aligned.tsv>

Stranded = yes
htseq-count --stranded=yes /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/ikes_1_2A_Aligned.out.sam /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/Mus_musculus.GRCm39.110.gtf > 1_2A_htseqoutStr_Aligned.tsv

htseq-count --stranded=yes /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/ikes_24_4A_Aligned.out.sam /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/Mus_musculus.GRCm39.110.gtf > 24_4A_htseqoutStr_Aligned.tsv

Stranded = reverse
htseq-count --stranded=reverse /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/ikes_1_2A_Aligned.out.sam /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/Mus_musculus.GRCm39.110.gtf > 1_2A_htseqoutUnstr_Aligned.tsv

htseq-count --stranded=reverse /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/ikes_24_4A_Aligned.out.sam /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/Mus_musculus.GRCm39.110.gtf > 24_4A_htseqoutUnstr_Aligned.tsv


Pull down files to HD:

scp ikes@login.talapas.uoregon.edu:/projects/bgmp/ikes/bioinfo/Bi623/QAA/*htseqout*.tsv .

#############################################################

## Stats homework #1 Mon Sep 11, 2023

Note on reminder about how to get stat summary:
call the df (e.g. compounds) and then reference the column of data I want summarized
(e.g. Retention.index)

```
compounds %>% 
  summarise(
    mean = mean(Retention.index),
    sd = sd(Retention.index)
  )
```

I really am having a hard time with this assignment. My loops don't want to iterate. Like, ever.
I have so many syntax issues but it's hard to get an experienced eye to look at my 
code and teach me where I am wrong.

###############################################################################


## DESeq2 Assignment

Date: Tue Sep 12

Get files from 
/projects/bgmp/shared/Bi623/Assignment_DESeq2/
files found there are:
Gacu_gut_counts.tsv  
Gacu_gut_metadata.tsv  
Gasterosteus_aculeatus.BROADS1.85.gtf

pull down to /Bi623 using iTerm2:

scp ikes@login.talapas.uoregon.edu:/projects/bgmp/shared/Bi623/Assignment_DESeq2/*.* .

The counts data was generated using htseq-count (just like you generated in your 
QAA assignment) and is in the file Gacu_gut_counts.tsv. The metadata containing 
information about each sample (i.e. column in the counts file) can be found in 
Gacu_gut_metadata.tsv.

I Installed DESeq2 into R:
```
>BiocManager::install("DESeq2")

```
and I answered "a" to update all files when prompted

Read each of these into R dataframes named countdata and coldata respectively. Each column 
in the metadata table should be a factor (i.e. category). Hint: What happened to the column names in countdata?

Vocab learned:
When the expected amount of variance is approximately the same across different mean values, 
the data is said to be homoskedastic

#3: Note from code: 
rlog() may take a few minutes with 30 or more samples,
vst() is a much faster transformation

Copied from workflow: 
https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#starting-from-count-matrices

As a solution, DESeq2 offers two transformations for count data that stabilize the variance across the mean: the variance stabilizing transformation (VST) for negative binomial data with a dispersion-mean trend (Anders and Huber 2010), implemented in the vst function, and the regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014).

Which transformation to choose? The VST is much faster to compute and is less sensitive to high count outliers than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples (an order of magnitude difference). We therefore recommend the VST for medium-to-large datasets (n > 30).

