# Workflow for processing frog 2bRAD parentage data in collaboration with Brandon Guell and Warkentin Lab
# Author: Hannah E. Aichelman, edited by Isabel Novick


#------------------------------GET DATA FROM TUFTS
# Go into directory that you want to move files to 

[inovick@scc1 tucf_data]$ pwd
/projectnb/mullenl/novick/global_pop_gen/tucf_data



# Run the wget command to pull files from Tufts:

wget -r -nH --cut-dirs=1 -nc ftp://sean.mullen:jidXxGK@130.64.74.72//231108-0165_Sean_Mullen-7878-rename/
wget -r -nH --cut-dirs=1 -nc ftp://sean.mullen:jidXxGK@130.64.74.72//231108-0164_Sean_Mullen-7878-rename/


#------------------------------CHECK THAT DATA DOWNLOADED CORRECTLY

# Run md5sum to check that all of the files transferred properly, separately within each of the original, zipped folders for 64 and 65.

drwxr-sr-x 5 inovick mullenl 8.0K Jan  8 16:00 231108-0164_Sean_Mullen-7878-rename
drwxr-sr-x 5 inovick mullenl 4.0K Jan  8 14:35 231108-0165_Sean_Mullen-7878-rename

md5sum *.fastq* > 64_rename_md5sum.md5 #within 231108-0164_Sean_Mullen-7878-rename
md5sum *.fastq* > 65_rename_md5sum.md5 #within 231108-0165_Sean_Mullen-7878-rename

# Copy and paste into excel file with one column from the md5sum.txt file and the other from the file we just made. In middle, put =[cell 1]=[cell 2] and it will say true or false if they match or don't match

# Compare the outputs to the md5sum.txt file from the sequencing facility - they match!


#------------------------------UNZIP FILES

# Copy and unzip files into new directory

[inovick@scc1 rad_workflow]$ pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow

# Make new directories for unzipped data to go into

mkdir 64_unzipped
mkdir 65_unzipped

# Commands for copying and unzipping data found in unzip_data.qsub

cp ../tucf_data/231108-0164_Sean_Mullen-7878-rename/*.fastq.gz /projectnb/mullenl/novick/global_pop_gen/rad_workflow/64_unzipped
cp ../tucf_data/231108-0165_Sean_Mullen-7878-rename/*.fastq.gz /projectnb/mullenl/novick/global_pop_gen/rad_workflow/65_unzipped

gunzip /projectnb/mullenl/novick/global_pop_gen/rad_workflow/64_unzipped *.fastq.gz
gunzip /projectnb/mullenl/novick/global_pop_gen/rad_workflow/65_unzipped *.fastq.gz


#------------------------------DEMULTIPLEXING, TRIMMING AND MAPPING 

#Download scripts from https://github.com/z0on/2bRAD_denovo

[inovick@scc1 global_pop_gen]$ pwd
/projectnb/mullenl/novick/global_pop_gen

[inovick@scc1 global_pop_gen]$ git clone https://github.com/z0on/2bRAD_denovo.git


# Make everything executable in the 2bRAD_denovo folder:

chmod +x *.pl
chmod +x *.py
chmod +x *.R

# Install modules - we will need bowtie2, samtools, and picard. 

module load perl
module load bowtie2
module load samtools
module load picard


# Using the reference genome from Jasmine

cd rad_workflow/
pwd 
/projectnb/mullenl/novick/global_pop_gen/rad_workflow
cp /projectnb/mullenl/alqassar/wcm_annotation/wcm_pseudochromosomes_renamed.fasta .

ls -lh

drwxr-sr-x 2 inovick mullenl 4.0K Jan  9 11:48 64_unzipped
drwxr-sr-x 2 inovick mullenl 4.0K Jan  9 12:01 65_unzipped
-rw-r--r-- 1 inovick mullenl 2.4K Jan  9 12:02 copy_and_unzip.qlog
-rw-r--r-- 1 inovick mullenl 1.1K Jan  8 16:04 unzip_data.qsub
-rw-r--r-- 1 inovick mullenl 236M Jan 11 12:05 wcm_pseudochromosomes_renamed.fasta

# this assigns the value of the file to this variable so I don't have to write it out every time

export GENOME_FASTA=wcm_pseudochromosomes_renamed.fasta

# indexing genome for bowtie2 mapper
bowtie2-build $GENOME_FASTA $GENOME_FASTA

# It seems like it's not finding any A's for some reason? Will this be a problem?

samtools faidx $GENOME_FASTA

#==================
# Step 1: Splitting by in-read barcode, deduplicating and quality-filtering the reads

# Making a new directory called trimming in rad_workflow, and copying all fastq files from 64_unzipped and 65_unzipped into that one directory

[inovick@scc1 rad_workflow]$ mkdir trimming
[inovick@scc1 rad_workflow]$ cd 64_unzipped/
[inovick@scc1 64_unzipped]$ cp *fastq ../trimming
[inovick@scc1 65_unzipped]$ cp *fastq ../trimming

[inovick@scc1 trimming]$ pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/trimming

# creating a file of commands to run (assuming reads are in fastq files, one file per sample.)
 
../../2bRAD_denovo/2bRAD_trim_launch_dedup.pl fastq > trims

# Modify this trims file to be able to submit it as a job on the SCC, designate where the perl script is, and add '&'s to the end of the line
# It looks like this:

cat trims

#Copy and paste readout into a new qsub file, naming it split_by_barcode.qsub

touch split_by_barcode.qsub

#Copy the following into new qsub file:

#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N split_by_barcodes # job name, anything you want
#$ -o split_by_barcodes.qlog #Here is the logfile
#$ -l h_rt=24:00:00 #maximum run time
#$ -M inovick@bu.edu #your email
#$ -m be

module load perl

../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L2_4_S12_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L3_6_S6_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L4_2_S10_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L1_5_S5_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L1_1_S1_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L3_2_S2_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L2_6_S14_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L3_1_S1_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L1_2_S2_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L4_3_S11_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L3_4_S4_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L2_5_S13_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=65_Undetermined_S0_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L2_1_S9_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L2_8_S16_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L2_3_S11_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L4_1_S9_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L4_4_S12_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L1_3_S3_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L1_8_S8_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=Undetermined_S0_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L3_3_S3_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L3_8_S8_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L3_7_S7_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L2_2_S10_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L3_5_S5_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L1_6_S6_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=65_Undetermined_S0_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L1_7_S7_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L2_7_S15_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=IAN_2brad_L1_4_S4_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=Undetermined_S0_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &

wait


# do we have expected number of *.tr0 files created?
ls -l *.tr0 | wc -l

#No, this was probably the wrong script. Going to now try 

../../2bRAD_denovo/2bRAD_trim_launch_dedup_old.pl fastq > trims

cat trims

#Copy and paste into old qsub, replacing previous trims

../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L2_4_S12_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L3_6_S6_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L4_2_S10_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L1_5_S5_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L1_1_S1_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L3_2_S2_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L2_6_S14_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L3_1_S1_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L1_2_S2_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L4_3_S11_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L3_4_S4_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L2_5_S13_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=65_Undetermined_S0_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L2_1_S9_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L2_8_S16_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L2_3_S11_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L4_1_S9_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L4_4_S12_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L1_3_S3_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L1_8_S8_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=Undetermined_S0_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L3_3_S3_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L3_8_S8_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L3_7_S7_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L2_2_S10_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L3_5_S5_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L1_6_S6_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=65_Undetermined_S0_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L1_7_S7_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L2_7_S15_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=IAN_2brad_L1_4_S4_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../../2bRAD_denovo/trim2bRAD_2barcodes_dedup_old.pl input=Undetermined_S0_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &

# This still didn't work, no tr0 files are being created, qlog file is empty, wallclock runtime is 0
#Ahhh it's because you have to put "wait" at the end and I was ommitting that, lets see if it works better now
#Just ran the original, not old, script again with "wait" at the end
#Ran for longer when I put wait at the end, but still didnt create any tr0 files and in the qlog, everything came up 0s?

###Now when using the "old" script it worked! 

# do we have expected number of *.tr0 files created?
ls -l *.tr0 | wc -l

373

# I originally started with 327 samples so there are more tr0 files than samples
# Subtracting those that are from undetermined, I have 325 tr0 files

# for reference-based analysis: trimming poor quality bases off ends:

>trimse
for file in *.tr0; do
echo "cutadapt -q 15,15 -m 25 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse_2;
done

# add job header to the top of the file and qsub
# looks like this:

cat trimse_2

# copy and paste output into new qsub


#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N trimse # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M inovick@bu.edu #your email
#$ -m be

module load cutadapt

cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim IAN_2brad_L1_1_S1_L001_R1_001_ACCA.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_ACCA.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim IAN_2brad_L1_1_S1_L001_R1_001_AGAC.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_AGAC.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim IAN_2brad_L1_1_S1_L001_R1_001_AGTG.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_AGTG.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim IAN_2brad_L1_1_S1_L001_R1_001_CATC.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_CATC.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim IAN_2brad_L1_1_S1_L001_R1_001_CTAC.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_CTAC.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim IAN_2brad_L1_1_S1_L001_R1_001_GACT.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_GACT.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim IAN_2brad_L1_1_S1_L001_R1_001_GCTT.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_GCTT.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim IAN_2brad_L1_1_S1_L001_R1_001_GTGA.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_GTGA.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim IAN_2brad_L1_1_S1_L001_R1_001_TCAC.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_TCAC.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim IAN_2brad_L1_1_S1_L001_R1_001_TCAG.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_TCAG.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_TGGT.trim IAN_2brad_L1_1_S1_L001_R1_001_TGGT.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_TGGT.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_1_S1_L001_R1_001_TGTC.trim IAN_2brad_L1_1_S1_L001_R1_001_TGTC.tr0 > IAN_2brad_L1_1_S1_L001_R1_001_TGTC.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_2_S2_L001_R1_001_ACCA.trim IAN_2brad_L1_2_S2_L001_R1_001_ACCA.tr0 > IAN_2brad_L1_2_S2_L001_R1_001_ACCA.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_2_S2_L001_R1_001_AGAC.trim IAN_2brad_L1_2_S2_L001_R1_001_AGAC.tr0 > IAN_2brad_L1_2_S2_L001_R1_001_AGAC.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_2_S2_L001_R1_001_AGTG.trim IAN_2brad_L1_2_S2_L001_R1_001_AGTG.tr0 > IAN_2brad_L1_2_S2_L001_R1_001_AGTG.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_2_S2_L001_R1_001_CATC.trim IAN_2brad_L1_2_S2_L001_R1_001_CATC.tr0 > IAN_2brad_L1_2_S2_L001_R1_001_CATC.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_2_S2_L001_R1_001_CTAC.trim IAN_2brad_L1_2_S2_L001_R1_001_CTAC.tr0 > IAN_2brad_L1_2_S2_L001_R1_001_CTAC.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_2_S2_L001_R1_001_GACT.trim IAN_2brad_L1_2_S2_L001_R1_001_GACT.tr0 > IAN_2brad_L1_2_S2_L001_R1_001_GACT.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_2_S2_L001_R1_001_GCTT.trim IAN_2brad_L1_2_S2_L001_R1_001_GCTT.tr0 > IAN_2brad_L1_2_S2_L001_R1_001_GCTT.tr0_trimlog.txt
cutadapt --format fastq -q 15,15 -m 25 -o IAN_2brad_L1_2_S2_L001_R1_001_GTGA.trim IAN_2brad_L1_2_S2_L001_R1_001_GTGA.tr0 > IAN_2brad_L1_2_S2_L001_R1_001_GTGA.tr0_trimlog.txt

# do we have expected number of *.trim files created?
ls -l *.trim | wc -l

# It made 0, it only made trimlog.txt. So I repeated this step

for file in *.tr0; do
echo "cutadapt -q 15,15 -m 25 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse_2;
done

# but removed "format"
# now added header to trimse_2, running trimse_2 hoping that it makes .trim files, YOU MUST RUN IN DIRECTORY WITH .TR0 FILES (/trimming)

/projectnb/mullenl/novick/global_pop_gen/rad_workflow/trimming

# do we have expected number of *.trim files created?
ls -l *.trim | wc -l

325

#yes

# for reference-based: 

pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/trimming


# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
# made qsub for this called bowtie2.qsub in /trimming

touch bowtie2.qsub

# nano into the file, add header and the following line:

export GENOME_FASTA=wcm_pseudochromosomes_renamed.fasta
/projectnb/mullenl/novick/global_pop_gen/2bRAD_denovo/2bRAD_bowtie2_launch.pl '\.trim$' $GENOME_FASTA > maps

# Once it ran, check out the output

head -20 maps

#Make sure the genome is bowtie readable using the following command, make into a qsub file with header:

touch make_genome_readable.qsub

#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N genome_readable # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M inovick@bu.edu #your email
#$ -m be

module load bowtie2
module load samtools

export GENOME_FASTA=wcm_pseudochromosomes_renamed.fasta

bowtie2-build $GENOME_FASTA $GENOME_FASTA

samtools faidx $GENOME_FASTA


# same as you did above for trimse, nano into maps once you created it and add your header and your 'module load line'.
# add header and then submit maps job

#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N maps # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M inovick@bu.edu #your email
#$ -m be

module load bowtie2

bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x wcm_pseudochromosomes_renamed.fasta -U IAN_2brad_L3_6_S6_L001_R1_001_GACT.trim -S IAN_2brad_L3_6_S6_L001_R1_001_GACT.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x wcm_pseudochromosomes_renamed.fasta -U IAN_2brad_L1_6_S6_L001_R1_001_ACCA.trim -S IAN_2brad_L1_6_S6_L001_R1_001_ACCA.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x wcm_pseudochromosomes_renamed.fasta -U IAN_2brad_L2_3_S11_L002_R1_001_GTGA.trim -S IAN_2brad_L2_3_S11_L002_R1_001_GTGA.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x wcm_pseudochromosomes_renamed.fasta -U IAN_2brad_L2_3_S11_L002_R1_001_TCAG.trim -S IAN_2brad_L2_3_S11_L002_R1_001_TCAG.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x wcm_pseudochromosomes_renamed.fasta -U IAN_2brad_L2_7_S15_L002_R1_001_CTAC.trim -S IAN_2brad_L2_7_S15_L002_R1_001_CTAC.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x wcm_pseudochromosomes_renamed.fasta -U IAN_2brad_L2_4_S12_L002_R1_001_TGGT.trim -S IAN_2brad_L2_4_S12_L002_R1_001_TGGT.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x wcm_pseudochromosomes_renamed.fasta -U IAN_2brad_L2_2_S10_L002_R1_001_CATC.trim -S IAN_2brad_L2_2_S10_L002_R1_001_CATC.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x wcm_pseudochromosomes_renamed.fasta -U IAN_2brad_L1_3_S3_L001_R1_001_TGTC.trim -S IAN_2brad_L1_3_S3_L001_R1_001_TGTC.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x wcm_pseudochromosomes_renamed.fasta -U IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim -S IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x wcm_pseudochromosomes_renamed.fasta -U IAN_2brad_L2_4_S12_L002_R1_001_AGAC.trim -S IAN_2brad_L2_4_S12_L002_R1_001_AGAC.trim.bt2.sam


# submit job
qsub maps

ls *trim.bt2.sam > sams
cat sams | wc -l  # number should match number of trim files (325)
# 325 - matches

# check alignment rates:
>alignmentRates
for F in `ls *trim`; do 
M=`grep -E '^[ATGCN]+$' $F | wc -l | grep -f - maps.e* -A 4 | tail -1 | perl -pe 's/maps\.e\d+-|% overall alignment rate//g'` ;
echo "$F.sam $M">>alignmentRates;
done

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:

#Hannah's:
module load samtools
cat sams | perl -pe 's/(\S+)\.sam/samtools view -bS $1\.sam >$1\.unsorted\.bam && samtools sort $1\.unsorted\.bam -o $1\.bam && samtools index $1\.bam/' >s2b


# James and Misha's:
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done



head -20 s2b

# add header
#The following example was made using Hannah's s2b script


#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N s2b # job name, anything you want
#$ -o s2b.qlog #Here is the logfile
#$ -l h_rt=24:00:00 #maximum run time
#$ -M inovick@bu.edu #your email
#$ -m be

module load samtools

samtools view -bS IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim.bt2.sam >IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim.bt2.unsorted.bam && samtools sort IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim.bt2.unsorted.bam -o IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim.bt2.bam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim.bt2.bam
samtools view -bS IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim.bt2.sam >IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim.bt2.unsorted.bam && samtools sort IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim.bt2.unsorted.bam -o IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim.bt2.bam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim.bt2.bam
samtools view -bS IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim.bt2.sam >IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim.bt2.unsorted.bam && samtools sort IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim.bt2.unsorted.bam -o IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim.bt2.bam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim.bt2.bam
samtools view -bS IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim.bt2.sam >IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim.bt2.unsorted.bam && samtools sort IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim.bt2.unsorted.bam -o IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim.bt2.bam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim.bt2.bam
samtools view -bS IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.sam >IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.unsorted.bam && samtools sort IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.unsorted.bam -o IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam
samtools view -bS IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim.bt2.sam >IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim.bt2.unsorted.bam && samtools sort IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim.bt2.unsorted.bam -o IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim.bt2.bam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim.bt2.bam
samtools view -bS IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim.bt2.sam >IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim.bt2.unsorted.bam && samtools sort IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim.bt2.unsorted.bam -o IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim.bt2.bam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim.bt2.bam
samtools view -bS IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim.bt2.sam >IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim.bt2.unsorted.bam && samtools sort IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim.bt2.unsorted.bam -o IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim.bt2.bam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim.bt2.bam

# and then submit s2b job
qsub s2b

# remove everything unsorted
rm -f *unsorted*



ls *bam | wc -l  # should be the same number as number of .fastq files
# I did not have enough, I got 79. Going to delete all bams and recreate s2b using Misha's script this time

for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done


#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N s2b # job name, anything you want
#$ -o s2b.qlog #Here is the logfile
#$ -l h_rt=24:00:00 #maximum run time
#$ -M inovick@bu.edu #your email
#$ -m be

module load samtools

samtools sort -O bam -o IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim.bt2.bam IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim.bt2.sam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim.bt2.bam
samtools sort -O bam -o IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim.bt2.bam IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim.bt2.sam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim.bt2.bam
samtools sort -O bam -o IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim.bt2.bam IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim.bt2.sam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim.bt2.bam
samtools sort -O bam -o IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim.bt2.bam IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim.bt2.sam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim.bt2.bam
samtools sort -O bam -o IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.sam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam
samtools sort -O bam -o IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim.bt2.bam IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim.bt2.sam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim.bt2.bam
samtools sort -O bam -o IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim.bt2.bam IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim.bt2.sam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim.bt2.bam
samtools sort -O bam -o IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim.bt2.bam IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim.bt2.sam && samtools index IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim.bt2.bam

# This worked. 325 bams

# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis. Archive them.

# To archive the bams:

[inovick@scc1 trimming]$ pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/trimming

cd ../
mkdir bams
cd trimming
cp *.bam ../bams

cd ../

[inovick@scc1 rad_workflow]$ pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow

tar -zcvf bams.tar.gz bams

# Now delete original nonzipped bams folder
rm -rf bams


#------------------------------ANGSD (fuzzy genotyping)
# this first run-through of ANGSD is just to take a look at base qualities and coverage depth, 
# we will run angsd again with filters informed by this first run-through

# Make new ANGSD directory
[inovick@scc1 rad_workflow]$ pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow

[inovick@scc1 rad_workflow]$ mkdir ANGSD
[inovick@scc1 rad_workflow]$ cd trimming/
[inovick@scc1 trimming]$ cp *.bam ../ANGSD/

# "FUZZY genotyping" with ANGSD - without calling actual genotypes but working with genotype likelihoods at each SNP. Optimal for low-coverage data (<10x).

# listing all bam filenames 
[inovick@scc1 ANGSD]$ pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/ANGSD


ls *bam >bams

module load angsd

#----------- assessing base qualities and coverage depth

# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping =< 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd : the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 3250 -minInd 163"

# T O   D O : 
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

# in the following line, -r argument is one chromosome or contig to work with (no need to do this for whole genome as long as the chosen chromosome or contig is long enough, ~1 Mb. When mapping to a real genome, consider chr1:1-1000000 )
# (look up lengths of your contigs in the header of *.sam files)
# my pseudochromosomes are very long so no need to use -r argument
# i'm not using here and just running on entire reference

angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out dd

# Submit this job as a qsub

touch angsd1.qsub
nano angsd1.qsub

#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ANGSD1 # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M inovick@bu.edu
#$ -m be

module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 3250 -minInd 163"
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out dd

# output file says:
-> Total number of sites analyzed: 3019468
-> Number of sites retained after filtering: 1142751

# summarizing results (using modified script by Matteo Fumagalli)
module load R
Rscript ../../2bRAD_denovo/plotQC.R prefix=dd

# proportion of sites covered at >5x:
cat quality.txt

IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam 0.0110273010437435
IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam 0.0193986356835845
IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam 0.0416131656037242
IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam 0.0508641375683166
IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam 0.0993919746965222
IAN_2brad_L1_7_S7_L001_R1_001_GACT.trim.bt2.bam 0.11513661706499
IAN_2brad_L1_4_S4_L001_R1_001_AGAC.trim.bt2.bam 0.124160484766929
IAN_2brad_L1_7_S7_L001_R1_001_TCAC.trim.bt2.bam 0.276104952113743
IAN_2brad_L1_5_S5_L001_R1_001_TGGT.trim.bt2.bam 0.441787876477245
IAN_2brad_L2_6_S14_L002_R1_001_CTAC.trim.bt2.bam 0.504444609153927
IAN_2brad_L3_8_S8_L001_R1_001_GTGA.trim.bt2.bam 0.532844586778205
IAN_2brad_L1_6_S6_L001_R1_001_ACCA.trim.bt2.bam 0.552969362312665
IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim.bt2.bam 0.578282219992728
IAN_2brad_L1_3_S3_L001_R1_001_GCTT.trim.bt2.bam 0.581139473163754
IAN_2brad_L3_4_S4_L001_R1_001_TGGT.trim.bt2.bam 0.584522207539318
IAN_2brad_L3_6_S6_L001_R1_001_ACCA.trim.bt2.bam 0.602030252968755
IAN_2brad_L1_4_S4_L001_R1_001_GACT.trim.bt2.bam 0.639252637863226
IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim.bt2.bam 0.639323960324049
IAN_2brad_L4_1_S9_L002_R1_001_TGTC.trim.bt2.bam 0.644452545606899
IAN_2brad_L2_2_S10_L002_R1_001_GACT.trim.bt2.bam 0.646494563700585
IAN_2brad_L1_8_S8_L001_R1_001_ACCA.trim.bt2.bam 0.647020584775678
IAN_2brad_L3_6_S6_L001_R1_001_CATC.trim.bt2.bam 0.654766781818008
IAN_2brad_L4_2_S10_L002_R1_001_AGAC.trim.bt2.bam 0.657654975400666
IAN_2brad_L1_2_S2_L001_R1_001_GCTT.trim.bt2.bam 0.657695908614573
IAN_2brad_L1_2_S2_L001_R1_001_GACT.trim.bt2.bam 0.660735063435546
IAN_2brad_L1_3_S3_L001_R1_001_TCAC.trim.bt2.bam 0.661056867371821
IAN_2brad_L1_8_S8_L001_R1_001_GACT.trim.bt2.bam 0.667101846593469
IAN_2brad_L3_3_S3_L001_R1_001_CATC.trim.bt2.bam 0.674444859792815
IAN_2brad_L2_2_S10_L002_R1_001_TCAG.trim.bt2.bam 0.675896416483859
IAN_2brad_L2_6_S14_L002_R1_001_AGAC.trim.bt2.bam 0.68427294071103
IAN_2brad_L2_5_S13_L002_R1_001_ACCA.trim.bt2.bam 0.687815370801389
IAN_2brad_L1_7_S7_L001_R1_001_TGTC.trim.bt2.bam 0.68972398334533
IAN_2brad_L1_3_S3_L001_R1_001_GTGA.trim.bt2.bam 0.699020910889597
IAN_2brad_L3_8_S8_L001_R1_001_AGTG.trim.bt2.bam 0.707982439887556
IAN_2brad_L2_4_S12_L002_R1_001_ACCA.trim.bt2.bam 0.711162189378311
IAN_2brad_L1_5_S5_L001_R1_001_TCAC.trim.bt2.bam 0.713490026236135
IAN_2brad_L2_2_S10_L002_R1_001_CATC.trim.bt2.bam 0.716000290088041
IAN_2brad_L2_2_S10_L002_R1_001_GTGA.trim.bt2.bam 0.723647043581287
IAN_2brad_L3_6_S6_L001_R1_001_TCAG.trim.bt2.bam 0.726568095557084
IAN_2brad_L3_1_S1_L001_R1_001_AGAC.trim.bt2.bam 0.727340376446712
IAN_2brad_L1_8_S8_L001_R1_001_CTAC.trim.bt2.bam 0.728438522324538
IAN_2brad_L1_8_S8_L001_R1_001_AGAC.trim.bt2.bam 0.733690878526896
IAN_2brad_L2_2_S10_L002_R1_001_GCTT.trim.bt2.bam 0.735576780528577
IAN_2brad_L3_8_S8_L001_R1_001_AGAC.trim.bt2.bam 0.737423772798394
IAN_2brad_L1_8_S8_L001_R1_001_TGTC.trim.bt2.bam 0.738474778350037
IAN_2brad_L3_2_S2_L001_R1_001_TCAG.trim.bt2.bam 0.744746084124699
IAN_2brad_L2_4_S12_L002_R1_001_TGGT.trim.bt2.bam 0.745909124343603
IAN_2brad_L1_3_S3_L001_R1_001_TCAG.trim.bt2.bam 0.747517497620083
IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim.bt2.bam 0.747560731003873
IAN_2brad_L1_6_S6_L001_R1_001_CTAC.trim.bt2.bam 0.748134730878088
IAN_2brad_L3_3_S3_L001_R1_001_ACCA.trim.bt2.bam 0.748630240328278
IAN_2brad_L1_7_S7_L001_R1_001_AGTG.trim.bt2.bam 0.755176113813629
IAN_2brad_L1_5_S5_L001_R1_001_GACT.trim.bt2.bam 0.756432232805663
IAN_2brad_L3_8_S8_L001_R1_001_TCAC.trim.bt2.bam 0.756452438765934
IAN_2brad_L1_1_S1_L001_R1_001_TGGT.trim.bt2.bam 0.757360828451143
IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim.bt2.bam 0.759153727991264
IAN_2brad_L3_1_S1_L001_R1_001_GTGA.trim.bt2.bam 0.759704792725249
IAN_2brad_L2_6_S14_L002_R1_001_GCTT.trim.bt2.bam 0.759852542625128
IAN_2brad_L2_2_S10_L002_R1_001_CTAC.trim.bt2.bam 0.760924187227391
IAN_2brad_L2_2_S10_L002_R1_001_ACCA.trim.bt2.bam 0.764790814329325
IAN_2brad_L1_6_S6_L001_R1_001_AGTG.trim.bt2.bam 0.765089911752552
IAN_2brad_L2_1_S9_L002_R1_001_AGTG.trim.bt2.bam 0.766981411821085
IAN_2brad_L3_1_S1_L001_R1_001_GCTT.trim.bt2.bam 0.76711633029978
IAN_2brad_L3_2_S2_L001_R1_001_GCTT.trim.bt2.bam 0.771416382498508
IAN_2brad_L3_3_S3_L001_R1_001_TCAG.trim.bt2.bam 0.774792472170357
IAN_2brad_L4_4_S12_L002_R1_001_TCAG.trim.bt2.bam 0.775374018869478
IAN_2brad_L2_4_S12_L002_R1_001_AGTG.trim.bt2.bam 0.776167220307916
IAN_2brad_L1_4_S4_L001_R1_001_GCTT.trim.bt2.bam 0.776212107202393
IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim.bt2.bam 0.777155463185937
IAN_2brad_L1_7_S7_L001_R1_001_GTGA.trim.bt2.bam 0.777601418830858
IAN_2brad_L2_2_S10_L002_R1_001_AGAC.trim.bt2.bam 0.777839069947774
IAN_2brad_L1_5_S5_L001_R1_001_GCTT.trim.bt2.bam 0.778331428002199
IAN_2brad_L1_3_S3_L001_R1_001_AGTG.trim.bt2.bam 0.7790208782967
IAN_2brad_L4_2_S10_L002_R1_001_AGTG.trim.bt2.bam 0.779380506998269
IAN_2brad_L2_2_S10_L002_R1_001_AGTG.trim.bt2.bam 0.780787193463993
IAN_2brad_L3_3_S3_L001_R1_001_AGTG.trim.bt2.bam 0.781444006930413
IAN_2brad_L3_2_S2_L001_R1_001_AGTG.trim.bt2.bam 0.782048077635692
IAN_2brad_L2_6_S14_L002_R1_001_ACCA.trim.bt2.bam 0.784322198742386
IAN_2brad_L2_7_S15_L002_R1_001_ACCA.trim.bt2.bam 0.785091660644642
IAN_2brad_L3_2_S2_L001_R1_001_CTAC.trim.bt2.bam 0.785224123746442
IAN_2brad_L4_1_S9_L002_R1_001_CTAC.trim.bt2.bam 0.785534562043991
IAN_2brad_L2_4_S12_L002_R1_001_CTAC.trim.bt2.bam 0.786540271323507
IAN_2brad_L2_2_S10_L002_R1_001_TGGT.trim.bt2.bam 0.787050542190371
IAN_2brad_L1_8_S8_L001_R1_001_GTGA.trim.bt2.bam 0.787172982998558
IAN_2brad_L3_5_S5_L001_R1_001_AGTG.trim.bt2.bam 0.787244638541563
IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim.bt2.bam 0.787464561137643
IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim.bt2.bam 0.789268628533152
IAN_2brad_L3_8_S8_L001_R1_001_TGGT.trim.bt2.bam 0.789341114970005
IAN_2brad_L3_4_S4_L001_R1_001_GTGA.trim.bt2.bam 0.789760183905279
IAN_2brad_L3_5_S5_L001_R1_001_GCTT.trim.bt2.bam 0.790095306296866
IAN_2brad_L3_8_S8_L001_R1_001_ACCA.trim.bt2.bam 0.790964517475671
IAN_2brad_L3_4_S4_L001_R1_001_GCTT.trim.bt2.bam 0.792297061714829
IAN_2brad_L3_5_S5_L001_R1_001_TGTC.trim.bt2.bam 0.793518604085517
IAN_2brad_L2_8_S16_L002_R1_001_AGTG.trim.bt2.bam 0.794234230857369
IAN_2brad_L2_8_S16_L002_R1_001_TGGT.trim.bt2.bam 0.794267442578908
IAN_2brad_L1_3_S3_L001_R1_001_TGTC.trim.bt2.bam 0.794767906930782
IAN_2brad_L3_8_S8_L001_R1_001_TCAG.trim.bt2.bam 0.794940532618209
IAN_2brad_L1_3_S3_L001_R1_001_CTAC.trim.bt2.bam 0.795466849558716
IAN_2brad_L1_7_S7_L001_R1_001_ACCA.trim.bt2.bam 0.796456690926511
IAN_2brad_L2_6_S14_L002_R1_001_TGTC.trim.bt2.bam 0.796809368233921
IAN_2brad_L3_3_S3_L001_R1_001_TCAC.trim.bt2.bam 0.798105637375904
IAN_2brad_L1_4_S4_L001_R1_001_TGGT.trim.bt2.bam 0.798910989165661
IAN_2brad_L3_7_S7_L001_R1_001_TCAG.trim.bt2.bam 0.799945714252206
IAN_2brad_L1_3_S3_L001_R1_001_GACT.trim.bt2.bam 0.801049087264774
IAN_2brad_L1_3_S3_L001_R1_001_CATC.trim.bt2.bam 0.801401298461091
IAN_2brad_L3_7_S7_L001_R1_001_GCTT.trim.bt2.bam 0.801873455637415
IAN_2brad_L3_1_S1_L001_R1_001_TGTC.trim.bt2.bam 0.802388118119948
IAN_2brad_L2_6_S14_L002_R1_001_AGTG.trim.bt2.bam 0.802741805984143
IAN_2brad_L2_1_S9_L002_R1_001_TGGT.trim.bt2.bam 0.803877249993702
IAN_2brad_L2_5_S13_L002_R1_001_GCTT.trim.bt2.bam 0.804680869438689
IAN_2brad_L3_3_S3_L001_R1_001_CTAC.trim.bt2.bam 0.805019058916063
IAN_2brad_L1_8_S8_L001_R1_001_AGTG.trim.bt2.bam 0.806256142176614
IAN_2brad_L1_6_S6_L001_R1_001_TCAG.trim.bt2.bam 0.807647277875169
IAN_2brad_L4_3_S11_L002_R1_001_GCTT.trim.bt2.bam 0.808024637378249
IAN_2brad_L1_7_S7_L001_R1_001_TGGT.trim.bt2.bam 0.808196496008806
IAN_2brad_L1_2_S2_L001_R1_001_TGTC.trim.bt2.bam 0.808486474966526
IAN_2brad_L4_1_S9_L002_R1_001_TCAG.trim.bt2.bam 0.809592462132682
IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim.bt2.bam 0.809954876802667
IAN_2brad_L1_1_S1_L001_R1_001_TGTC.trim.bt2.bam 0.811703898324009
IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim.bt2.bam 0.812523939912288
IAN_2brad_L1_7_S7_L001_R1_001_GCTT.trim.bt2.bam 0.812692433838053
IAN_2brad_L3_6_S6_L001_R1_001_GCTT.trim.bt2.bam 0.816763342960996
IAN_2brad_L2_3_S11_L002_R1_001_AGTG.trim.bt2.bam 0.817048625315689
IAN_2brad_L3_4_S4_L001_R1_001_TGTC.trim.bt2.bam 0.817632527464471
IAN_2brad_L2_8_S16_L002_R1_001_GCTT.trim.bt2.bam 0.818614636168374
IAN_2brad_L3_7_S7_L001_R1_001_GACT.trim.bt2.bam 0.819005979901528
IAN_2brad_L3_1_S1_L001_R1_001_AGTG.trim.bt2.bam 0.819570480450389
IAN_2brad_L4_3_S11_L002_R1_001_AGTG.trim.bt2.bam 0.819959050908614
IAN_2brad_L1_6_S6_L001_R1_001_TGTC.trim.bt2.bam 0.820555171653506
IAN_2brad_L3_4_S4_L001_R1_001_TCAG.trim.bt2.bam 0.821316506165452
IAN_2brad_L2_1_S9_L002_R1_001_CATC.trim.bt2.bam 0.821406310848283
IAN_2brad_L3_1_S1_L001_R1_001_CATC.trim.bt2.bam 0.821681088059649
IAN_2brad_L3_5_S5_L001_R1_001_ACCA.trim.bt2.bam 0.822017376094755
IAN_2brad_L2_1_S9_L002_R1_001_GCTT.trim.bt2.bam 0.823556348446225
IAN_2brad_L3_4_S4_L001_R1_001_AGTG.trim.bt2.bam 0.824164932710329
IAN_2brad_L1_7_S7_L001_R1_001_TCAG.trim.bt2.bam 0.824455231102641
IAN_2brad_L2_5_S13_L002_R1_001_GACT.trim.bt2.bam 0.824554241856667
IAN_2brad_L3_2_S2_L001_R1_001_GTGA.trim.bt2.bam 0.824827782798071
IAN_2brad_L2_1_S9_L002_R1_001_TGTC.trim.bt2.bam 0.82491149150727
IAN_2brad_L2_1_S9_L002_R1_001_GACT.trim.bt2.bam 0.824924306347169
IAN_2brad_L2_3_S11_L002_R1_001_ACCA.trim.bt2.bam 0.825010595802468
IAN_2brad_L1_4_S4_L001_R1_001_AGTG.trim.bt2.bam 0.825477865890807
IAN_2brad_L2_3_S11_L002_R1_001_TCAG.trim.bt2.bam 0.825731890965465
IAN_2brad_L3_3_S3_L001_R1_001_GCTT.trim.bt2.bam 0.826289179803111
IAN_2brad_L3_2_S2_L001_R1_001_AGAC.trim.bt2.bam 0.826972460434667
IAN_2brad_L2_1_S9_L002_R1_001_AGAC.trim.bt2.bam 0.828130617274276
IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim.bt2.bam 0.829918438907725
IAN_2brad_L2_3_S11_L002_R1_001_GCTT.trim.bt2.bam 0.83001341361634
IAN_2brad_L4_2_S10_L002_R1_001_TGTC.trim.bt2.bam 0.830200529963598
IAN_2brad_L2_1_S9_L002_R1_001_TCAG.trim.bt2.bam 0.830454286702142
IAN_2brad_L1_4_S4_L001_R1_001_TCAG.trim.bt2.bam 0.83121044766326
IAN_2brad_L3_2_S2_L001_R1_001_GACT.trim.bt2.bam 0.831612477958383
IAN_2brad_L2_6_S14_L002_R1_001_GTGA.trim.bt2.bam 0.831620763345803
IAN_2brad_L2_4_S12_L002_R1_001_TCAG.trim.bt2.bam 0.83177523936719
IAN_2brad_L4_1_S9_L002_R1_001_TCAC.trim.bt2.bam 0.831800047817462
IAN_2brad_L2_1_S9_L002_R1_001_TCAC.trim.bt2.bam 0.832249567774342
IAN_2brad_L4_4_S12_L002_R1_001_TGTC.trim.bt2.bam 0.83263159278109
IAN_2brad_L2_5_S13_L002_R1_001_GTGA.trim.bt2.bam 0.832749464173488
IAN_2brad_L3_6_S6_L001_R1_001_TGTC.trim.bt2.bam 0.832857096090801
IAN_2brad_L4_1_S9_L002_R1_001_TGGT.trim.bt2.bam 0.83319062971572
IAN_2brad_L2_7_S15_L002_R1_001_AGAC.trim.bt2.bam 0.833809249232787
IAN_2brad_L2_4_S12_L002_R1_001_CATC.trim.bt2.bam 0.835077824131168
IAN_2brad_L2_4_S12_L002_R1_001_GCTT.trim.bt2.bam 0.835116089735594
IAN_2brad_L2_5_S13_L002_R1_001_AGTG.trim.bt2.bam 0.835430453862239
IAN_2brad_L2_5_S13_L002_R1_001_TGTC.trim.bt2.bam 0.835659631948136
IAN_2brad_L2_1_S9_L002_R1_001_GTGA.trim.bt2.bam 0.835986952291889
IAN_2brad_L3_3_S3_L001_R1_001_GTGA.trim.bt2.bam 0.836227738258032
IAN_2brad_L3_7_S7_L001_R1_001_CTAC.trim.bt2.bam 0.837310234263963
IAN_2brad_L2_8_S16_L002_R1_001_TGTC.trim.bt2.bam 0.837450298304326
IAN_2brad_L2_2_S10_L002_R1_001_TGTC.trim.bt2.bam 0.837618100939983
IAN_2brad_L2_7_S15_L002_R1_001_CATC.trim.bt2.bam 0.837682536037982
IAN_2brad_L2_7_S15_L002_R1_001_GCTT.trim.bt2.bam 0.837800720169152
IAN_2brad_L3_2_S2_L001_R1_001_TCAC.trim.bt2.bam 0.838320390537017
IAN_2brad_L1_2_S2_L001_R1_001_AGTG.trim.bt2.bam 0.838383330267077
IAN_2brad_L3_8_S8_L001_R1_001_GCTT.trim.bt2.bam 0.838997979643908
IAN_2brad_L1_2_S2_L001_R1_001_TCAC.trim.bt2.bam 0.839046973061432
IAN_2brad_L4_4_S12_L002_R1_001_CATC.trim.bt2.bam 0.839065102180439
IAN_2brad_L1_5_S5_L001_R1_001_AGTG.trim.bt2.bam 0.839894470969795
IAN_2brad_L2_1_S9_L002_R1_001_ACCA.trim.bt2.bam 0.840147592726652
IAN_2brad_L3_6_S6_L001_R1_001_TCAC.trim.bt2.bam 0.840213501079947
IAN_2brad_L4_2_S10_L002_R1_001_GTGA.trim.bt2.bam 0.840520220771317
IAN_2brad_L3_3_S3_L001_R1_001_TGTC.trim.bt2.bam 0.840984943328217
IAN_2brad_L3_1_S1_L001_R1_001_ACCA.trim.bt2.bam 0.841205916991613
IAN_2brad_L4_4_S12_L002_R1_001_AGTG.trim.bt2.bam 0.841517395202823
IAN_2brad_L1_8_S8_L001_R1_001_CATC.trim.bt2.bam 0.842579160181557
IAN_2brad_L2_1_S9_L002_R1_001_CTAC.trim.bt2.bam 0.842820653192173
IAN_2brad_L1_6_S6_L001_R1_001_AGAC.trim.bt2.bam 0.843571314263178
IAN_2brad_L1_4_S4_L001_R1_001_TCAC.trim.bt2.bam 0.843617990770353
IAN_2brad_L4_2_S10_L002_R1_001_ACCA.trim.bt2.bam 0.844840393217958
IAN_2brad_L1_5_S5_L001_R1_001_GTGA.trim.bt2.bam 0.844888873532847
IAN_2brad_L4_3_S11_L002_R1_001_TCAC.trim.bt2.bam 0.845929889764036
IAN_2brad_L2_7_S15_L002_R1_001_GTGA.trim.bt2.bam 0.845956779950013
IAN_2brad_L3_7_S7_L001_R1_001_AGAC.trim.bt2.bam 0.847460236533745
IAN_2brad_L3_2_S2_L001_R1_001_TGTC.trim.bt2.bam 0.847609768293281
IAN_2brad_L3_2_S2_L001_R1_001_ACCA.trim.bt2.bam 0.847860808221749
IAN_2brad_L1_2_S2_L001_R1_001_CTAC.trim.bt2.bam 0.848559061560047
IAN_2brad_L3_4_S4_L001_R1_001_CATC.trim.bt2.bam 0.848650000323687
IAN_2brad_L3_5_S5_L001_R1_001_GTGA.trim.bt2.bam 0.848955146801563
IAN_2brad_L2_6_S14_L002_R1_001_TGGT.trim.bt2.bam 0.849126929558105
IAN_2brad_L3_7_S7_L001_R1_001_TCAC.trim.bt2.bam 0.84964364512412
IAN_2brad_L1_2_S2_L001_R1_001_GTGA.trim.bt2.bam 0.850061983659821
IAN_2brad_L4_1_S9_L002_R1_001_ACCA.trim.bt2.bam 0.850531440294803
IAN_2brad_L3_6_S6_L001_R1_001_GTGA.trim.bt2.bam 0.850737522450599
IAN_2brad_L3_4_S4_L001_R1_001_TCAC.trim.bt2.bam 0.852296667561759
IAN_2brad_L3_5_S5_L001_R1_001_AGAC.trim.bt2.bam 0.852317095357564
IAN_2brad_L1_4_S4_L001_R1_001_ACCA.trim.bt2.bam 0.85298726433251
IAN_2brad_L1_4_S4_L001_R1_001_TGTC.trim.bt2.bam 0.853192209071881
IAN_2brad_L1_5_S5_L001_R1_001_TCAG.trim.bt2.bam 0.853344865989239
IAN_2brad_L1_2_S2_L001_R1_001_TGGT.trim.bt2.bam 0.853414495683906
IAN_2brad_L3_4_S4_L001_R1_001_ACCA.trim.bt2.bam 0.853803618225402
IAN_2brad_L4_4_S12_L002_R1_001_ACCA.trim.bt2.bam 0.853945778358937
IAN_2brad_L4_1_S9_L002_R1_001_GCTT.trim.bt2.bam 0.854395373278307
IAN_2brad_L3_8_S8_L001_R1_001_GACT.trim.bt2.bam 0.854632646340963
IAN_2brad_L2_7_S15_L002_R1_001_TGGT.trim.bt2.bam 0.854809237333268
IAN_2brad_L1_4_S4_L001_R1_001_GTGA.trim.bt2.bam 0.855035905657032
IAN_2brad_L3_7_S7_L001_R1_001_AGTG.trim.bt2.bam 0.855609659080986
IAN_2brad_L2_3_S11_L002_R1_001_TCAC.trim.bt2.bam 0.855746214160661
IAN_2brad_L3_5_S5_L001_R1_001_CTAC.trim.bt2.bam 0.855860700028714
IAN_2brad_L3_5_S5_L001_R1_001_TCAG.trim.bt2.bam 0.857022170212236
IAN_2brad_L1_7_S7_L001_R1_001_CATC.trim.bt2.bam 0.857859959975451
IAN_2brad_L3_6_S6_L001_R1_001_TGGT.trim.bt2.bam 0.858142035066255
IAN_2brad_L1_6_S6_L001_R1_001_GCTT.trim.bt2.bam 0.85851778989276
IAN_2brad_L1_7_S7_L001_R1_001_CTAC.trim.bt2.bam 0.858880304907915
IAN_2brad_L3_2_S2_L001_R1_001_TGGT.trim.bt2.bam 0.859329824779859
IAN_2brad_L2_3_S11_L002_R1_001_TGTC.trim.bt2.bam 0.85939784929641
IAN_2brad_L3_4_S4_L001_R1_001_CTAC.trim.bt2.bam 0.859805775456498
IAN_2brad_L1_2_S2_L001_R1_001_TCAG.trim.bt2.bam 0.860091700434524
IAN_2brad_L1_5_S5_L001_R1_001_TGTC.trim.bt2.bam 0.861115947039031
IAN_2brad_L3_7_S7_L001_R1_001_TGTC.trim.bt2.bam 0.86158067674336
IAN_2brad_L1_7_S7_L001_R1_001_AGAC.trim.bt2.bam 0.861781798501517
IAN_2brad_L2_7_S15_L002_R1_001_AGTG.trim.bt2.bam 0.863851267882443
IAN_2brad_L3_3_S3_L001_R1_001_GACT.trim.bt2.bam 0.863961851394135
IAN_2brad_L3_3_S3_L001_R1_001_AGAC.trim.bt2.bam 0.866191988141522
IAN_2brad_L1_2_S2_L001_R1_001_AGAC.trim.bt2.bam 0.866207960092944
IAN_2brad_L3_1_S1_L001_R1_001_TGGT.trim.bt2.bam 0.866500738163363
IAN_2brad_L3_8_S8_L001_R1_001_TGTC.trim.bt2.bam 0.866624609229207
IAN_2brad_L4_2_S10_L002_R1_001_GCTT.trim.bt2.bam 0.867525702095713
IAN_2brad_L1_5_S5_L001_R1_001_CTAC.trim.bt2.bam 0.868252246019066
IAN_2brad_L1_2_S2_L001_R1_001_ACCA.trim.bt2.bam 0.868434418087071
IAN_2brad_L4_4_S12_L002_R1_001_GCTT.trim.bt2.bam 0.869211556096953
IAN_2brad_L3_5_S5_L001_R1_001_CATC.trim.bt2.bam 0.869345526255664
IAN_2brad_L4_1_S9_L002_R1_001_GACT.trim.bt2.bam 0.869891461573046
IAN_2brad_L2_8_S16_L002_R1_001_CTAC.trim.bt2.bam 0.870308731837973
IAN_2brad_L2_7_S15_L002_R1_001_GACT.trim.bt2.bam 0.871373940886877
IAN_2brad_L3_7_S7_L001_R1_001_GTGA.trim.bt2.bam 0.872351671535506
IAN_2brad_L1_5_S5_L001_R1_001_ACCA.trim.bt2.bam 0.872382321518071
IAN_2brad_L3_7_S7_L001_R1_001_ACCA.trim.bt2.bam 0.872385586460343
IAN_2brad_L2_8_S16_L002_R1_001_TCAC.trim.bt2.bam 0.872528621410213
IAN_2brad_L2_4_S12_L002_R1_001_TGTC.trim.bt2.bam 0.873391275630399
IAN_2brad_L1_5_S5_L001_R1_001_CATC.trim.bt2.bam 0.873406944528517
IAN_2brad_L3_7_S7_L001_R1_001_TGGT.trim.bt2.bam 0.875422017200618
IAN_2brad_L1_8_S8_L001_R1_001_GCTT.trim.bt2.bam 0.8788637165167
IAN_2brad_L3_7_S7_L001_R1_001_CATC.trim.bt2.bam 0.879457664005642
IAN_2brad_L2_4_S12_L002_R1_001_TCAC.trim.bt2.bam 0.879591898445001
IAN_2brad_L3_6_S6_L001_R1_001_CTAC.trim.bt2.bam 0.881002727333486
IAN_2brad_L1_6_S6_L001_R1_001_GTGA.trim.bt2.bam 0.881391669358336
IAN_2brad_L2_8_S16_L002_R1_001_TCAG.trim.bt2.bam 0.88148486133099
IAN_2brad_L1_6_S6_L001_R1_001_TGGT.trim.bt2.bam 0.88221645136032
IAN_2brad_L4_2_S10_L002_R1_001_GACT.trim.bt2.bam 0.882320111281639
IAN_2brad_L3_6_S6_L001_R1_001_AGAC.trim.bt2.bam 0.882566037856187
IAN_2brad_L2_5_S13_L002_R1_001_CTAC.trim.bt2.bam 0.882853625621227
IAN_2brad_L2_3_S11_L002_R1_001_AGAC.trim.bt2.bam 0.88338355516545
IAN_2brad_L4_4_S12_L002_R1_001_GACT.trim.bt2.bam 0.883677781571885
IAN_2brad_L1_6_S6_L001_R1_001_CATC.trim.bt2.bam 0.88375865046856
IAN_2brad_L2_6_S14_L002_R1_001_GACT.trim.bt2.bam 0.88444884816571
IAN_2brad_L2_5_S13_L002_R1_001_TCAG.trim.bt2.bam 0.885284112102071
IAN_2brad_L2_3_S11_L002_R1_001_CTAC.trim.bt2.bam 0.885465090355467
IAN_2brad_L3_1_S1_L001_R1_001_CTAC.trim.bt2.bam 0.88566658119024
IAN_2brad_L1_2_S2_L001_R1_001_CATC.trim.bt2.bam 0.885674406092544
IAN_2brad_L2_3_S11_L002_R1_001_GTGA.trim.bt2.bam 0.885802763342345
IAN_2brad_L1_5_S5_L001_R1_001_AGAC.trim.bt2.bam 0.886270133168454
IAN_2brad_L2_8_S16_L002_R1_001_GTGA.trim.bt2.bam 0.886384399872379
IAN_2brad_L1_6_S6_L001_R1_001_TCAC.trim.bt2.bam 0.886768272741949
IAN_2brad_L2_5_S13_L002_R1_001_TCAC.trim.bt2.bam 0.888485085961348
IAN_2brad_L2_3_S11_L002_R1_001_TGGT.trim.bt2.bam 0.888795501292728
IAN_2brad_L3_1_S1_L001_R1_001_TCAG.trim.bt2.bam 0.888832999941921
IAN_2brad_L3_6_S6_L001_R1_001_GACT.trim.bt2.bam 0.889730162989513
IAN_2brad_L1_3_S3_L001_R1_001_AGAC.trim.bt2.bam 0.890114329651949
IAN_2brad_L3_5_S5_L001_R1_001_TGGT.trim.bt2.bam 0.890764598323234
IAN_2brad_L2_8_S16_L002_R1_001_GACT.trim.bt2.bam 0.891100568561606
IAN_2brad_L3_8_S8_L001_R1_001_CTAC.trim.bt2.bam 0.891426892081747
IAN_2brad_L2_4_S12_L002_R1_001_GTGA.trim.bt2.bam 0.891653800519051
IAN_2brad_L1_4_S4_L001_R1_001_CTAC.trim.bt2.bam 0.892128160600378
IAN_2brad_L1_6_S6_L001_R1_001_GACT.trim.bt2.bam 0.892928521191008
IAN_2brad_L1_8_S8_L001_R1_001_TCAG.trim.bt2.bam 0.893472311278226
IAN_2brad_L4_3_S11_L002_R1_001_TGGT.trim.bt2.bam 0.894408048716826
IAN_2brad_L4_2_S10_L002_R1_001_TCAG.trim.bt2.bam 0.894447085373808
IAN_2brad_L4_1_S9_L002_R1_001_AGAC.trim.bt2.bam 0.894878326780434
IAN_2brad_L1_4_S4_L001_R1_001_CATC.trim.bt2.bam 0.895975856922867
IAN_2brad_L2_3_S11_L002_R1_001_GACT.trim.bt2.bam 0.89606541630444
IAN_2brad_L2_7_S15_L002_R1_001_TCAG.trim.bt2.bam 0.896120058292067
IAN_2brad_L2_5_S13_L002_R1_001_TGGT.trim.bt2.bam 0.896297371916647
IAN_2brad_L3_5_S5_L001_R1_001_GACT.trim.bt2.bam 0.896441867405172
IAN_2brad_L3_1_S1_L001_R1_001_TCAC.trim.bt2.bam 0.896847767925839
IAN_2brad_L4_4_S12_L002_R1_001_AGAC.trim.bt2.bam 0.897353955106795
IAN_2brad_L4_2_S10_L002_R1_001_TCAC.trim.bt2.bam 0.898855519728835
IAN_2brad_L4_2_S10_L002_R1_001_CATC.trim.bt2.bam 0.9006180584626
IAN_2brad_L2_5_S13_L002_R1_001_AGAC.trim.bt2.bam 0.901706375606679
IAN_2brad_L2_6_S14_L002_R1_001_TCAC.trim.bt2.bam 0.901954094126074
IAN_2brad_L2_8_S16_L002_R1_001_CATC.trim.bt2.bam 0.902796402987024
IAN_2brad_L4_4_S12_L002_R1_001_TCAC.trim.bt2.bam 0.90318331444054
IAN_2brad_L4_3_S11_L002_R1_001_CTAC.trim.bt2.bam 0.903402429226111
IAN_2brad_L4_3_S11_L002_R1_001_TGTC.trim.bt2.bam 0.903449501683172
IAN_2brad_L2_3_S11_L002_R1_001_CATC.trim.bt2.bam 0.904307409040996
IAN_2brad_L2_6_S14_L002_R1_001_CATC.trim.bt2.bam 0.905382220898098
IAN_2brad_L3_5_S5_L001_R1_001_TCAC.trim.bt2.bam 0.905498352666172
IAN_2brad_L3_1_S1_L001_R1_001_GACT.trim.bt2.bam 0.906223778931315
IAN_2brad_L2_7_S15_L002_R1_001_CTAC.trim.bt2.bam 0.908059032425832
IAN_2brad_L2_2_S10_L002_R1_001_TCAC.trim.bt2.bam 0.908160872856419
IAN_2brad_L4_4_S12_L002_R1_001_CTAC.trim.bt2.bam 0.908500329702256
IAN_2brad_L4_1_S9_L002_R1_001_GTGA.trim.bt2.bam 0.909061330870176
IAN_2brad_L2_7_S15_L002_R1_001_TCAC.trim.bt2.bam 0.911129523815861
IAN_2brad_L4_3_S11_L002_R1_001_CATC.trim.bt2.bam 0.911345566596929
IAN_2brad_L2_5_S13_L002_R1_001_CATC.trim.bt2.bam 0.912280001735468
IAN_2brad_L4_1_S9_L002_R1_001_CATC.trim.bt2.bam 0.912557871866605
IAN_2brad_L4_3_S11_L002_R1_001_TCAG.trim.bt2.bam 0.91356173393828
IAN_2brad_L2_4_S12_L002_R1_001_AGAC.trim.bt2.bam 0.91457409061086
IAN_2brad_L2_6_S14_L002_R1_001_TCAG.trim.bt2.bam 0.914654704094381
IAN_2brad_L2_7_S15_L002_R1_001_TGTC.trim.bt2.bam 0.915812114097841
IAN_2brad_L1_3_S3_L001_R1_001_TGGT.trim.bt2.bam 0.917192859694387
IAN_2brad_L4_2_S10_L002_R1_001_CTAC.trim.bt2.bam 0.919097743760604
IAN_2brad_L1_3_S3_L001_R1_001_ACCA.trim.bt2.bam 0.922202105709176
IAN_2brad_L2_4_S12_L002_R1_001_GACT.trim.bt2.bam 0.926040696804029
IAN_2brad_L1_8_S8_L001_R1_001_TCAC.trim.bt2.bam 0.926567571933589
IAN_2brad_L4_2_S10_L002_R1_001_TGGT.trim.bt2.bam 0.926912677093424


# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds, and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ,  -minIndDepth and -minInd filters in subsequent ANGSD runs
# we have a few samples with very low coverage, but mostly looks like high quality sequencing of what we do have. 

#I just use cyberduck to get dd.pdf onto my local machine into my Chapter_2_RAD folder

pwd /projectnb/mullenl/novick/global_pop_gen/rad_workflow/ANGSD

#which will move the pdf onto your local computer into your current directory. 

#----------- clones detection (**minInd ~80% of samples)

# the -doBcf flag produces a vcf file at this stage.
# set minInd to 80% of number of samples --> 260

# Making a new qsub for this run called angsd2.qsub

#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ANGSD1 # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M inovick@bu.edu
#$ -m be

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 260 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"


angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out myresult2

NSITES=`zcat myresult2.mafs.gz | wc -l`
echo $NSITES

# 6145

# Next change mindepthind to 5, minind 80%

# Add the -setMinDepthInd filter to see what happens, this filter was really important in Nicola's project but not included in Misha's pipeline. 
# Starting with -setMinDepthInd 2 and see what happens when we increase it

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 260 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 2"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out myresult3

NSITES=`zcat myresult3.mafs.gz | wc -l`
echo $NSITES
# 5267

# Next change mindepthind to 5, minind 80%

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 260 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 5"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out myresult4

NSITES=`zcat myresult4.mafs.gz | wc -l`
echo $NSITES
# 4193

##For the above run, I accidentally overwrote all the previous myresult3 runs because I forgot to change the name. So this is all myresult3 now, everything from the previous run is gone except the ibsMat that I scp to my local computer



#------------------------------ remove clones and re-run angsd

#Run 3 analyses, 1 where we remove all the bads/ everything below 50%, one where we remove everyone below 10%, then one where we keep everything. Remove tech reps for all, do minindepth 3 for all, 80% minind
# To choose which replicate to use, you want to choose the one with the most number of reads

##Deciding which replicate to use, which to get rid of, and removing those with a proportion of sites covered at 5x under 50%:

#Reads with a proportion of sites covered at 5x under 50%:

IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam 0.0110273010437435   -> IAN6-2-22-4 Kwara state, Nigeria #remove
IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam 0.0193986356835845   -> IAN6-2-22-3 Kwara state, Nigeria #remove
IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam 0.0416131656037242   -> IAN6-2-22-4 Kwara state, Nigeria #keep*
IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam 0.0508641375683166   -> IAN6-2-22-3 Kwara state, Nigeria #keep*
IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam 0.0993919746965222    -> IAN3-17-22-14 Belmont, MA    #remove
IAN_2brad_L1_7_S7_L001_R1_001_GACT.trim.bt2.bam 0.11513661706499      -> IAN7-18-22-12 Bethlehem, Palestine #remove?
IAN_2brad_L1_4_S4_L001_R1_001_AGAC.trim.bt2.bam 0.124160484766929     -> IAN6-2-22-1 Kwara state, Nigeria #keep
IAN_2brad_L1_7_S7_L001_R1_001_TCAC.trim.bt2.bam 0.276104952113743     -> IAN7-18-22-11 Bethlehem, Palestine #keep
IAN_2brad_L1_5_S5_L001_R1_001_TGGT.trim.bt2.bam 0.441787876477245     -> IAN5-19-22-25 Adelaide, Australia   #remove

#Then use this script to count total number of reads in bam file:
module load samtools
samtools view -c SAMPLE.bam

#List of tech reps:												 #no. of reads			#location

IAN_2brad_L1_3_S3_L001_R1_001_TGGT.trim.bt2.bam -- IAN6-13-22-35 #1917668 #keep			#Albertinaplatz, Vienna
IAN_2brad_L3_1_S1_L001_R1_001_TGGT.trim.bt2.bam -- IAN6-13-22-35 #1106907 #remove		#Albertinaplatz, Vienna

IAN_2brad_L3_3_S3_L001_R1_001_AGAC.trim.bt2.bam -- IAN6-13-22-83 #1391998 #remove		#Albertinaplatz, Vienna
IAN_2brad_L1_3_S3_L001_R1_001_AGAC.trim.bt2.bam -- IAN6-13-22-83 #2029866 #keep			#Albertinaplatz, Vienna

IAN_2brad_L3_4_S4_L001_R1_001_ACCA.trim.bt2.bam -- IAN6-20-22-26 #1075228 #remove		#Mariahilf, Vienna
IAN_2brad_L1_3_S3_L001_R1_001_ACCA.trim.bt2.bam -- IAN6-20-22-26 #2187474 #keep			#Mariahilf, Vienna

IAN_2brad_L1_2_S2_L001_R1_001_AGAC.trim.bt2.bam -- IAN3-17-22-26 #1927691 #keep			#Belmont, MA
IAN_2brad_L1_8_S8_L001_R1_001_ACCA.trim.bt2.bam -- IAN3-17-22-26 #637744 #remove		#Belmont, MA

IAN_2brad_L1_8_S8_L001_R1_001_GTGA.trim.bt2.bam -- IAN3-17-22-29 #1259273 #remove		#Belmont, MA
IAN_2brad_L1_2_S2_L001_R1_001_ACCA.trim.bt2.bam -- IAN3-17-22-29 #1763040 #keep			#Belmont, MA

IAN_2brad_L1_2_S2_L001_R1_001_AGTG.trim.bt2.bam -- IAN3-17-22-34 #1903828 #keep			#Belmont, MA
IAN_2brad_L1_8_S8_L001_R1_001_CTAC.trim.bt2.bam -- IAN3-17-22-34 #868540 #remove		#Belmont, MA

IAN_2brad_L1_8_S8_L001_R1_001_TGTC.trim.bt2.bam -- IAN11-22-21-13 #1049390 #remove		#Cambridge, MA
IAN_2brad_L1_5_S5_L001_R1_001_CATC.trim.bt2.bam -- IAN11-22-21-13 #1807356 #keep		#Cambridge, MA

IAN_2brad_L1_5_S5_L001_R1_001_GTGA.trim.bt2.bam -- IAN11-22-21-16 #1865569 #keep		#Cambridge, MA	
IAN_2brad_L1_8_S8_L001_R1_001_GACT.trim.bt2.bam -- IAN11-22-21-16 #755920 #remove		#Cambridge, MA

IAN_2brad_L2_4_S12_L002_R1_001_ACCA.trim.bt2.bam -- IAN3-24-22-9 #646561 #keep			#Boston, MA
IAN_2brad_L2_5_S13_L002_R1_001_ACCA.trim.bt2.bam -- IAN3-24-22-9 #593865 #remove		#Boston, MA

IAN_2brad_L2_6_S14_L002_R1_001_TCAG.trim.bt2.bam -- IAN6-13-22-25 #2938140 #keep		#Albertinaplatz, Vienna
IAN_2brad_L2_8_S16_L002_R1_001_TCAG.trim.bt2.bam -- IAN6-13-22-25 #2008032 #remove		#Albertinaplatz, Vienna
IAN_2brad_L3_1_S1_L001_R1_001_TCAG.trim.bt2.bam -- IAN6-13-22-25 #1669833 #remove		#Albertinaplatz, Vienna

IAN_2brad_L2_6_S14_L002_R1_001_ACCA.trim.bt2.bam -- IAN5-28-22-20 #829202 #remove		#Richmond, UK
IAN_2brad_L2_7_S15_L002_R1_001_ACCA.trim.bt2.bam -- IAN5-28-22-20 #843613 #keep			#Richmond, UK

IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam -- IAN6-2-22-3 #76558 #keep			#Kwara State, Nigeria
IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam -- IAN6-2-22-3 #55854 #remove			#Kwara State, Nigeria

IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam -- IAN6-2-22-4 #60245 #keep			#Kwara State, Nigeria
IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam -- IAN6-2-22-4 #41847 #remove			#Kwara State, Nigeria

IAN_2brad_L4_3_S11_L002_R1_001_AGTG.trim.bt2.bam -- IAN8-18-22-1 #1415677 #keep			#Canberra, Australia
IAN_2brad_L2_8_S16_L002_R1_001_AGTG.trim.bt2.bam -- IAN8-18-22-1 #1382630 #remove		#Canberra, Australia
IAN_2brad_L3_8_S8_L001_R1_001_GACT.trim.bt2.bam -- IAN8-18-22-1 #1217871 #remove		#Canberra, Australia
IAN_2brad_L4_1_S9_L002_R1_001_TGGT.trim.bt2.bam -- IAN8-18-22-1 #1256891 #remove		#Canberra, Australia

IAN_2brad_L4_3_S11_L002_R1_001_CATC.trim.bt2.bam -- IAN8-18-22-4 #2794577 #keep			#Canberra, Australia
IAN_2brad_L2_8_S16_L002_R1_001_CATC.trim.bt2.bam -- IAN8-18-22-4 #2696478 #remove		#Canberra, Australia

IAN_2brad_L3_1_S1_L001_R1_001_CATC.trim.bt2.bam -- IAN6-13-22-44 #646763 #keep			#Albertinaplatz, Vienna	
IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim.bt2.bam -- IAN6-13-22-44 #546451 #remove		#Albertinaplatz, Vienna

IAN_2brad_L3_4_S4_L001_R1_001_GCTT.trim.bt2.bam -- IAN6-20-22-45 #1261414 #keep			#Mariahilf, Vienna
IAN_2brad_L3_5_S5_L001_R1_001_GCTT.trim.bt2.bam -- IAN6-20-22-45 #1221530 #remove		#Mariahilf, Vienna

IAN_2brad_L3_6_S6_L001_R1_001_CATC.trim.bt2.bam -- IAN11-22-21-08 #394493 #remove		#Cambridge, MA
IAN_2brad_L1_5_S5_L001_R1_001_ACCA.trim.bt2.bam -- IAN11-22-21-08 #1969009 #keep		#Cambridge, MA

IAN_2brad_L3_7_S7_L001_R1_001_CTAC.trim.bt2.bam -- IAN8-18-22-7 #1035518 #remove		#Canberra, Australia
IAN_2brad_L3_8_S8_L001_R1_001_CTAC.trim.bt2.bam -- IAN8-18-22-7 #1761498 #keep			#Canberra, Australia

IAN_2brad_L2_3_S11_L002_R1_001_TCAG.trim.bt2.bam -- IAN2-1-22-28 #1579647 #keep			#Cambridge, MA
IAN_2brad_L3_7_S7_L001_R1_001_GACT.trim.bt2.bam -- IAN2-1-22-28 #1153178 #remove		#Cambridge, MA

IAN_2brad_L3_8_S8_L001_R1_001_TGTC.trim.bt2.bam -- IAN5-19-22-18 #1732258 #remove		#Adelaide, Australia		#Can't find in dataset?
IAN_2brad_L4_2_S10_L002_R1_001_TCAC.trim.bt2.bam -- IAN5-19-22-18 #2451836 #remove		#Adelaide, Australia
IAN_2brad_L4_4_S12_L002_R1_001_TCAC.trim.bt2.bam -- IAN5-19-22-18 #2681971 #keep		#Adelaide, Australia

IAN_2brad_L4_1_S9_L002_R1_001_CATC.trim.bt2.bam -- IAN8-23-22-15 #1895564 #keep			#Albuquerque, NM
IAN_2brad_L4_2_S10_L002_R1_001_CATC.trim.bt2.bam -- IAN8-23-22-15 #1796299 #remove		#Albuquerque, NM

IAN_2brad_L4_2_S10_L002_R1_001_TCAG.trim.bt2.bam -- IAN3-17-22-8 #2875644 #keep			#Belmont, MA
IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim.bt2.bam -- IAN3-17-22-8 #1684653 #remove		#Belmont, MA

IAN_2brad_L4_2_S10_L002_R1_001_TGGT.trim.bt2.bam -- IAN3-22-22-7 #2371379 #keep			#Boston, MA
IAN_2brad_L4_3_S11_L002_R1_001_TGGT.trim.bt2.bam -- IAN3-22-22-7 #1726486 #remove		#Boston, MA

IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim.bt2.bam -- IAN3-17-22-20 #1697793 #remove		#Belmont, MA
IAN_2brad_L4_3_S11_L002_R1_001_CTAC.trim.bt2.bam -- IAN3-17-22-20 #2041673 #keep		#Belmont, MA

IAN_2brad_L4_3_S11_L002_R1_001_TGTC.trim.bt2.bam -- IAN6-13-22-82 #2657134 #keep		#Albertinaplatz, Vienna	
IAN_2brad_L1_4_S4_L001_R1_001_GTGA.trim.bt2.bam -- IAN6-13-22-82 #2027699 #remove		#Albertinaplatz, Vienna	

IAN_2brad_L1_6_S6_L001_R1_001_AGAC.trim.bt2.bam -- IAN3-31-22-3 #1608716 #remove		#Suitland, MD
IAN_2brad_L4_3_S11_L002_R1_001_TCAC.trim.bt2.bam -- IAN3-31-22-3 #2402001 #keep			#Suitland, MD

IAN_2brad_L4_2_S10_L002_R1_001_GCTT.trim.bt2.bam -- IAN5-19-22-8 #2262861 #remove		#Adelaide, Australia
IAN_2brad_L4_4_S12_L002_R1_001_GCTT.trim.bt2.bam -- IAN5-19-22-8 #2393476 #keep			#Adelaide, Australia

IAN_2brad_L4_2_S10_L002_R1_001_CTAC.trim.bt2.bam -- IAN5-19-22-13 #1698417 #remove		#Adelaide, Australia
IAN_2brad_L4_4_S12_L002_R1_001_CTAC.trim.bt2.bam -- IAN5-19-22-13 #1782172 #keep		#Adelaide, Australia

IAN_2brad_L4_4_S12_L002_R1_001_TGTC.trim.bt2.bam -- IAN5-19-22-17 #1636926 #remove		#Adelaide, Australia
IAN_2brad_L4_2_S10_L002_R1_001_TGTC.trim.bt2.bam -- IAN5-19-22-17 #1660442 #keep		#Adelaide, Australia

IAN_2brad_L4_4_S12_L002_R1_001_GACT.trim.bt2.bam -- IAN5-19-22-19 #2096125 #keep		#Adelaide, Australia
IAN_2brad_L4_2_S10_L002_R1_001_GACT.trim.bt2.bam -- IAN5-19-22-19 #1767896 #remove		#Adelaide, Australia

#### new dups that I forgot!!!! GO BACK REDO EVERYTHING AND ADD THESE

IAN_2brad_L4_3_S11_L002_R1_001_GCTT.trim.bt2.bam -- IAN3-31-22-4 #1911856 #keep			#Suitland, MD
IAN_2brad_L1_2_S2_L001_R1_001_CTAC.trim.bt2.bam -- IAN3-31-22-4 #1708109 #remove		#Suitland, MD

IAN_2brad_L3_5_S5_L001_R1_001_GACT.trim.bt2.bam -- IAN2-1-22-6 #2072796 #keep			#Cambridge, MA
IAN_2brad_L1_5_S5_L001_R1_001_CTAC.trim.bt2.bam -- IAN2-1-22-06 #1414710 #remove		#Cambridge, MA

# List of bams to remove for run removing everything under 50% coverage and TECH REPS. This includes all Nigeria and all Bethlehem-- alldups_50p_bams_nodups
IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam 
IAN_2brad_L1_7_S7_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L1_4_S4_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L1_7_S7_L001_R1_001_TCAC.trim.bt2.bam 
IAN_2brad_L1_5_S5_L001_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L3_1_S1_L001_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L3_3_S3_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L3_4_S4_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L2_5_S13_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L3_1_S1_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L2_6_S14_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L4_1_S9_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_5_S5_L001_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_CATC.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L1_4_S4_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_6_S6_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L4_4_S12_L002_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_GACT.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_CTAC.trim.bt2.bam


# Making a new bams file "alldups_50p_bams_nodups" with all above bams removed:
grep -v -e IAN_2brad_L1_2_S2_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_5_S5_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_7_S7_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L1_4_S4_L001_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L1_7_S7_L001_R1_001_TCAC.trim.bt2.bam  -e IAN_2brad_L1_5_S5_L001_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L3_1_S1_L001_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L3_3_S3_L001_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L3_4_S4_L001_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_GTGA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L2_5_S13_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L3_1_S1_L001_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L2_6_S14_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_AGTG.trim.bt2.bam -e IAN_2brad_L3_8_S8_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L4_1_S9_L002_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_5_S5_L001_R1_001_GCTT.trim.bt2.bam -e IAN_2brad_L3_6_S6_L001_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_7_S7_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L3_7_S7_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L3_8_S8_L001_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_TCAC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim.bt2.bam -e IAN_2brad_L1_4_S4_L001_R1_001_GTGA.trim.bt2.bam -e IAN_2brad_L1_6_S6_L001_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_GCTT.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L4_4_S12_L002_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_GACT.trim.bt2.bam bams > alldups_50p_bams_nodups

# List of bams to remove for final 50p run, keeping tech reps in, but removing everything under 50% coverage at 5x
IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_7_S7_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L1_4_S4_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L1_7_S7_L001_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_TGGT.trim.bt2.bam

# Making a new bams file for this, 50p_final_with_techreps
grep -v -e IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_7_S7_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L1_4_S4_L001_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L1_7_S7_L001_R1_001_TCAC.trim.bt2.bam -e IAN_2brad_L1_5_S5_L001_R1_001_TGGT.trim.bt2.bam bams > 50p_final_with_techreps

# List of bams to remove for final 12p run, keeping tech reps in, but removing everything under 12% coverage at 5x
IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_7_S7_L001_R1_001_GACT.trim.bt2.bam

# Making new bams file for this, 12p_final_with_techreps
grep -v -e IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_7_S7_L001_R1_001_GACT.trim.bt2.bam bams > 12p_final_with_techreps

# List of bams to remove for run removing everything under 10% coverage and TECH REPS. This includes all Nigeria except the one above 10%-- alldups_10p_bams_nodups
IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam 
IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L3_1_S1_L001_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L3_3_S3_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L3_4_S4_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L2_5_S13_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L3_1_S1_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L2_6_S14_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L4_1_S9_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_5_S5_L001_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_CATC.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L1_4_S4_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_6_S6_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L4_4_S12_L002_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_GACT.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_CTAC.trim.bt2.bam

# Making a new bams file "alldups_10p_bams_nodups" with all above bams removed:
grep -v -e IAN_2brad_L1_2_S2_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_5_S5_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L3_1_S1_L001_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L3_3_S3_L001_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L3_4_S4_L001_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_GTGA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L2_5_S13_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L3_1_S1_L001_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L2_6_S14_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_AGTG.trim.bt2.bam -e IAN_2brad_L3_8_S8_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L4_1_S9_L002_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_5_S5_L001_R1_001_GCTT.trim.bt2.bam -e IAN_2brad_L3_6_S6_L001_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_7_S7_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L3_7_S7_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L3_8_S8_L001_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_TCAC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim.bt2.bam -e IAN_2brad_L1_4_S4_L001_R1_001_GTGA.trim.bt2.bam -e IAN_2brad_L1_6_S6_L001_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_GCTT.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L4_4_S12_L002_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_GACT.trim.bt2.bam bams > alldups_10p_bams_nodups


# List of bams to remove for run removing only tech reps, including everything no matter how poor the coverage -- alldups_allbams_nodups
IAN_2brad_L3_1_S1_L001_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L3_3_S3_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L3_4_S4_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L2_5_S13_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L3_1_S1_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L2_6_S14_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L4_1_S9_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_5_S5_L001_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_CATC.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L1_4_S4_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_6_S6_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L4_4_S12_L002_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_GACT.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_CTAC.trim.bt2.bam

# Making a new bams file "alldups_allbams_nodups" with all above bams removed:
grep -v -e IAN_2brad_L1_2_S2_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_5_S5_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L3_1_S1_L001_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L3_3_S3_L001_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L3_4_S4_L001_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_GTGA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L2_5_S13_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L3_1_S1_L001_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L2_6_S14_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_AGTG.trim.bt2.bam -e IAN_2brad_L3_8_S8_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L4_1_S9_L002_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_5_S5_L001_R1_001_GCTT.trim.bt2.bam -e IAN_2brad_L3_6_S6_L001_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_7_S7_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L3_7_S7_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L3_8_S8_L001_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_TCAC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim.bt2.bam -e IAN_2brad_L1_4_S4_L001_R1_001_GTGA.trim.bt2.bam -e IAN_2brad_L1_6_S6_L001_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_GCTT.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L4_4_S12_L002_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_GACT.trim.bt2.bam bams > alldups_allbams_nodups

#I'm also going to make a file "alldups_12p_nodups" removing dups and everything below 12%, so I'm just getting rid of that weird palestine one and keeping the good Palestine one and the best Nigeria one

IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam 
IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L3_1_S1_L001_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L3_3_S3_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L3_4_S4_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L2_5_S13_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L3_1_S1_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L2_6_S14_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L4_1_S9_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L2_8_S16_L002_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_5_S5_L001_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_CATC.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L1_4_S4_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_6_S6_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L4_4_S12_L002_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_GACT.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_7_S7_L001_R1_001_GACT.trim.bt2.bam

# Making a new bams file "alldups_12p_bams_nodups" with all above bams removed:
grep -v -e IAN_2brad_L1_7_S7_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L1_2_S2_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_5_S5_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L3_1_S1_L001_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L3_3_S3_L001_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L3_4_S4_L001_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_GTGA.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L1_8_S8_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L2_5_S13_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L3_1_S1_L001_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L2_6_S14_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_ACCA.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_AGTG.trim.bt2.bam -e IAN_2brad_L3_8_S8_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L4_1_S9_L002_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L2_8_S16_L002_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_2_S2_L001_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_5_S5_L001_R1_001_GCTT.trim.bt2.bam -e IAN_2brad_L3_6_S6_L001_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L3_7_S7_L001_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L3_7_S7_L001_R1_001_GACT.trim.bt2.bam -e IAN_2brad_L3_8_S8_L001_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_TCAC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_CATC.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim.bt2.bam -e IAN_2brad_L4_3_S11_L002_R1_001_TGGT.trim.bt2.bam -e IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim.bt2.bam -e IAN_2brad_L1_4_S4_L001_R1_001_GTGA.trim.bt2.bam -e IAN_2brad_L1_6_S6_L001_R1_001_AGAC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_GCTT.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_CTAC.trim.bt2.bam -e IAN_2brad_L4_4_S12_L002_R1_001_TGTC.trim.bt2.bam -e IAN_2brad_L4_2_S10_L002_R1_001_GACT.trim.bt2.bam bams > alldups_12p_bams_nodups


# All of this is being done in this directory:
pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/ANGSD

# Also just moved all actual bam files to new directory 
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/ANGSD/bams_dir

# Let's count how many bams are in each to make sure we got them all
cat alldups_50p_bams_nodups | wc -l
#281

cat alldups_10p_bams_nodups | wc -l
#285

cat alldups_allbams_nodups | wc -l
#288

cat alldups_12p_bams_nodups | wc -l
#284

####--- new runs ----####
cat 50p_final_with_techreps | wc -l

#316

# All set

#----------------------------- Re-run angsd with all samples but without duplicates
# First doing alldups_50p_bams_nodups: Removed clones and anything under 50% with 5x coverage, including all Palestine and Nigerian samples

pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/ANGSD

# Made new qsub here called 50p_nodups.qsub
# Using minind 80%, minindepth 3

module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 260 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b alldups_50p_bams_nodups -GL 1 $FILTERS $TODO -P 1 -out myresult_alldups_50p_nodups

NSITES=`zcat myresult_alldups_50p_nodups.mafs.gz | wc -l`
echo $NSITES
# 3917

# Now doing alldups_10p_bams_nodups: Removed clones and anything under 10% with 5x coverage, including most of the Nigeria samples. 1 Nigeria sample and both Palestine samples remain in this run

# Made new qsub here called 10p_nodups.qsub
# Using minind 80%, minindepth 3

module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 260 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b alldups_10p_bams_nodups -GL 1 $FILTERS $TODO -P 1 -out myresult_alldups_10p_nodups

NSITES=`zcat myresult_alldups_10p_nodups.mafs.gz | wc -l`
echo $NSITES
# 3995

# Now doing alldups_allbams_nodups, which includes everything except the clones 

# Made new qsub here called allbams_nodups.qsub
# Using minind 80%, minindepth 3

module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 260 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b alldups_allbams_nodups -GL 1 $FILTERS $TODO -P 1 -out myresult_alldups_allbams_nodups

NSITES=`zcat myresult_alldups_allbams_nodups.mafs.gz | wc -l`
echo $NSITES
# 4016

# Now doing alldups_12p_bams_nodups: Removed clones and anything under 12% with 5x coverage, including most of the Nigeria samples. 1 Nigeria sample and one good Palestine sample remain in this run

# Made new qsub here called 12p_nodups.qsub
# Using minind 80%, minindepth 3
# This is actually technically a minind of 91.5% (260/284) because of my mistake, but I like how the dendrogram comes out so I think I will use it for future analyses


module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 260 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b alldups_12p_bams_nodups -GL 1 $FILTERS $TODO -P 1 -out myresult_alldups_12p_nodups

NSITES=`zcat myresult_alldups_12p_nodups.mafs.gz | wc -l`
echo $NSITES
# 3984

# Realized I made a mistake here and didn't change the mindInd to actually 80% of 284. I made that mistake with all of the above ones as well omgggg kill me
# This run is going to be v2_alldups_12p_bams_nodups: Removed clones and anything under 12% with 5x coverage, including most of the Nigeria samples. 1 Nigeria sample and one good Palestine sample remain in this run

# Made new qsub here called v2_12p_nodups.qsub
# Using minind 80%, minindepth 3

module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 227 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"



angsd -b alldups_12p_bams_nodups -GL 1 $FILTERS $TODO -P 1 -out v2_myresult_alldups_12p_nodups

NSITES=`zcat v2_myresult_alldups_12p_nodups.mafs.gz | wc -l`
echo $NSITES
# 5109

############----------- New Stuff (Final Analyses) ---------------#############################################

# Prefixes for final global files are:
final_50p_nodups
final_50p_withdups
final_50_nodups_fis
final_12p_nodups
final_12p_withdups


### Now doing a run for later analyses (2/28/25) where keeping everything, and doing 80% minind and mindindepth of 5 and KEEPING TECHNICAL REPLICATES

cd /projectnb/mullenl/novick/global_pop_gen/rad_workflow/ANGSD

### Removing everything under 50% at 5x coverage, no tech reps, 80% minind, mindepthind 3
#50p_nodups.qsub

module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 225 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b alldups_50p_bams_nodups -GL 1 $FILTERS $TODO -P 1 -out final_50p_nodups

NSITES=`zcat final_50p_nodups.mafs.gz | wc -l`
echo $NSITES
# 5134

### Removing everything under 50% at 5x coverage, KEEPING TECH REPS, 80% minind, mindepthind 3
#50p_withdups.qsub

module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 316 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b 50p_final_with_techreps -GL 1 $FILTERS $TODO -P 1 -out final_50p_withdups

NSITES=`zcat final_50p_withdups.mafs.gz | wc -l`
echo $NSITES
# 1551

### Removing everything under 50% at 5x coverage, no tech reps, 80% minind, mindepthind 3, hard called genotypes for Fis. Changed doGeno to 2 for hard calling, postCutoff 0.95 for confidence in hard calls
# qsub called 50p_nodups_FIS.qsub 
module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 225 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 2 -doBcf 1 -doPost 1 -postCutoff 0.95 -doGlf 2"

angsd -b alldups_50p_bams_nodups -GL 1 $FILTERS $TODO -P 1 -out final_50p_nodups_fis

NSITES=`zcat final_50p_nodups_fis.mafs.gz | wc -l`
echo $NSITES
# 4598
# For FIS, going to use the final_50p_nodups_fis.geno.gz file

### Removing everything under 12% at 5x coverage, removing tech reps, 80% minind, mindepthind 3
# qsub is 12p_nodups.qsub
module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 227 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b alldups_12p_bams_nodups -GL 1 $FILTERS $TODO -P 1 -out final_12p_nodups

NSITES=`zcat final_12p_nodups.mafs.gz | wc -l`
echo $NSITES
# 5109

### Removing everything under 12% at 5x coverage, KEEPING TECH REPS, 80% minind, mindepthind 3
# qsub is 12p_withdups.qsub
module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 255 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b 12p_final_with_techreps -GL 1 $FILTERS $TODO -P 1 -out final_12p_withdups

NSITES=`zcat final_12p_withdups.mafs.gz | wc -l`
echo $NSITES
# 4922

######## Going to do one more ANGSD run of only guys who are from Massachusetts
#Copy all MA bams from alldups_allbams_nodups into new file mass_bams

pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/ANGSD

#Mass bams list:
IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L1_2_S2_L001_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L1_4_S4_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L1_4_S4_L001_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L2_1_S9_L002_R1_001_GACT.trim.bt2.bam
IAN_2brad_L2_1_S9_L002_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L2_1_S9_L002_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L2_2_S10_L002_R1_001_GACT.trim.bt2.bam
IAN_2brad_L2_2_S10_L002_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L2_3_S11_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L2_3_S11_L002_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L2_3_S11_L002_R1_001_CATC.trim.bt2.bam
IAN_2brad_L2_3_S11_L002_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L2_3_S11_L002_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L2_3_S11_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L2_4_S12_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L2_4_S12_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L2_4_S12_L002_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L2_4_S12_L002_R1_001_CATC.trim.bt2.bam
IAN_2brad_L2_4_S12_L002_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L2_4_S12_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L3_5_S5_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_GCTT.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L3_8_S8_L001_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_AGAC.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_AGTG.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L4_2_S10_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_CTAC.trim.bt2.bam

# Make new file for them called mass_bams
touch mass_bams

# Copy them over

ls IAN_2brad_L1_1_S1_L001_R1_001_ACCA.trim.bt2.bam \
IAN_2brad_L1_1_S1_L001_R1_001_AGAC.trim.bt2.bam \
IAN_2brad_L1_1_S1_L001_R1_001_AGTG.trim.bt2.bam \
IAN_2brad_L1_1_S1_L001_R1_001_CATC.trim.bt2.bam \
IAN_2brad_L1_1_S1_L001_R1_001_CTAC.trim.bt2.bam \
IAN_2brad_L1_1_S1_L001_R1_001_GACT.trim.bt2.bam \
IAN_2brad_L1_1_S1_L001_R1_001_GCTT.trim.bt2.bam \
IAN_2brad_L1_1_S1_L001_R1_001_GTGA.trim.bt2.bam \
IAN_2brad_L1_1_S1_L001_R1_001_TGGT.trim.bt2.bam \
IAN_2brad_L1_1_S1_L001_R1_001_TGTC.trim.bt2.bam \
IAN_2brad_L1_2_S2_L001_R1_001_ACCA.trim.bt2.bam \
IAN_2brad_L1_2_S2_L001_R1_001_AGAC.trim.bt2.bam \
IAN_2brad_L1_2_S2_L001_R1_001_AGTG.trim.bt2.bam \
IAN_2brad_L1_2_S2_L001_R1_001_CATC.trim.bt2.bam \
IAN_2brad_L1_2_S2_L001_R1_001_GCTT.trim.bt2.bam \
IAN_2brad_L1_2_S2_L001_R1_001_GTGA.trim.bt2.bam \
IAN_2brad_L1_2_S2_L001_R1_001_TCAG.trim.bt2.bam \
IAN_2brad_L1_2_S2_L001_R1_001_TGGT.trim.bt2.bam \
IAN_2brad_L1_4_S4_L001_R1_001_GACT.trim.bt2.bam \
IAN_2brad_L1_4_S4_L001_R1_001_TCAC.trim.bt2.bam \
IAN_2brad_L1_5_S5_L001_R1_001_ACCA.trim.bt2.bam \
IAN_2brad_L1_5_S5_L001_R1_001_AGAC.trim.bt2.bam \
IAN_2brad_L1_5_S5_L001_R1_001_AGTG.trim.bt2.bam \
IAN_2brad_L1_5_S5_L001_R1_001_CATC.trim.bt2.bam \
IAN_2brad_L1_5_S5_L001_R1_001_GACT.trim.bt2.bam \
IAN_2brad_L1_5_S5_L001_R1_001_GCTT.trim.bt2.bam \
IAN_2brad_L1_5_S5_L001_R1_001_GTGA.trim.bt2.bam \
IAN_2brad_L1_5_S5_L001_R1_001_TCAC.trim.bt2.bam \
IAN_2brad_L1_5_S5_L001_R1_001_TCAG.trim.bt2.bam \
IAN_2brad_L1_5_S5_L001_R1_001_TGTC.trim.bt2.bam \
IAN_2brad_L1_8_S8_L001_R1_001_AGAC.trim.bt2.bam \
IAN_2brad_L1_8_S8_L001_R1_001_AGTG.trim.bt2.bam \
IAN_2brad_L1_8_S8_L001_R1_001_CATC.trim.bt2.bam \
IAN_2brad_L1_8_S8_L001_R1_001_GCTT.trim.bt2.bam \
IAN_2brad_L1_8_S8_L001_R1_001_TCAC.trim.bt2.bam \
IAN_2brad_L1_8_S8_L001_R1_001_TCAG.trim.bt2.bam \
IAN_2brad_L2_1_S9_L002_R1_001_GACT.trim.bt2.bam \
IAN_2brad_L2_1_S9_L002_R1_001_TCAC.trim.bt2.bam \
IAN_2brad_L2_1_S9_L002_R1_001_TGTC.trim.bt2.bam \
IAN_2brad_L2_2_S10_L002_R1_001_GACT.trim.bt2.bam \
IAN_2brad_L2_2_S10_L002_R1_001_TCAC.trim.bt2.bam \
IAN_2brad_L2_3_S11_L002_R1_001_AGAC.trim.bt2.bam \
IAN_2brad_L2_3_S11_L002_R1_001_AGTG.trim.bt2.bam \
IAN_2brad_L2_3_S11_L002_R1_001_CATC.trim.bt2.bam \
IAN_2brad_L2_3_S11_L002_R1_001_GTGA.trim.bt2.bam \
IAN_2brad_L2_3_S11_L002_R1_001_TCAG.trim.bt2.bam \
IAN_2brad_L2_3_S11_L002_R1_001_TGGT.trim.bt2.bam \
IAN_2brad_L2_4_S12_L002_R1_001_ACCA.trim.bt2.bam \
IAN_2brad_L2_4_S12_L002_R1_001_AGAC.trim.bt2.bam \
IAN_2brad_L2_4_S12_L002_R1_001_AGTG.trim.bt2.bam \
IAN_2brad_L2_4_S12_L002_R1_001_CATC.trim.bt2.bam \
IAN_2brad_L2_4_S12_L002_R1_001_GTGA.trim.bt2.bam \
IAN_2brad_L2_4_S12_L002_R1_001_TGGT.trim.bt2.bam \
IAN_2brad_L3_5_S5_L001_R1_001_GACT.trim.bt2.bam \
IAN_2brad_L3_6_S6_L001_R1_001_ACCA.trim.bt2.bam \
IAN_2brad_L3_6_S6_L001_R1_001_AGAC.trim.bt2.bam \
IAN_2brad_L3_6_S6_L001_R1_001_CTAC.trim.bt2.bam \
IAN_2brad_L3_6_S6_L001_R1_001_GACT.trim.bt2.bam \
IAN_2brad_L3_6_S6_L001_R1_001_GCTT.trim.bt2.bam \
IAN_2brad_L3_6_S6_L001_R1_001_GTGA.trim.bt2.bam \
IAN_2brad_L3_6_S6_L001_R1_001_TCAC.trim.bt2.bam \
IAN_2brad_L3_6_S6_L001_R1_001_TCAG.trim.bt2.bam \
IAN_2brad_L3_6_S6_L001_R1_001_TGGT.trim.bt2.bam \
IAN_2brad_L3_6_S6_L001_R1_001_TGTC.trim.bt2.bam \
IAN_2brad_L3_7_S7_L001_R1_001_ACCA.trim.bt2.bam \
IAN_2brad_L3_7_S7_L001_R1_001_AGAC.trim.bt2.bam \
IAN_2brad_L3_7_S7_L001_R1_001_AGTG.trim.bt2.bam \
IAN_2brad_L3_7_S7_L001_R1_001_CATC.trim.bt2.bam \
IAN_2brad_L3_7_S7_L001_R1_001_GCTT.trim.bt2.bam \
IAN_2brad_L3_7_S7_L001_R1_001_GTGA.trim.bt2.bam \
IAN_2brad_L3_7_S7_L001_R1_001_TCAC.trim.bt2.bam \
IAN_2brad_L3_7_S7_L001_R1_001_TCAG.trim.bt2.bam \
IAN_2brad_L3_7_S7_L001_R1_001_TGGT.trim.bt2.bam \
IAN_2brad_L3_7_S7_L001_R1_001_TGTC.trim.bt2.bam \
IAN_2brad_L3_8_S8_L001_R1_001_ACCA.trim.bt2.bam \
IAN_2brad_L3_8_S8_L001_R1_001_AGAC.trim.bt2.bam \
IAN_2brad_L3_8_S8_L001_R1_001_AGTG.trim.bt2.bam \
IAN_2brad_L3_8_S8_L001_R1_001_TGGT.trim.bt2.bam \
IAN_2brad_L4_2_S10_L002_R1_001_ACCA.trim.bt2.bam \
IAN_2brad_L4_2_S10_L002_R1_001_AGAC.trim.bt2.bam \
IAN_2brad_L4_2_S10_L002_R1_001_AGTG.trim.bt2.bam \
IAN_2brad_L4_2_S10_L002_R1_001_TCAG.trim.bt2.bam \
IAN_2brad_L4_2_S10_L002_R1_001_TGGT.trim.bt2.bam \
IAN_2brad_L4_3_S11_L002_R1_001_CTAC.trim.bt2.bam > mass_bams

#Count them
cat mass_bams | wc -l
# 84 good

### For new analysis #####
# Make new file with all these massbams, that will also include mass tech reps, called mass_bams_with_techreps
# Tech reps to add:
IAN_2brad_L1_8_S8_L001_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_GTGA.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_CTAC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_TGTC.trim.bt2.bam
IAN_2brad_L1_8_S8_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L2_5_S13_L002_R1_001_ACCA.trim.bt2.bam
IAN_2brad_L3_6_S6_L001_R1_001_CATC.trim.bt2.bam
IAN_2brad_L3_7_S7_L001_R1_001_GACT.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TCAG.trim.bt2.bam
IAN_2brad_L4_3_S11_L002_R1_001_TGGT.trim.bt2.bam
IAN_2brad_L1_1_S1_L001_R1_001_TCAC.trim.bt2.bam
IAN_2brad_L1_5_S5_L001_R1_001_CTAC.trim.bt2.bam


# Made new qsub here called mass_bams.qsub
# Using minind 80%, minindepth 3

module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 67 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b mass_bams -GL 1 $FILTERS $TODO -P 1 -out v2_myresult_massbams

NSITES=`zcat v2_myresult_massbams.mafs.gz | wc -l`
echo $NSITES
# 6611

# Using minind 80%, minindepth 3, everything for fis so changing dogeno to 2
#qsub is mass_bams_fis.qsub

module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 67 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 2 -doBcf 1 -postCutoff 0.95 -doPost 1 -doGlf 2"

angsd -b mass_bams -GL 1 $FILTERS $TODO -P 1 -out mass_bams_fis

NSITES=`zcat mass_bams_fis.mafs.gz | wc -l`
echo $NSITES
# 

## For mass bams with tech reps:
# Made new qsub here called mass_bams_techreps.qsub
# Using minind 80%, minindepth 3

module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 3"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b mass_bams_with_techreps -GL 1 $FILTERS $TODO -P 1 -out mass_bams_techreps

NSITES=`zcat mass_bams_techreps.mafs.gz | wc -l`
echo $NSITES
# 6209

# Final for mass bams
#v2_myresult_massbams
#mass_bams_techreps

# cyberduck to computer to look at who is who


### Checking final_50p_nodups_fis.geno (for FIS) for missing data
# Check across dataset
grep -o "\-1" final_50p_nodups_fis.geno | wc -l

# 74038

# Total proportion of missing data
total_sites=$(awk '{count += NF} END {print count}' final_50p_nodups_fis.geno)
missing_sites=$(grep -o "\-1" final_50p_nodups_fis.geno | wc -l)
echo "scale=4; $missing_sites / $total_sites" | bc

#0.0569

# Missing data across SNPs
awk '{for (i=1; i<=NF; i++) if ($i == -1) count[i]++ } END {for (i in count) print i, count[i]}' final_50p_nodups_fis.geno | sort -k2,2nr | head
#222 1634
#247 1264
#144 1214
#209 1153
#39 1135
#55 1071
#16 1034
#6 961
#29 950
#17 904

# ---------------------------------------- Maximum Likelihood Analysis using RaxML

# Make a new directory called raxml, located here:
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/raxml

# Copy the myresult_alldups_12p_nodups.bcf file as well as the raxml.qsub file from the UCE analysis over

cp ANGSD/myresult_alldups_12p_nodups.bcf raxml/
cp /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/raxml.qsub raxml

# Going to convert the bcf file to a vcf file

module load htslib/1.16
module load bcftools/1.16
bcftools view myresult_alldups_12p_nodups.bcf -O v -o myresult_alldups_12p_nodups.vcf

# Now convert the new vcf to phylip so it can be fed to raxml (takes phylip and nexus)
# Clone the python script from https://github.com/edgardomortiz/vcf2phylip.git

git clone https://github.com/edgardomortiz/vcf2phylip.git

# This made its own directory within raxml, so we're going to copy the .py script to the raxml directory

cd vcf2phylip
cp vcf2phylip.py ../

# Where are we?

[inovick@scc1 raxml]$ pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/raxml

# Use default parameters to create a PHYLIP matrix with a minimum of 4 samples per SNP
# Execute the command to get the vcf to phylip format

python vcf2phylip.py --input myresult_alldups_12p_nodups.vcf --output-prefix myresult_alldups_12p_nodups.phy

# For this conversion from vcf to phylip file, there are some differing codes for polymorphic sites. These are:
# '*' is a deletion in GATK

# AG -> R
# CT -> Y
# ACGT -> N
# *ACGT -> N
# AT -> W
# GT -> K

# There might be more, but these are the most common ones I've been seeing in the phylip file. To check other ones, look in the vcf2phylip.py file


# Now we have it and it's called myresult_alldups_12p_nodups.phy.min4.phy
# This tree doesn't have an outgroup so it's going to be unrooted
# Going to run raxml qsub:

module load raxml/8.2.9
raxmlHPC -f a -p 12345 -x 12345 -N 100 -m GTRGAMMAI -s myresult_alldups_12p_nodups.phy.min4.phy -n raxml_myresult_alldups_12p_nodups

# All done
# Now cyberduck RAxML_bipartitions.raxml_myresult_alldups_12p_nodups to local machine to look at in figtree
# Going to have to change the names to city/state to actually understand what we're looking at
# Using FindAndReplace_uce_global_samples.sh script, but going to put in all bam names with city, state, and country

# ---------------------------------------- Admixture analyses: NGSadmix and PCANGSD

# Starting with NGSadmix
# New directory for admixture using NGSadmix called 
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/admixture/NGSadmix

# Copy over relevant files from ANGSD directory
# Going to do the Mass samples and the global samples

cp v2_myresult_massbams.beagle.gz ../admixture/NGSadmix/
cp myresult_alldups_12p_nodups.beagle.gz ../admixture/NGSadmix/


# Make a bin directory
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/bin

# Install and download to bin
wget popgen.dk/software/download/NGSadmix/ngsadmix32.cpp 

g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix


# All set, downloaded in bin
# I also copied it and NGSadmix to my working directory

# For global samples:
# NgsAdmix for K from 2 to 11 first : do not run if the dataset contains clones or genotyping replicates!

for K in `seq 2 11` ;
do
NGSadmix  -likes myresult_alldups_12p_nodups.beagle.gz -K $K -P 10 -o mydata.allcountries_k${K};
done

# cyberduck files to local machine, use admixturePlotting_v5.R to plot in Rstudio
pwd 
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/admixture/NGSadmix

# On local machine, path is:
/Desktop/Tineoid_Moth_Project/Chapter_2_RAD/Admixture/global_samples

# For Mass samples
# NgsAdmix for K from 2 to 5 first : do not run if the dataset contains clones or genotyping replicates!

for K in `seq 2 5` ;
do
NGSadmix  -likes v2_myresult_massbams.beagle.gz -K $K -P 10 -o mydata.mass_k${K};
done

# cyberduck files to local machine, use admixturePlotting_v5.R to plot in Rstudio
pwd 
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/admixture/NGSadmix

# On local machine, path is:
/Desktop/Tineoid_Moth_Project/Chapter_2_RAD/Admixture/mass_samples

# Now going to select the optimal k using https://baylab.github.io/MarineGenomics/week-9-population-structure-using-ngsadmix.html section 10.4
# We should run NGSadmix 2 more times for both mass and global samples to get an accurate optimal k

# For global samples replicate 2:
for K in `seq 2 11` ;
do
NGSadmix  -likes myresult_alldups_12p_nodups.beagle.gz -K $K -P 10 -o mydata.allcountries_k${K}_rep2;
done

# For global samples replicate 2:
for K in `seq 2 11` ;
do
NGSadmix  -likes myresult_alldups_12p_nodups.beagle.gz -K $K -P 10 -o mydata.allcountries_k${K}_rep3;
done

# For mass samples replicate 2:
for K in `seq 2 5` ;
do
NGSadmix  -likes v2_myresult_massbams.beagle.gz -K $K -P 10 -o mydata.mass_k${K}_rep2;
done

# For mass samples replicate 3:
for K in `seq 2 5` ;
do
NGSadmix  -likes v2_myresult_massbams.beagle.gz -K $K -P 10 -o mydata.mass_k${K}_rep3;
done

# Optimal k global samples rep 1: optimal k is 5 with 52.135971 probability
# Optimal k global samples rep 2: optimal k is 5 with 52.135971 probability
# Optimal k global samples rep 3: optimal k is 5 with 52.135971 probability

# Optimal k mass samples rep 1: optimal k is 5 with 23.820069 probability
# Optimal k mass samples rep 2: optimal k is 5 with 23.82007 probability
# Optimal k mass samples rep 3: optimal k is 5 with 23.82007 probability


#------------------------------Analysis of Genetic Divergence (Pairwise Fst)
#### Global samples first
# For global samples, made a new directory:
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/Fst/global_samples

# Need to copy over global samples bams list and all global samples (no dups, 12p coverage that I've been using for all analyses)
# First copy over the bams list
cd ../../ANGSD
cp alldups_12p_bams_nodups ../Fst/global_samples/

# Cyberduck alldups_12p_bams_nodups to local machine
# Save as new file on local machine called:
/Desktop/Tineoid_Moth_Project/Chapter_2_RAD/Fst/alldups_12p_bams_nodups_cp_global_samps.txt

# Now copy over all relevant bams files by copying the content of that into the terminal to copy all relevant bams to global_samples dir

# Count how many bam files there are in the dir
find . -type f -name "*.trim.bt2.bam" | wc -l

# 284
# Great

# Now need to make 5 bams lists corresponding to the 5 populations identified by NGSadmix
# How to identify the 5 populations? Use Rscript found in ngsadmix.R line 260:

# recording cluster affiliations
cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.25),collapse=".")) })
cluster.admix.df = as.data.frame(cluster.admix)

# make the rownames a column so we can combine with metadata
cluster.admix.df = cluster.admix.df %>%
  rownames_to_column("ind")
  
# now combine with metadata file
cluster.admix.i2p = left_join(i2p, cluster.admix.df)
write.csv(cluster.admix.i2p, "/Users/izzynovick/Desktop/Tineoid_Moth_Project/Chapter_2_RAD/Admixture/global_samples/admixture_cluster_assignments.csv", row.names=FALSE)

# Now have file called 
/Users/izzynovick/Desktop/Tineoid_Moth_Project/Chapter_2_RAD/Admixture/global_samples/admixture_cluster_assignments.csv


# And if I'm looking at this file:
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/admixture/NGSadmix/mydata.allcountries_k5.qopt

# It has 5 columns for each sample with the proportion of admixture for each lineage. I'm going to make a new file called 
/Users/izzynovick/Desktop/Tineoid_Moth_Project/Chapter_2_RAD/Fst/global_samples/admixture_cluster_assignments_modified_basedon_k5_qopt

# Which is basically this file:
admixture_cluster_assignments.csv

# But I'm changing the assignments for those individuals who belong to 2 populations to be the one that they have a higher association with, found in the qopt file

# Making bams lists for each of the 5 lineages, found in:
/Users/izzynovick/Desktop/Tineoid_Moth_Project/Chapter_2_RAD/Fst/global_samples/

# Called:
lineage1_bams_list.txt
lineage2_bams_list.txt
lineage3_bams_list.txt
lineage4_bams_list.txt
lineage5_bams_list.txt

# These 5 lineages are not delineated based on geography. Some are, but most are a mix of different places. Esp. lineage 2, which has a lot of samples, contains samples from all over the world

# Cyberduck to scc directory:
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/Fst/global_samples

# With all clones removed and everything below 12% coverage removed, we have 284 samples

# filtering sites to work on - use only filters that do not distort allele frequency
# set minInd to 75-90% of the total number fo individuals in the project-- I'm going to use 90% of 284 which is 256 for all Fst runs
# if you are doing any other RAD than 2bRAD or GBS, remove '-sb_pval 1e-5' from FILTERS
# -maxHetFreq 0.5 is the lumped paralog filter - important for doing SFS analyses
# Make qsub for first angsd run called global_bams_filt.qsub 

module load angsd/0.923


FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 256"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doVcf 1"


angsd -b alldups_12p_bams_nodups -GL 1 $FILTERS $TODO -out global_bams_90p


# Count sites
NSITES=`zcat global_bams_90p.mafs.gz | wc -l`
echo $NSITES
# 725808

# -maxHetFreq 0.5 is the lumped paralog filter - important for doing SFS analyses. Going to run this also to see if it makes a difference for Fst values
# Make new qsub called global_bams_filt_maxhet.qsub

module load angsd


FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 256 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doBcf 1"


angsd -b alldups_12p_bams_nodups -GL 1 $FILTERS $TODO -out global_bams_90p_maxhet

# Count sites
NSITES=`zcat global_bams_90p_maxhet.mafs.gz | wc -l`
echo $NSITES
# 725803

# Make vcf only for running bayescan
# Qsub called bayescan.qsub

module load angsd


FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 256 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doBcf 1"

angsd -b alldups_12p_bams_nodups -GL 1 $FILTERS $TODO -out basescan_vcf_global_bams

# Collecting and indexing filter-passing sites
module load angsd 

zcat global_bams_90p.mafs.gz | cut -f 1,2 | tail -n +2 >allSites
angsd sites index allSites

# max het filter
zcat global_bams_90p_maxhet.mafs.gz | cut -f 1,2 | tail -n +2 >allSites_maxhet
angsd sites index allSites_maxhet

# Reindex genome in current working directory. First copy the genome:

cd /projectnb/mullenl/novick/global_pop_gen/rad_workflow/
cp wcm_pseudochromosomes_renamed.fasta Fst/global_samples/

# Then reindex here
cd /projectnb/mullenl/novick/global_pop_gen/rad_workflow/Fst/global_samples

module load samtools
samtools faidx wcm_pseudochromosomes_renamed.fasta

export GENOME_REF=wcm_pseudochromosomes_renamed.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"

# Assigning the populations 

# Without maxhet filter
angsd -sites allSites -b lineage1_bams_list.txt -GL 1 -P 1 $TODO -out L1
angsd -sites allSites -b lineage2_bams_list.txt -GL 1 -P 1 $TODO -out L2
angsd -sites allSites -b lineage3_bams_list.txt -GL 1 -P 1 $TODO -out L3
angsd -sites allSites -b lineage4_bams_list.txt -GL 1 -P 1 $TODO -out L4
angsd -sites allSites -b lineage5_bams_list.txt -GL 1 -P 1 $TODO -out L5

# max het filter
angsd -sites allSites_maxhet -b lineage1_bams_list.txt -GL 1 -P 1 $TODO -out L1_maxhet
angsd -sites allSites_maxhet -b lineage2_bams_list.txt -GL 1 -P 1 $TODO -out L2_maxhet
angsd -sites allSites_maxhet -b lineage3_bams_list.txt -GL 1 -P 1 $TODO -out L3_maxhet
angsd -sites allSites_maxhet -b lineage4_bams_list.txt -GL 1 -P 1 $TODO -out L4_maxhet
angsd -sites allSites_maxhet -b lineage5_bams_list.txt -GL 1 -P 1 $TODO -out L5_maxhet


# generating per-population SFS
realSFS L1.saf.idx >L1.sfs
realSFS L2.saf.idx >L2.sfs
realSFS L3.saf.idx >L3.sfs
realSFS L4.saf.idx >L4.sfs
realSFS L5.saf.idx >L5.sfs

# generating per-population SFS - max het filter
realSFS L1_maxhet.saf.idx >L1_maxhet.sfs
realSFS L2_maxhet.saf.idx >L2_maxhet.sfs
realSFS L3_maxhet.saf.idx >L3_maxhet.sfs
realSFS L4_maxhet.saf.idx >L4_maxhet.sfs
realSFS L5_maxhet.saf.idx >L5_maxhet.sfs

### Pairwise comparisons:

# writing down 2d-SFS priors - L1 vs L2 max het
realSFS L1_maxhet.saf.idx L2_maxhet.saf.idx -P 24 > p1_2_maxhet.sfs ; realSFS fst index L1_maxhet.saf.idx L2_maxhet.saf.idx -sfs p1_2_maxhet.sfs -fstout p1_2_maxhet

# writing down 2d-SFS priors - L1 vs L2 without max het
realSFS L1.saf.idx L2.saf.idx -P 24 > p1_2.sfs ; realSFS fst index L1.saf.idx L2.saf.idx -sfs p1_2.sfs -fstout p1_2

# global Fst between populations-- max het
realSFS fst stats p1_2_maxhet.fst.idx

# Output L1 and L2 max het:
	-> Assuming idxname:p1_2_maxhet.fst.idx
	-> Assuming .fst.gz file: p1_2_maxhet.fst.gz
	-> FST.Unweight[nObs:725774]:-0.000735 Fst.Weight:0.125999
-0.000735	0.125999

# global Fst between populations-- no max het
realSFS fst stats p1_2.fst.idx

# Output L1 and L2 no max het:
	-> Assuming idxname:p1_2.fst.idx
	-> Assuming .fst.gz file: p1_2.fst.gz
	-> FST.Unweight[nObs:725779]:-0.000734 Fst.Weight:0.125989
-0.000734	0.125989

# Making a qsub containing all max het and non max het comparisons called site_comparisons.qsub

# writing down 2d-SFS priors - L1 vs L3 with max het
realSFS L1_maxhet.saf.idx L3_maxhet.saf.idx -P 24 > p1_3_maxhet.sfs ; realSFS fst index L1_maxhet.saf.idx L3_maxhet.saf.idx -sfs p1_3_maxhet.sfs -fstout p1_3_maxhet

# global Fst between populations-- max het
realSFS fst stats p1_3_maxhet.fst.idx

# Output L1 and L3 max het:
	-> Assuming idxname:p1_3_maxhet.fst.idx
	-> Assuming .fst.gz file: p1_3_maxhet.fst.gz
	-> FST.Unweight[nObs:725794]:0.012211 Fst.Weight:0.337713
0.012211	0.337713

# writing down 2d-SFS priors - L1 vs L3 without max het
realSFS L1.saf.idx L3.saf.idx -P 24 > p1_3.sfs ; realSFS fst index L1.saf.idx L3.saf.idx -sfs p1_3.sfs -fstout p1_3

# global Fst between populations-- no max het
realSFS fst stats p1_3.fst.idx

# Output L1 and L3 no max het:
	-> Assuming idxname:p1_3.fst.idx
	-> Assuming .fst.gz file: p1_3.fst.gz
	-> FST.Unweight[nObs:725799]:0.012213 Fst.Weight:0.337659
0.012213	0.337659

# writing down 2d-SFS priors - L1 vs L4 with max het
realSFS L1_maxhet.saf.idx L4_maxhet.saf.idx -P 24 > p1_4_maxhet.sfs ; realSFS fst index L1_maxhet.saf.idx L4_maxhet.saf.idx -sfs p1_4_maxhet.sfs -fstout p1_4_maxhet

# global Fst between populations-- max het
realSFS fst stats p1_4_maxhet.fst.idx

# Output L1 and L4 max het:
	-> Assuming idxname:p1_4_maxhet.fst.idx
	-> Assuming .fst.gz file: p1_4_maxhet.fst.gz
	-> FST.Unweight[nObs:725795]:0.011435 Fst.Weight:0.269860
0.011435	0.269860

# writing down 2d-SFS priors - L1 vs L4 without max het
realSFS L1.saf.idx L4.saf.idx -P 24 > p1_4.sfs ; realSFS fst index L1.saf.idx L4.saf.idx -sfs p1_4.sfs -fstout p1_4

# global Fst between populations-- no max het
realSFS fst stats p1_4.fst.idx

# Output L1 and L4 no max het:
	-> Assuming idxname:p1_4.fst.idx
	-> Assuming .fst.gz file: p1_4.fst.gz
	-> FST.Unweight[nObs:725800]:0.011436 Fst.Weight:0.269701
0.011436	0.269701

# writing down 2d-SFS priors - L1 vs L5 with max het
realSFS L1_maxhet.saf.idx L5_maxhet.saf.idx -P 24 > p1_5_maxhet.sfs ; realSFS fst index L1_maxhet.saf.idx L5_maxhet.saf.idx -sfs p1_5_maxhet.sfs -fstout p1_5_maxhet

# global Fst between populations-- max het
realSFS fst stats p1_5_maxhet.fst.idx

# Output L1 and L5 max het:
	-> Assuming idxname:p1_5_maxhet.fst.idx
	-> Assuming .fst.gz file: p1_5_maxhet.fst.gz
	-> FST.Unweight[nObs:725797]:0.003452 Fst.Weight:0.437056
0.003452	0.437056

# writing down 2d-SFS priors - L1 vs L5 without max het
realSFS L1.saf.idx L5.saf.idx -P 24 > p1_5.sfs ; realSFS fst index L1.saf.idx L5.saf.idx -sfs p1_5.sfs -fstout p1_5

# global Fst between populations-- no max het
realSFS fst stats p1_5.fst.idx

# Output L1 and L5 no max het:
	-> Assuming idxname:p1_5.fst.idx
	-> Assuming .fst.gz file: p1_5.fst.gz
	-> FST.Unweight[nObs:725802]:0.003452 Fst.Weight:0.436996
0.003452	0.436996

# writing down 2d-SFS priors - L2 vs L3 with max het
realSFS L2_maxhet.saf.idx L3_maxhet.saf.idx -P 24 > p2_3_maxhet.sfs ; realSFS fst index L2_maxhet.saf.idx L3_maxhet.saf.idx -sfs p2_3_maxhet.sfs -fstout p2_3_maxhet

# global Fst between populations-- max het
realSFS fst stats p2_3_maxhet.fst.idx

# Output L2 and L3 max het:
	-> Assuming idxname:p2_3_maxhet.fst.idx
	-> Assuming .fst.gz file: p2_3_maxhet.fst.gz
	-> FST.Unweight[nObs:725774]:-0.002423 Fst.Weight:0.171561
-0.002423	0.171561

# writing down 2d-SFS priors - L2 vs L3 without max het
realSFS L2.saf.idx L3.saf.idx -P 24 > p2_3.sfs ; realSFS fst index L2.saf.idx L3.saf.idx -sfs p2_3.sfs -fstout p2_3

# global Fst between populations-- no max het
realSFS fst stats p2_3.fst.idx

# Output L2 and L3 no max het:
[inovick@scc1 global_samples]$ realSFS fst stats p2_3.fst.idx
	-> Assuming idxname:p2_3.fst.idx
	-> Assuming .fst.gz file: p2_3.fst.gz
	-> FST.Unweight[nObs:725779]:-0.002422 Fst.Weight:0.171439
-0.002422	0.171439

# writing down 2d-SFS priors - L2 vs L4 with max het
realSFS L2_maxhet.saf.idx L4_maxhet.saf.idx -P 24 > p2_4_maxhet.sfs ; realSFS fst index L2_maxhet.saf.idx L4_maxhet.saf.idx -sfs p2_4_maxhet.sfs -fstout p2_4_maxhet

# global Fst between populations-- max het
realSFS fst stats p2_4_maxhet.fst.idx

# Output L2 and L4 max het:
	-> Assuming idxname:p2_4_maxhet.fst.idx
	-> Assuming .fst.gz file: p2_4_maxhet.fst.gz
	-> FST.Unweight[nObs:725774]:-0.001597 Fst.Weight:0.126299
-0.001597	0.126299

# writing down 2d-SFS priors - L2 vs L4 without max het
realSFS L2.saf.idx L4.saf.idx -P 24 > p2_4.sfs ; realSFS fst index L2.saf.idx L4.saf.idx -sfs p2_4.sfs -fstout p2_4

# global Fst between populations-- no max het
realSFS fst stats p2_4.fst.idx

# Output L2 and L4 no max het:
	-> Assuming idxname:p2_4.fst.idx
	-> Assuming .fst.gz file: p2_4.fst.gz
	-> FST.Unweight[nObs:725779]:-0.001597 Fst.Weight:0.126181
-0.001597	0.126181

# writing down 2d-SFS priors - L2 vs L5 with max het
realSFS L2_maxhet.saf.idx L5_maxhet.saf.idx -P 24 > p2_5_maxhet.sfs ; realSFS fst index L2_maxhet.saf.idx L5_maxhet.saf.idx -sfs p2_5_maxhet.sfs -fstout p2_5_maxhet

# global Fst between populations-- max het
realSFS fst stats p2_5_maxhet.fst.idx

# Output L2 and L5 max het:
	-> Assuming idxname:p2_5_maxhet.fst.idx
	-> Assuming .fst.gz file: p2_5_maxhet.fst.gz
	-> FST.Unweight[nObs:725774]:-0.005721 Fst.Weight:0.260042
-0.005721	0.260042

# writing down 2d-SFS priors - L2 vs L5 without max het
realSFS L2.saf.idx L5.saf.idx -P 24 > p2_5.sfs ; realSFS fst index L2.saf.idx L5.saf.idx -sfs p2_5.sfs -fstout p2_5

# global Fst between populations-- no max het
realSFS fst stats p2_5.fst.idx

# Output L2 and L5 no max het:
	-> Assuming idxname:p2_5.fst.idx
	-> Assuming .fst.gz file: p2_5.fst.gz
	-> FST.Unweight[nObs:725779]:-0.005719 Fst.Weight:0.260033
-0.005719	0.260033

# writing down 2d-SFS priors - L3 vs L4 with max het
realSFS L3_maxhet.saf.idx L4_maxhet.saf.idx -P 24 > p3_4_maxhet.sfs ; realSFS fst index L3_maxhet.saf.idx L4_maxhet.saf.idx -sfs p3_4_maxhet.sfs -fstout p3_4_maxhet

# global Fst between populations-- max het
realSFS fst stats p3_4_maxhet.fst.idx

# Output L3 and L4 max het:
	-> Assuming idxname:p3_4_maxhet.fst.idx
	-> Assuming .fst.gz file: p3_4_maxhet.fst.gz
	-> FST.Unweight[nObs:725798]:0.012284 Fst.Weight:0.325756
0.012284	0.325756

# writing down 2d-SFS priors - L3 vs L4 without max het
realSFS L3.saf.idx L4.saf.idx -P 24 > p3_4.sfs ; realSFS fst index L3.saf.idx L4.saf.idx -sfs p3_4.sfs -fstout p3_4

# global Fst between populations-- no max het
realSFS fst stats p3_4.fst.idx

# Output L3 and L4 no max het:
	-> Assuming idxname:p3_4.fst.idx
	-> Assuming .fst.gz file: p3_4.fst.gz
	-> FST.Unweight[nObs:725803]:0.012285 Fst.Weight:0.325512
0.012285	0.325512

# writing down 2d-SFS priors - L3 vs L5 with max het
realSFS L3_maxhet.saf.idx L5_maxhet.saf.idx -P 24 > p3_5_maxhet.sfs ; realSFS fst index L3_maxhet.saf.idx L5_maxhet.saf.idx -sfs p3_5_maxhet.sfs -fstout p3_5_maxhet

# global Fst between populations-- max het
realSFS fst stats p3_5_maxhet.fst.idx

# Output L3 and L5 max het:
	-> Assuming idxname:p3_5_maxhet.fst.idx
	-> Assuming .fst.gz file: p3_5_maxhet.fst.gz
	-> FST.Unweight[nObs:725799]:0.009380 Fst.Weight:0.511383
0.009380	0.511383

# writing down 2d-SFS priors - L3 vs L5 without max het
realSFS L3.saf.idx L5.saf.idx -P 24 > p3_5.sfs ; realSFS fst index L3.saf.idx L5.saf.idx -sfs p3_5.sfs -fstout p3_5

# global Fst between populations-- no max het
realSFS fst stats p3_5.fst.idx

# Output L3 and L5 no max het:
	-> Assuming idxname:p3_5.fst.idx
	-> Assuming .fst.gz file: p3_5.fst.gz
	-> FST.Unweight[nObs:725804]:0.009379 Fst.Weight:0.511192
0.009379	0.511192

# writing down 2d-SFS priors - L4 vs L5 with max het
realSFS L4_maxhet.saf.idx L5_maxhet.saf.idx -P 24 > p4_5_maxhet.sfs ; realSFS fst index L4_maxhet.saf.idx L5_maxhet.saf.idx -sfs p4_5_maxhet.sfs -fstout p4_5_maxhet

# global Fst between populations-- max het
realSFS fst stats p4_5_maxhet.fst.idx

# Output L4 and L5 max het:
	-> Assuming idxname:p4_5_maxhet.fst.idx
	-> Assuming .fst.gz file: p4_5_maxhet.fst.gz
	-> FST.Unweight[nObs:725800]:0.006762 Fst.Weight:0.424046
0.006762	0.424046

# writing down 2d-SFS priors - L4 vs L5 without max het
realSFS L4.saf.idx L5.saf.idx -P 24 > p4_5.sfs ; realSFS fst index L4.saf.idx L5.saf.idx -sfs p4_5.sfs -fstout p4_5

# global Fst between populations-- no max het
realSFS fst stats p4_5.fst.idx

# Output L4 and L5 no max het:
	-> Assuming idxname:p4_5.fst.idx
	-> Assuming .fst.gz file: p4_5.fst.gz
	-> FST.Unweight[nObs:725805]:0.006764 Fst.Weight:0.423974
0.006764	0.423974

#### Massachusetts samples 
# For mass samples, made a new directory:
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/Fst/mass_samples

# Need to copy over global samples bams list and all global samples (no dups, 12p coverage that I've been using for all analyses)
# First copy over the bams list
cd ../../ANGSD
cp mass_bams ../Fst/mass_samples/

# Cyberduck alldups_12p_bams_nodups to local machine
# Made new file on local machine called:
/Desktop/Tineoid_Moth_Project/Chapter_2_RAD/Fst/mass_bams_cp_mass_samps.txt

# Now copy over all relevant bams files by copying the content of that into the terminal to copy all relevant bams to mass_samples dir

# Count how many bam files there are in the dir
find . -type f -name "*.trim.bt2.bam" | wc -l

# 84
# Awesome

# Now need to make 5 bams lists corresponding to the 5 towns in Mass since NGSadmix indicated that each town is its own distinct population
# Manually making each one based on the massachusetts_bams_list excel file

# Boston:
boston_pop_bams_list.txt

# Cambridge:
cambridge_pop_bams_list.txt

# Arlington:
arlington_pop_bams_list.txt

# Belmont:
belmont_pop_bams_list.txt

# Waltham:
waltham_pop_bams_list.txt

# Cyberduck all to 
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/Fst/mass_samples

# with all clones removed, we have 84 samples

module load angsd/0.923

# filtering sites to work on - use only filters that do not distort allele frequency
# set minInd to 75-90% of the total number fo individuals in the project-- going to do 90% (76 samples)
# if you are doing any other RAD than 2bRAD or GBS, remove '-sb_pval 1e-5' from FILTERS
# -maxHetFreq 0.5 is the lumped paralog filter - important for doing SFS analyses
# Made qsub called mass_bams_filt.qsub 

module load angsd/0.923

FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 76"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doVcf 1"


angsd -b mass_bams -GL 1 $FILTERS $TODO -out mass_bams_90p

# Count sites
NSITES=`zcat mass_bams_90p.mafs.gz | wc -l`
echo $NSITES

# 729167

# -maxHetFreq 0.5 is the lumped paralog filter - important for doing SFS analyses. seeing if it makes a difference for Fst values
# Made qsub called mass_bams_filt_maxhet.qsub

module load angsd # I removed the version because it wasn't recognizing the maxhetfreq filter in v.923 for some reason

FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 76 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doBcf 1"


angsd -b mass_bams -GL 1 $FILTERS $TODO -out mass_bams_90p_maxhet

# Count sites
NSITES=`zcat mass_bams_90p_maxhet.mafs.gz | wc -l`
echo $NSITES

# 729079

# make vcf only for running bayescan:
# Qsub called bayescan.qsub

FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 76 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doBcf 1"

angsd -b mass_bams -GL 1 $FILTERS $TODO -out basescan_vcf_mass_bams

# collecting and indexing filter-passing sites
zcat mass_bams_90p.mafs.gz | cut -f 1,2 | tail -n +2 >allSites
angsd sites index allSites

# max het filter
zcat mass_bams_90p_maxhet.mafs.gz | cut -f 1,2 | tail -n +2 >allSites_maxhet
angsd sites index allSites_maxhet

# Reindex genome in current working directory. First copy the genome:

cd /projectnb/mullenl/novick/global_pop_gen/rad_workflow/
cp wcm_pseudochromosomes_renamed.fasta Fst/mass_samples/

# Then reindex here
module load samtools
samtools faidx wcm_pseudochromosomes_renamed.fasta

export GENOME_REF=wcm_pseudochromosomes_renamed.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"



# Without maxhet
angsd -sites allSites -b boston_pop_bams_list.txt -GL 1 -P 1 $TODO -out boston_L1
angsd -sites allSites -b cambridge_pop_bams_list.txt -GL 1 -P 1 $TODO -out cambridge_L2
angsd -sites allSites -b arlington_pop_bams_list.txt -GL 1 -P 1 $TODO -out arlington_L3
angsd -sites allSites -b belmont_pop_bams_list.txt -GL 1 -P 1 $TODO -out belmont_L4
angsd -sites allSites -b waltham_pop_bams_list.txt -GL 1 -P 1 $TODO -out waltham_L5

# With maxhet filter
angsd -sites allSites_maxhet -b boston_pop_bams_list.txt -GL 1 -P 1 $TODO -out boston_L1_maxhet
angsd -sites allSites_maxhet -b cambridge_pop_bams_list.txt -GL 1 -P 1 $TODO -out cambridge_L2_maxhet
angsd -sites allSites_maxhet -b arlington_pop_bams_list.txt -GL 1 -P 1 $TODO -out arlington_L3_maxhet
angsd -sites allSites_maxhet -b belmont_pop_bams_list.txt -GL 1 -P 1 $TODO -out belmont_L4_maxhet
angsd -sites allSites_maxhet -b waltham_pop_bams_list.txt -GL 1 -P 1 $TODO -out waltham_L5_maxhet

# generating per-population SFS
realSFS boston_L1.saf.idx >boston_L1.sfs
realSFS cambridge_L2.saf.idx >cambridge_L2.sfs
realSFS arlington_L3.saf.idx >arlington_L3.sfs
realSFS belmont_L4.saf.idx >belmont_L4.sfs
realSFS waltham_L5.saf.idx >waltham_L5.sfs

# generating per-population SFS - max het filter
realSFS boston_L1_maxhet.saf.idx >boston_L1_maxhet.sfs
realSFS cambridge_L2_maxhet.saf.idx >cambridge_L2_maxhet.sfs
realSFS arlington_L3_maxhet.saf.idx >arlington_L3_maxhet.sfs
realSFS belmont_L4_maxhet.saf.idx >belmont_L4_maxhet.sfs
realSFS waltham_L5_maxhet.saf.idx >waltham_L5_maxhet.sfs

#### Pairwise comparisons

# writing down 2d-SFS priors - boston_L1 vs cambridge_L2 for maxhet
realSFS boston_L1_maxhet.saf.idx cambridge_L2_maxhet.saf.idx -P 24 > camb_bos_maxhet.sfs ; realSFS fst index boston_L1_maxhet.saf.idx cambridge_L2_maxhet.saf.idx -sfs camb_bos_maxhet.sfs -fstout camb_bos_maxhet

# global Fst between populations-- maxhet
realSFS fst stats camb_bos_maxhet.fst.idx

# Output:
	-> Assuming idxname:camb_bos_maxhet.fst.idx
	-> Assuming .fst.gz file: camb_bos_maxhet.fst.gz
	-> FST.Unweight[nObs:729076]:0.018530 Fst.Weight:0.191619
0.018530	0.191619

# writing down 2d-SFS priors - boston_L1 vs cambridge_L2 for without maxhet
realSFS boston_L1.saf.idx cambridge_L2.saf.idx -P 24 > camb_bos.sfs ; realSFS fst index boston_L1.saf.idx cambridge_L2.saf.idx -sfs camb_bos.sfs -fstout camb_bos

# global Fst between populations-- no maxhet
realSFS fst stats camb_bos.fst.idx

# Output:
	-> Assuming idxname:camb_bos.fst.idx
	-> Assuming .fst.gz file: camb_bos.fst.gz
	-> FST.Unweight[nObs:729164]:0.018536 Fst.Weight:0.189298
0.018536	0.189298

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].

# writing down 2d-SFS priors - boston_L1 vs arlington_L3--maxhet
realSFS boston_L1_maxhet.saf.idx arlington_L3_maxhet.saf.idx -P 24 > bos_arl_maxhet.sfs ; realSFS fst index boston_L1_maxhet.saf.idx arlington_L3_maxhet.saf.idx -sfs bos_arl_maxhet.sfs -fstout bos_arl_maxhet

# writing down 2d-SFS priors - boston_L1 vs arlington_L3--without maxhet
realSFS boston_L1.saf.idx arlington_L3.saf.idx -P 24 > bos_arl.sfs ; realSFS fst index boston_L1.saf.idx arlington_L3.saf.idx -sfs bos_arl.sfs -fstout bos_arl

# global Fst between populations-- maxhet
realSFS fst stats bos_arl_maxhet.fst.idx

# Output:
	-> Assuming idxname:bos_arl_maxhet.fst.idx
	-> Assuming .fst.gz file: bos_arl_maxhet.fst.gz
	-> FST.Unweight[nObs:729078]:0.029082 Fst.Weight:0.217944
0.029082	0.217944

# global Fst between populations-- no maxhet
realSFS fst stats bos_arl.fst.idx

# Output:
	-> Assuming idxname:bos_arl.fst.idx
	-> Assuming .fst.gz file: bos_arl.fst.gz
	-> FST.Unweight[nObs:729166]:0.029096 Fst.Weight:0.217091
0.029096	0.217091

# writing down 2d-SFS priors - boston_L1 vs belmont_L4--maxhet
realSFS boston_L1_maxhet.saf.idx belmont_L4_maxhet.saf.idx -P 24 > bos_bel_maxhet.sfs ; realSFS fst index boston_L1_maxhet.saf.idx belmont_L4_maxhet.saf.idx -sfs bos_bel_maxhet.sfs -fstout bos_bel_maxhet

# writing down 2d-SFS priors - boston_L1 vs belmont_L4---without maxhet
realSFS boston_L1.saf.idx belmont_L4.saf.idx -P 24 > bos_bel.sfs ; realSFS fst index boston_L1.saf.idx belmont_L4.saf.idx -sfs bos_bel.sfs -fstout bos_bel

# global Fst between populations-- maxhet
realSFS fst stats bos_bel_maxhet.fst.idx

# Output:
	-> Assuming idxname:bos_bel_maxhet.fst.idx
	-> Assuming .fst.gz file: bos_bel_maxhet.fst.gz
	-> FST.Unweight[nObs:729078]:0.019365 Fst.Weight:0.130837
0.019365	0.130837

# global Fst between populations-- no maxhet
realSFS fst stats bos_bel.fst.idx

# Output:
	-> Assuming idxname:bos_bel.fst.idx
	-> Assuming .fst.gz file: bos_bel.fst.gz
	-> FST.Unweight[nObs:729166]:0.019370 Fst.Weight:0.129753
0.019370	0.129753

# writing down 2d-SFS priors - boston_L1 vs waltham_L5--maxhet
realSFS boston_L1_maxhet.saf.idx waltham_L5_maxhet.saf.idx -P 24 > bos_walt_maxhet.sfs ; realSFS fst index boston_L1_maxhet.saf.idx waltham_L5_maxhet.saf.idx -sfs bos_walt_maxhet.sfs -fstout bos_walt_maxhet

# writing down 2d-SFS priors - boston_L1 vs belmont_L4---without maxhet
realSFS boston_L1.saf.idx waltham_L5.saf.idx -P 24 > bos_walt.sfs ; realSFS fst index boston_L1.saf.idx waltham_L5.saf.idx -sfs bos_walt.sfs -fstout bos_walt

# global Fst between populations-- maxhet
realSFS fst stats bos_walt_maxhet.fst.idx

# Output:
	-> Assuming idxname:bos_walt_maxhet.fst.idx
	-> Assuming .fst.gz file: bos_walt_maxhet.fst.gz
	-> FST.Unweight[nObs:729077]:0.016795 Fst.Weight:0.177571
0.016795	0.177571

# global Fst between populations-- no maxhet
realSFS fst stats bos_walt.fst.idx

# Output:
	-> Assuming idxname:bos_walt.fst.idx
	-> Assuming .fst.gz file: bos_walt.fst.gz
	-> FST.Unweight[nObs:729165]:0.016801 Fst.Weight:0.175617
0.016801	0.175617

# writing down 2d-SFS priors - cambridge_L2 vs arlington_L3--maxhet
realSFS cambridge_L2_maxhet.saf.idx arlington_L3_maxhet.saf.idx -P 24 > camb_arl_maxhet.sfs ; realSFS fst index cambridge_L2_maxhet.saf.idx arlington_L3_maxhet.saf.idx -sfs camb_arl_maxhet.sfs -fstout camb_arl_maxhet

# writing down 2d-SFS priors - cambridge_L2 vs arlington_L3---without maxhet
realSFS cambridge_L2.saf.idx arlington_L3.saf.idx -P 24 > camb_arl.sfs ; realSFS fst index cambridge_L2.saf.idx arlington_L3.saf.idx -sfs camb_arl.sfs -fstout camb_arl

# global Fst between populations-- maxhet
realSFS fst stats camb_arl_maxhet.fst.idx

# Output:
	-> Assuming idxname:camb_arl_maxhet.fst.idx
	-> Assuming .fst.gz file: camb_arl_maxhet.fst.gz
	-> FST.Unweight[nObs:729076]:0.134751 Fst.Weight:0.279415
0.134751	0.279415

# global Fst between populations-- no maxhet
realSFS fst stats camb_arl.fst.idx

# Output:
	-> Assuming idxname:camb_arl.fst.idx
	-> Assuming .fst.gz file: camb_arl.fst.gz
	-> FST.Unweight[nObs:729164]:0.134754 Fst.Weight:0.277617
0.134754	0.277617

# writing down 2d-SFS priors - cambridge_L2 vs belmont_L4--maxhet
realSFS cambridge_L2_maxhet.saf.idx belmont_L4_maxhet.saf.idx -P 24 > camb_bel_maxhet.sfs ; realSFS fst index cambridge_L2_maxhet.saf.idx belmont_L4_maxhet.saf.idx -sfs camb_bel_maxhet.sfs -fstout camb_bel_maxhet

# writing down 2d-SFS priors - cambridge_L2 vs belmont_L4---without maxhet
realSFS cambridge_L2.saf.idx belmont_L4.saf.idx -P 24 > camb_bel.sfs ; realSFS fst index cambridge_L2.saf.idx belmont_L4.saf.idx -sfs camb_bel.sfs -fstout camb_bel

# global Fst between populations-- maxhet
realSFS fst stats camb_bel_maxhet.fst.idx

# Output:
	-> Assuming idxname:camb_bel_maxhet.fst.idx
	-> Assuming .fst.gz file: camb_bel_maxhet.fst.gz
	-> FST.Unweight[nObs:729076]:0.022896 Fst.Weight:0.192759
0.022896	0.192759

# global Fst between populations-- no maxhet
realSFS fst stats camb_bel.fst.idx

# Output:
	-> Assuming idxname:camb_bel.fst.idx
	-> Assuming .fst.gz file: camb_bel.fst.gz
	-> FST.Unweight[nObs:729164]:0.022901 Fst.Weight:0.190546
0.022901	0.190546

# writing down 2d-SFS priors - cambridge_L2 vs waltham_L5--maxhet
realSFS cambridge_L2_maxhet.saf.idx waltham_L5_maxhet.saf.idx -P 24 > camb_walt_maxhet.sfs ; realSFS fst index cambridge_L2_maxhet.saf.idx waltham_L5_maxhet.saf.idx -sfs camb_walt_maxhet.sfs -fstout camb_walt_maxhet

# writing down 2d-SFS priors - cambridge_L2 vs waltham_L5---without maxhet
realSFS cambridge_L2.saf.idx waltham_L5.saf.idx -P 24 > camb_walt.sfs ; realSFS fst index cambridge_L2.saf.idx waltham_L5.saf.idx -sfs camb_walt.sfs -fstout camb_walt

# global Fst between populations-- maxhet
realSFS fst stats camb_walt_maxhet.fst.idx

# Output:
	-> Assuming idxname:camb_walt_maxhet.fst.idx
	-> Assuming .fst.gz file: camb_walt_maxhet.fst.gz
	-> FST.Unweight[nObs:729076]:0.010410 Fst.Weight:0.230253
0.010410	0.230253

# global Fst between populations-- no maxhet
realSFS fst stats camb_walt.fst.idx

# Output:
	-> Assuming idxname:camb_walt.fst.idx
	-> Assuming .fst.gz file: camb_walt.fst.gz
	-> FST.Unweight[nObs:729164]:0.010416 Fst.Weight:0.227082
0.010416	0.227082

# writing down 2d-SFS priors - arlington_L3 vs belmont_L4--maxhet
realSFS arlington_L3_maxhet.saf.idx belmont_L4_maxhet.saf.idx -P 24 > arl_bel_maxhet.sfs ; realSFS fst index arlington_L3_maxhet.saf.idx belmont_L4_maxhet.saf.idx -sfs arl_bel_maxhet.sfs -fstout arl_bel_maxhet

# writing down 2d-SFS priors - arlington_L3 vs belmont_L4---without maxhet
realSFS arlington_L3.saf.idx belmont_L4.saf.idx -P 24 > arl_bel.sfs ; realSFS fst index arlington_L3.saf.idx belmont_L4.saf.idx -sfs arl_bel.sfs -fstout arl_bel

# global Fst between populations-- maxhet
realSFS fst stats arl_bel_maxhet.fst.idx

# Output:
	-> Assuming idxname:arl_bel_maxhet.fst.idx
	-> Assuming .fst.gz file: arl_bel_maxhet.fst.gz
	-> FST.Unweight[nObs:729078]:0.001874 Fst.Weight:0.210582
0.001874	0.210582

# global Fst between populations-- no maxhet
realSFS fst stats arl_bel.fst.idx

# Output:
	-> Assuming idxname:arl_bel.fst.idx
	-> Assuming .fst.gz file: arl_bel.fst.gz
	-> FST.Unweight[nObs:729166]:0.001891 Fst.Weight:0.209951
0.001891	0.209951

# writing down 2d-SFS priors - arlington_L3 vs waltham_L5--maxhet
realSFS arlington_L3_maxhet.saf.idx waltham_L5_maxhet.saf.idx -P 24 > arl_walt_maxhet.sfs ; realSFS fst index arlington_L3_maxhet.saf.idx waltham_L5_maxhet.saf.idx -sfs arl_walt_maxhet.sfs -fstout arl_walt_maxhet

# writing down 2d-SFS priors - arlington_L3 vs waltham_L5---without maxhet
realSFS arlington_L3.saf.idx waltham_L5.saf.idx -P 24 > arl_walt.sfs ; realSFS fst index arlington_L3.saf.idx waltham_L5.saf.idx -sfs arl_walt.sfs -fstout arl_walt

# global Fst between populations-- maxhet
realSFS fst stats arl_walt_maxhet.fst.idx

# Output:
	-> Assuming idxname:arl_walt_maxhet.fst.idx
	-> Assuming .fst.gz file: arl_walt_maxhet.fst.gz
	-> FST.Unweight[nObs:729077]:0.067162 Fst.Weight:0.263810
0.067162	0.263810

# global Fst between populations-- no maxhet
realSFS fst stats arl_walt.fst.idx

# Output:
	-> Assuming idxname:arl_walt.fst.idx
	-> Assuming .fst.gz file: arl_walt.fst.gz
	-> FST.Unweight[nObs:729165]:0.067172 Fst.Weight:0.261947
0.067172	0.261947

# writing down 2d-SFS priors - belmont_L4 vs waltham_L5--maxhet
realSFS belmont_L4_maxhet.saf.idx waltham_L5_maxhet.saf.idx -P 24 > bel_walt_maxhet.sfs ; realSFS fst index belmont_L4_maxhet.saf.idx waltham_L5_maxhet.saf.idx -sfs bel_walt_maxhet.sfs -fstout bel_walt_maxhet

# writing down 2d-SFS priors - belmont_L4 vs waltham_L5---without maxhet
realSFS belmont_L4.saf.idx waltham_L5.saf.idx -P 24 > bel_walt.sfs ; realSFS fst index belmont_L4.saf.idx waltham_L5.saf.idx -sfs bel_walt.sfs -fstout bel_walt

# global Fst between populations-- maxhet
realSFS fst stats bel_walt_maxhet.fst.idx

# Output:
	-> Assuming idxname:bel_walt_maxhet.fst.idx
	-> Assuming .fst.gz file: bel_walt_maxhet.fst.gz
	-> FST.Unweight[nObs:729077]:0.026733 Fst.Weight:0.175961
0.026733	0.175961

# global Fst between populations-- no maxhet
realSFS fst stats bel_walt.fst.idx

# Output:
	-> Assuming idxname:bel_walt.fst.idx
	-> Assuming .fst.gz file: bel_walt.fst.gz
	-> FST.Unweight[nObs:729165]:0.026739 Fst.Weight:0.174107
0.026739	0.174107


#------------------------------NGSRelate-- Estimating pairwise relatedness, inbreeding coefficients, kinship coefficient theta

# First make new directory called ngsRelate

pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/ngsRelate

# Inside are two more dirs called global_samples and mass_samples containing all the filtered bams (284) and the file lists (for global) and all mass bams and file list (for mass)

ls -lh
global_samples
mass_samples

# Download and install from https://github.com/ANGSD/NgsRelate

git clone --recursive https://github.com/SAMtools/htslib
git clone https://github.com/ANGSD/ngsRelate
cd htslib/;make -j2;cd ../ngsRelate;make HTSSRC=../htslib/

####---- Let's start with the global samples using the default filters

# First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).

# Make a qsub for this, it's too big, called allele_freqs_ngsrelate.qsub

module load angsd
angsd -b alldups_12p_bams_nodups -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3

### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat angsdput.mafs.gz | cut -f5 |sed 1d >freq

### run NgsRelate
# Make a qsub for this too called ngsrelate.qsub

../ngsRelate/ngsRelate  -g angsdput.glf.gz -n 284 -f freq  -O newres_global

# changing the name of angsdput.glf.gz to angsdput_12p.glf.gz to avoid confusion

### Doing again for the 50p dataset
# Changing the allele_freqs_ngsrelate.qsub to be this:

module load angsd
angsd -b alldups_50p_bams_nodups -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3

### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs). "angsdput.mafs.gz" is the output for the 50p run
zcat angsdput.mafs.gz | cut -f5 |sed 1d >freq_50p

### run NgsRelate for 50p
# Make a qsub for this too called ngsrelate.qsub

../ngsRelate/ngsRelate  -g angsdput.glf.gz -n 281 -f freq_50p  -O newres_global_50p

# Everything from this run has "angsdput" as the prefix and it made a file called newres (changed on local machine to newres_global_12p)
# Look at the file called newres_global_12p for everything

####---- Now going to do mass samples also with default filters
cd ../mass_samples/

pwd
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/ngsRelate/mass_samples

# Let's start with the global samples using the default filters

# First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).

# Make a qsub for this, it's too big, called allele_freqs_ngsrelate.qsub

module load angsd
angsd -b mass_bams -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3

### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat angsdput.mafs.gz | cut -f5 |sed 1d >freq

### run NgsRelate
# Make a qsub for this too called ngsrelate.qsub

../ngsRelate/ngsRelate  -g angsdput.glf.gz -n 84 -f freq  -O newres_mass

# Everything from this run has "angsdput" as the prefix and it made a file called newres_mass
# Look at the file called newres_mass for everything


# Now run NGSrelate one more time with a technical replicate added in for relatedness visualization. This tech rep is:
IAN_2brad_L3_3_S3_L001_R1_001_AGAC.trim.bt2.bam

# Make a new config file called alldups_12p_bams_nodups_1techrep and put the tech rep between

IAN_2brad_L3_3_S3_L001_R1_001_ACCA.trim.bt2.bam

# And

IAN_2brad_L3_3_S3_L001_R1_001_AGTG.trim.bt2.bam

# Copy it over from 
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/ANGSD/bams_dir

# Run this qsub with tech rep to make new freq file:
allele_freqs_ngsrelate.qsub

# With these parameters:

module load angsd

angsd -b alldups_12p_bams_nodups_1techrep -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3


### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat angsdput.mafs.gz | cut -f5 |sed 1d >freq_techreps

# Run NGSrelate one more time with these parameters:
../ngsRelate/ngsRelate  -g angsdput.glf.gz -n 285 -f freq_techreps  -O newres_techreps

#------------------------ Treemix ML tree of admixture events

# First going to run through the tutorial from the Speciation and Population Genomics physalia course: https://speciationgenomics.github.io/Treemix/, this specific tutorial is found here: https://speciationgenomics.github.io/Treemix/

# Need to create a conda environment called treemix to do some of this stuff in and to load certain modules up

# ---------- Set Up Conda Environment ---------- #

module load miniconda/24.5.0

conda config --add channels defaults
# Warning: 'defaults' already in 'channels' list, moving to the top
conda config --add channels bioconda
conda config --add channels conda-forge

# create new conda environment for treemix stuff
conda create --name treemix

# Where is was installed:
 environment location: /projectnb/mullenl/inovick/.conda/envs/treemix

# activate conda environment
conda activate treemix

# install packages necessary for genome assembly
# for example, install the program "treemix"
conda install bioconda::treemix

# When done, to deactivate conda environment:

conda deactivate

# First need to ldprune to get rid of samples with high linkage
# Going to use this script found here: https://github.com/speciationgenomics/scripts/blob/master/ldPruning.sh

# Input file has to be vcf.gz so need to get my vcf file from 
/projectnb/mullenl/novick/global_pop_gen/rad_workflow/ANGSD/

# Using the minindc 91.5%, copied from when I did raxml on this earlier, found in:

/projectnb/mullenl/novick/global_pop_gen/rad_workflow/raxml/myresult_alldups_12p_nodups.vcf

# Copy over to /projectnb/mullenl/novick/global_pop_gen/rad_workflow/treemix

cd /projectnb/mullenl/novick/global_pop_gen/rad_workflow/treemix

cp /projectnb/mullenl/novick/global_pop_gen/rad_workflow/raxml/myresult_alldups_12p_nodups.vcf ./

# Make qsub for ldpruning called ldpruning.qsub

#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -o ldpruning.qlog #Here is the logfile
#$ -N ldpruning # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M inovick@bu.edu #your email
#$ -m be
#$ -pe omp 4

module load plink2/2.0

# perform linkage pruning - i.e. identify prune sites
plink --vcf myresult_alldups_12p_nodups.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out ldpruned_12p_nodups

# This makes 2 files, one called 
ldpruned_12p_nodups.prune.in

# and
ldpruned_12p_nodups.prune.out

# We want the ldpruned_12p_nodups.prune.in file because this contains the sites that fell below the linkage threshold (those we should retain)

# Also going to linkage prune again using this script to get the right usable input files for treemix:
module load vcftools
vcftools --gzvcf myresult_alldups_12p_nodups.vcf --max-missing 1 --recode --stdout | gzip > myresult_alldups_12p_nodups.vcf.gz

module load plink
ldPruning.sh myresult_alldups_12p_nodups.vcf.gz
# Makes this: myresult_alldups_12p_nodups.LDpruned.vcf.gz

# This should be the input file for later, the myresult_alldups_12p_nodups.vcf.gz


# Going to generate the clust file
bcftools query -l myresult_alldups_12p_nodups.LDpruned.vcf.gz | awk '{split($1,pop,"."); print $1"\t"$1"\t"pop[2]}' > moths.clust

# Run the conversion script
vcf2treemix.sh myresult_alldups_12p_nodups.LDpruned.vcf.gz moths.clust

# Trying to remake the clust file so it shows lineages as countries in the names and a separate column for lineages, following /Users/izzynovick/Desktop/Tineoid_Moth_Project/Chapter_2_RAD/Migration/Treemix/clust_file_gen.R
# New clust file called moths.clust

# Rerun this:
vcf2treemix.sh myresult_alldups_12p_nodups.LDpruned.vcf.gz moths.clust

# It worked!

# Now to run treemix. Going to use 0-11 potential migration edges because of the number of countries, bootstrapping 500 replicates

for i in {0..11}
do
 treemix -i myresult_alldups_12p_nodups.LDpruned.treemix.frq.gz -m $i -o moth_output_treemix.$i -bootstrap -k 500 > treemix_${i}_log &
done

# Now go here to treemix data page: https://bitbucket.org/nygcresearch/treemix/downloads/
# And download treemix-1.13.tar.gz, go to src, and download plotting_funcs.R

# This resulted in 11 "treemix_[number]_log" files that are all empty. After emailing the makers of the treemix course I've been using, and being told this could be because I have multiallelic sites or indels, I checked again and don't have that. Not sure why this is not running

# This could be because it thinks I have 12 populations. Going to change it to 0...10

for i in {0..10}
do
 treemix -i myresult_alldups_12p_nodups.LDpruned.treemix.frq.gz -m $i -o moth_output_treemix.$i -bootstrap -k 500 > treemix_${i}_log &
done

# Still didn't work. Going to try and edit the moths.clust file and rerun conversion script
# Making a new clust file called moths2.clust where it's just going to be 2 columns instead of 3
# Also copied the myresult_alldups_12p_nodups.LDpruned.vcf.gz file to be called 

myresult_alldups_12p_nodups_2.LDpruned.vcf.gz

# Now going to run this 
vcf2treemix.sh myresult_alldups_12p_nodups_2.LDpruned.vcf.gz moths2.clust





# Got this series of errors:
Unrecognized values used for CHROM: pseudochr_OX411273.1_29 - Replacing with 0.
Done.
After filtering, kept 539 out of a possible 539 Sites
Run Time = 3.00 seconds
PLINK v1.90b6.27 64-bit (10 Dec 2022)          www.cog-genomics.org/plink/1.9/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to myresult_alldups_12p_nodups_2.LDpruned.log.
Options in effect:
  --allow-extra-chr 0
  --allow-no-sex
  --file myresult_alldups_12p_nodups_2.LDpruned
  --make-bed
  --out myresult_alldups_12p_nodups_2.LDpruned

256039 MB RAM detected; reserving 128019 MB for main workspace.
Allocated 17087 MB successfully, after larger attempt(s) failed.
.ped scan complete (for binary autoconversion).
Performing single-pass .bed write (539 variants, 284 people).
--file: myresult_alldups_12p_nodups_2.LDpruned-temporary.bed +
myresult_alldups_12p_nodups_2.LDpruned-temporary.bim +
myresult_alldups_12p_nodups_2.LDpruned-temporary.fam written.
539 variants loaded from .bim file.
284 people (0 males, 0 females, 284 ambiguous) loaded from .fam.
Ambiguous sex IDs written to myresult_alldups_12p_nodups_2.LDpruned.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 284 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is exactly 1.
539 variants and 284 people pass filters and QC.
Note: No phenotypes present.
--make-bed to myresult_alldups_12p_nodups_2.LDpruned.bed +
myresult_alldups_12p_nodups_2.LDpruned.bim +
myresult_alldups_12p_nodups_2.LDpruned.fam ... done.
PLINK v1.90b6.27 64-bit (10 Dec 2022)          www.cog-genomics.org/plink/1.9/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to myresult_alldups_12p_nodups_2.LDpruned.log.
Options in effect:
  --allow-extra-chr 0
  --allow-no-sex
  --bfile myresult_alldups_12p_nodups_2.LDpruned
  --freq
  --missing
  --out myresult_alldups_12p_nodups_2.LDpruned
  --within moths2.clust

256039 MB RAM detected; reserving 128019 MB for main workspace.
Allocated 17087 MB successfully, after larger attempt(s) failed.
539 variants loaded from .bim file.
284 people (0 males, 0 females, 284 ambiguous) loaded from .fam.
Ambiguous sex IDs written to myresult_alldups_12p_nodups_2.LDpruned.nosex .
Warning: No samples named in --within file remain in the current analysis.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 284 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is exactly 1.
--freq: Allele frequencies (founders only) written to
myresult_alldups_12p_nodups_2.LDpruned.frq .
--missing: Sample missing data report written to
myresult_alldups_12p_nodups_2.LDpruned.imiss, and variant-based missing data
report written to myresult_alldups_12p_nodups_2.LDpruned.lmiss.
gzip: myresult_alldups_12p_nodups_2.LDpruned.frq.strat: No such file or directory
Traceback (most recent call last):
  File "./plink2treemix.py", line 9, in <module>
    infile = gzip.open(sys.argv[1])
  File "/usr/lib64/python2.7/gzip.py", line 34, in open
    return GzipFile(filename, mode, compresslevel)
  File "/usr/lib64/python2.7/gzip.py", line 94, in __init__
    fileobj = self.myfileobj = __builtin__.open(filename, mode or 'rb')
IOError: [Errno 2] No such file or directory: 'myresult_alldups_12p_nodups_2.LDpruned.frq.strat.gz'
gzip: myresult_alldups_12p_nodups_2.LDpruned.treemix.frq.gz: No such file or directory
gzip: myresult_alldups_12p_nodups_2.LDpruned.frq.strat.gz: No such file or directory
paste: myresult_alldups_12p_nodups_2.LDpruned.treemix.frq: No such file or directory
gzip: myresult_alldups_12p_nodups_2.LDpruned.treemix.frq: No such file or directory

# Still not fully working. I need the .frq.strat file so I'm going to try this:
plink --vcf myresult_alldups_12p_nodups_2.LDpruned.vcf.gz --freq --within moths.clust --out myresult_alldups_12p_nodups_2.LDpruned.strat --double-id --allow-extra-chr

# I'm going to make a qsub for this and try it, maybe it's too much ram
# It's called
vcf2treemix.qsub 

# It didn't work, now I'm going to actually zip the freq.strat file instead of being an idiot
# It still didn't work and I've been continually getting this error: Error: All --within Categories Are Null
# I'm going to try and compare names from vcf file to names in clust file and see if they are the same?? Which they should be?

# Extract sample names from vcf file:
module load htslib/1.16
module load bcftools
bcftools query -l myresult_alldups_12p_nodups_2.LDpruned.vcf.gz > vcf_samples.txt

# Compare sample names by making a new file with the names that match and then counting the lines in that file
grep -F -f vcf_samples.txt moths2.clust > matched_samples.txt

# They have the same number

# Goingto see if syntax etc is the same using
awk 'NR==FNR{samples[$1]; next} $1 in samples' vcf_samples.txt moths2.clust

# They all match
# Going to try and get rid of trailing spaces in the clust file by making a new clust file 

sed 's/[[:space:]]*$//' moths2.clust > moths2_clean.clust

# Going to try to linkage prune again and make new input files
module load vcftools
vcftools --gzvcf myresult_alldups_12p_nodups.vcf --max-missing 1 --recode --stdout | gzip > myresult_alldups_12p_nodups_3.vcf.gz

module load plink
ldPruning.sh myresult_alldups_12p_nodups_3.vcf.gz

# Going to run the conversion script again now
# Didn't work

# Going to edit the myresult_alldups_12p_nodups_3.LDpruned.vcf and change the header and the names of the chromosomes to try and match the dog one
# Wait I realized I need to redo the linkage pruning so I made a new copy of the original vcf called 
myresult_alldups_12p_nodups_edited.vcf 

# Editing the header and going to change the names of the chromosomes to 1, 2, etc
sed -i 's/pseudochr_OX411244.1_Z/30/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 30 even tho its chromosome z
sed -i 's/pseudochr_OX411245.1_1/1/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 1
sed -i 's/pseudochr_OX411246.1_2/2/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 2
sed -i 's/pseduochr_OX411247.1_3/3/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 3
sed -i 's/pseudochr_OX411248.1_4/4/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 4
sed -i 's/pseudochr_OX411249.1_5/5/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 5
sed -i 's/pseudochr_OX411250.1_6/6/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 6
sed -i 's/pseudochr_OX411251.1_7/7/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 7
sed -i 's/pseudochr_OX411252.1_8/8/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 8
sed -i 's/pseudochr_OX411253.1_9/9/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 9
sed -i 's/pseudochr_OX411254.1_10/10/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 10
sed -i 's/pseudochr_OX411255.1_11/11/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 11
sed -i 's/pseudochr_OX411256.1_12/12/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 12
sed -i 's/pseudochr_OX411257.1_13/13/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 13
sed -i 's/pseudochr_OX411258.1_14/14/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 14
sed -i 's/pseudochr_OX411259.1_15/15/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 15
sed -i 's/pseudochr_OX411260.1_16/16/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 16
sed -i 's/pseudochr_OX411261.1_17/17/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 17
sed -i 's/pseudochr_OX411262.1_18/18/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 18
sed -i 's/pseudochr_OX411263.1_19/19/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 19
sed -i 's/pseudochr_OX411264.1_20/20/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 20
sed -i 's/pseudochr_OX411265.1_21/21/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 21
sed -i 's/pseudochr_OX411266.1_22/22/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 22
sed -i 's/pseudochr_OX411267.1_23/23/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 23
sed -i 's/pseudochr_OX411268.1_24/24/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 24
sed -i 's/pseudochr_OX411269.1_25/25/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 25
sed -i 's/pseudochr_OX411270.1_26/26/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 26
sed -i 's/pseudochr_OX411271.1_27/27/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 27
sed -i 's/pseudochr_OX411272.1_28/28/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 28
sed -i 's/pseudochr_OX411273.1_29/29/g' myresult_alldups_12p_nodups_edited.vcf # replace this chrom with 29
sed -i 's/z/30/g' myresult_alldups_12p_nodups_edited.vcf
# Going to try to linkage prune again and make new input files

module load plink2/2.0

# perform linkage pruning - i.e. identify prune sites (method 1)
plink --vcf myresult_alldups_12p_nodups_edited.vcf --double-id --allow-extra-chr --chr-set 30 no-xy no-mt\
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out edited_ldpruned_12p_nodups

# I don't think it did what I wanted it to do so I'm going to try this method now (method 2):
# Making a new qsub for it to do this called ldpruning_2.qsub
module load vcftools
vcftools --gzvcf myresult_alldups_12p_nodups_edited.vcf --max-missing 1 --recode --stdout | gzip > myresult_alldups_12p_nodups_edited.vcf.gz

module load plink
ldPruning.sh myresult_alldups_12p_nodups_edited.vcf.gz

# ok i think the first method did work b ut it results in a .prune.in file like so (this is the one I used with the chr-set flag):
edited_3_ldpruned_12p_nodups.prune.in
# and i'm not sure how to continue with the pipeline with  this

###### Estimated Effective Migration Surface (eems)
# Looking at the documentation here: https://github.com/dipetkov/eems/blob/master/Documentation/EEMS-doc.pdf
# First download eigen- go to software directory 

cd /projectnb/mullenl/software
wget https://gitlab.com/libeigen/eigen/-/archive/3.2.10/eigen-3.2.10.tar.gz
tar -xvzf eigen-3.2.10.tar.gz

# Then install boost in same directory
wget http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.gz
tar -xvzf boost_1_57_0.tar.gz

# Make a new directory called boost, copy path
# Then go into the boost_1_57 directory and type this:

module load python2/2.7.13
./bootstrap.sh --prefix=/projectnb/mullenl/software/boost
./b2 install

# Now your software/boost/ directory should have a lib/ and an include/ directory
# Made a new directory in /projectnb/mullenl/novick/global_pop_gen/rad_workflow/ called eems
git clone https://github.com/dipetkov/eems.git
cd eems/runeems_snps/src 

# Edit Makefile as follows: 
EIGEN_INC = /projectnb/mullenl/software/eigen-3.2.10
BOOST_LIB = /projectnb/mullenl/software/boost/lib
BOOST_INC = /projectnb/mullenl/software/boost/include

# Exit the Makefile, then run:
make linux

# Need .bed, .bim and .fam files from VCF
# Use vcf tools:
vcftools --vcf my_data.vcf --out my_data --plink

# or plink:
plink --vcf myfile.vcf.gz --recode --out myfile

