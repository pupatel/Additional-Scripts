#Created by Parth Patel, DBI @ University of Delaware, Newark, Delaware 19717
#Date created: 12/12/2017 

#This is an example script for processing PACBIO ISOSEQ PROTOCOL

# Step 1. Create temporary directory for output
#BASEDIR=$(dirname "$0")
#echo "$BASEDIR"
#mkdir /home/parth/PACBIO

cd /efs/home/pupate/CEW
# Step 2. Download bax.h5 files from S3

aws s3 cp ‐‐sse AES256 s3://Analysis_Results/mxxxx_1_xxxx__bax.h5
aws s3 cp ‐‐sse AES256 s3://Analysis_Results/mxxxx_2_xxxx__bax.h5
aws s3 cp ‐‐sse AES256 s3://genome‐analytics‐sequencing‐legacy/pacbio/2015/Run349_578/A04_1/Analysis_Results/m150217_021725_00117_c100776102550000001823164607301524_bax.h5


# Step 3. For each movie, use 3 bax.h5 files to create bam file containing subreads

bax2bam ‐o mxxxx_1_xxxx__bax.h5 mxxxx_2_xxxx.bam
bax2bam ‐o mxxxx_2_xxxx__bax.h5 mxxxx_2_xxxx.bam

# Step 4. CCS STEP: For the subreads from each movie, build circular consensus sequences
ccs ‐‐numThreads 15 ‐‐minLength=300 ‐‐minPasses=0 ‐‐minZScore=‐999 ‐‐maxDropFraction=0.8 ‐‐minPredictedAccuracy=0.75 ‐‐minSnr=3.75 ccs_1.bam
ccs ‐‐numThreads 15 ‐‐minLength=300 ‐‐minPasses=0 ‐‐minZScore=‐999 ‐‐maxDropFraction=0.8 ‐‐minPredictedAccuracy=0.75 ‐‐minSnr=3.75 ccs_2.bam

# Step 5. Define subread and circular consensus data sets containing files associated with the 4 BBM movies
dataset create ‐‐type ConsensusReadSet ccs.xml ccs_1.bam ccs_2.bam
dataset create ‐‐type SubreadSet ccs.subreads.xml ccs_1.subreads.bam ccs_2.subreads.bam


# Step 6. CLASSIFY STEP:
# ‐ classify sequences as full length/non full length.
# ‐ identify and remove polyA/T tails, remove primers, and identify read strandedness
pbtranscript classify ‐‐flnc isoseq_flnc.fasta ‐‐nfl isoseq_nfl.fasta ‐d classifyOut ‐‐cpus 15 ‐‐min_seq_len 100 ccs.xml isoseq.pacbio.fasta


# Step 7. Seperate reads into length bins
separate_flnc.py ‐‐bin_size_kb 2 ‐‐max_base_limit 300 isoseq_flnc.fasta clusterOut out.pickle


# Step 9. Upload results to S3
aws s3 cp ‐‐sse AES256 isoseq.pacbio.fasta s3://Output


# Step 10. Collapse similar sequences based on % ID using cdhit
cdhit‐est ‐i pacbio.fasta ‐o pacbio.CDHIT100.fasta ‐c 1.0 ‐T 15 ‐G 0 ‐aL 0.90 ‐AL 100 ‐aS 0.99 ‐AS 30
