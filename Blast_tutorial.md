# This document describes how to use the Dacini COI Database to identify species identity using Maximum likelihood, BLASTn and Bold searches.

## How to install programs

We use a BLAST and BOLD approaches to molecularly identify sequences using both local and public databases and then we use a phylogenetic approach using Iqtree's fast phylogenetic tree searching algorithm

This pipeline requires installing several freely available software packages and was tested in OSX but should work in Linux systems. 

Several of the output manipulations steps may require some modification to run on Windows but should run on Windows subsystem linux.

The easiest way to install these programs is to use first install anaconda on to your OS.

Windows  (https://docs.anaconda.com/anaconda/install/windows/)
Mac OSX  (https://docs.anaconda.com/anaconda/install/mac-os/)
Linux    (https://docs.anaconda.com/anaconda/install/linux/)

Follow instructions and once finished we will use anaconda and make a new environment to install the following programs.
```
$ conda env -n IQ_blast anaconda
$ conda activate IQ_blast
$ conda install -c bioconda iqtree
$ conda install -c bioconda muscle
$ conda install -c bioconda blast
$ conda install -c bioconda seqtk
$ conda install -c bioconda ncbi-acc-download
```
This should install muscle, blast, seqtk, iqtree and ncbi-acc-download into your path in the conda environment IQ_blast.

You will need to use software to visualize your phylogeny. We recommend Figtree (http://tree.bio.ed.ac.uk/software/figtree/)

This markdown will show the step by step process of how to molecularly identify Dacini specimens using a highly curated local Dacini Database and downloaded NCBI database.

We have also included a script "Script.txt" which will do all the steps in one command 
```
sh Script.txt Unique_haplotypes_COI_Dacini.fasta good_larvae_updated.fasta NCBI_ref_brachycera_COI.fasta
```
Dacini_ID_pipeline.sh contains script 
database.fasta contains all sequences to be included in local database
sequence_of_interest.fasta contains all sequences we want to identify
NCBI_ref_brachycera_COI.fasta contains all COI sequences of brachycera from genebank

This script will out put several files. These are some of the important files for molecularly identifying sequences of interest.

```
Local_DB_top_hit_Seq_ID.txt
```
Contains tab formated table with the top Blast hit for each sequence from the local blast search
| query ID  | Subject ID  | evalue  | Percent identical   |
|---|---|---|---|
| BX141023-001  |  ms07466_Bactrocera_occipitalis_FFxxMa002ME_MYS_COI_1493 | 0  | 99.665  |

```
Local_DB_uncertain_ID.txt
```
Contains tab formated table with the top Blast hit for sequences with less than 97% P-ID from the local blast search.

```
NCBI_BLAST_uncertain_top_hit_Seq_ID.out
```
Contains tab formated table with the top Blast hit for each sequence from the Downloaded NCBI database blast search 

```
NCBI_BLAST_uncertain_Seq_ID.txt
```
Contains tab formated table with the top Blast hit for sequences with less than 97% P-ID from the Downloaded NCBI database blast search.


```
NCBI_Blast_uncertain_ID_genebank_accessions.txt
```
Contains the extracted genebank accessions from `NCBI_BLAST_uncertain_Seq_ID.txt`

```
NCBI_Blast_uncertain_genebank_info.txt
```
Contains the header information from the genebank accessions for which had a top blast hit.

BOLD is not as computer friendly as NCBI. 

You will need to copy and paste the fasta file `NCBI_Blast_uncertain_ID.fasta` onto (http://www.boldsystems.org/index.php/IDS_OpenIdEngine) to use the BOLD database.

You will also need to have an account to lookup large numbers of sequences but having an account will allow you to get email results for your search.

With this information you can use these information to now molecularly identify your sequences. 

## Step by step process

## STEP ONE local blast searches

## Create local databases
```
makeblastdb -in Unique_haplotypes_COI_Dacini.fas -dbtype nucl  -out Unique_haplotypes_COI_Dacini_DB -title "Unique_haplotypes_COI_Dacini_DB"
```
This command makes local blast db from a fasta file within the folder.

## Search highly curated local database to conduct sequence id

```
blastn -task blastn -db Unique_haplotypes_COI_Dacini_DB -query good_larvae_updated.fasta -dust no -max_target_seqs 5 -outfmt "6 qseqid sseqid evalue pident " -out Seq_ID.txt -num_threads 4
```

This command searches local database using BLASTn for sequences within the good_larvae fasta outputing the top 5 hits in a tabular table outputed as Seq_ID.txt

The table will have 4 columns example below

| query ID  | Subject ID  | evalue  | Percent identical   |
|---|---|---|---|
| BX141023-001  |  ms07466_Bactrocera_occipitalis_FFxxMa002ME_MYS_COI_1493 | 0  | 99.665  |



## Sort BLASTn output and export only top hit by highest percent identical

Due to the flag ```-max_target_seqs 5 ``` we will have up to 5 hits per query.

We will sort and pull only the top hit using command below.

```
sort -k1,1 -k4,4nr Seq_ID.txt | sort -u -k1,1 > top_hit_Seq_ID.txt
```
This first sorts by query ID (k1) and then sorts in desending order (Percent identical). 

This output is piped to another 'sort' which takes the first unique query ID (-u k1,1) in the tabular txt file and outputs it into top_hits_seq_ID.txt

We will use awk to extract sequences with pident less than 97%

```
awk '$4 < 97' top_hit_Seq_ID.txt > uncertain_ID.txt
awk '$4 < 97 {print $1}' top_hit_Seq_ID.txt > uncertain_ID_voucher.txt
```
These commands outputs 2 text files. One is a table with only sequences with lest than 97% P-ID
and a list of these sequences.

We then use seqtk to extract sequences from these samples from the input fasta "good_larvae_updated.fasta" 
```
seqtk subseq good_larvae_updated.fasta uncertain_ID_voucher.txt > uncertain_ID.fasta
```
Now we have a fasta containing only sequences that we have uncertain molecular ID so we use remote blastn  against the NCBI NR database to see if the NCBI database has matching sequences.
```
blastn -db nt -query uncertain_ID.fasta -taxids 7147 -max_target_seqs 5 -outfmt "6 qseqid sseqid evalue pident" -out results_taxid.out -remote 

```
This may take a while due to low priority

If you have many uncertain sequences this can take a long time. 

It is faster to download all COI sequences from genebank for your taxa of interest. 

To do this go to NCBI(https://www.ncbi.nlm.nih.gov/nuccore) and search Brachycera AND (COI OR CO1 OR COX1 OR COXI) restricting length to 200-2000bp.

This should download +400000 sequences. 

Once downloaded you will need to format NCBI database to search 
```
makeblastdb -in NCBI_ref_brachycera_COI.fasta -dbtype nucl  -out NCBI_ref_brachycera_COI_DB -title "NCBI_ref_brachycera_COI_DB"

```

alternative is to run them one at a time using web based blast

Processing results from NCBI Blastn search should follow same steps as your highly curated local DB search 
```
sort -k1,1 -k4,4nr NCBI_uncertain_Seq_ID.txt | sort -u -k1,1 > NCBI_BLAST_uncertain_top_hit_Seq_ID.out
awk '$4 < 97' NCBI_BLAST_uncertain_top_hit_Seq_ID.out > NCBI_BLAST_uncertain_Seq_ID.txt
awk '$4 < 97 {print $1}' NCBI_BLAST_uncertain_top_hit_Seq_ID.out > NCBI_Blast_uncertain_ID_voucher.txt
awk '{print $2}' NCBI_BLAST_uncertain_top_hit_Seq_ID.out | sort | unique > NCBI_Blast_uncertain_ID_genebank_accessions.txt
seqtk subseq $2 NCBI_Blast_uncertain_ID_voucher.txt > NCBI_Blast_uncertain_ID.fasta
```

Then we download Genebank fastas for all uncertain ID's from the NCBI database to get information

```
file="NCBI_Blast_uncertain_genebank.txt"
while read line; do 
    ncbi-acc-download -m nucleotide --out /dev/stdout -F fasta $line >> NCBI_Blast_uncertain_genebank.fasta
done < $file
```
extract all header information

```
grep -e ">" NCBI_Blast_uncertain_genebank.fasta | awk 'sub(/^>/, "")' > NCBI_Blast_uncertain_genebank_info.txt
```

If you have a BOLD account you can copy an paste the fasta sequences from NCBI_Blast_uncertain_ID.fasta into bold to get bold ID's and compare to Genebank

We have run through the BLAST/BOLD searches and we may still have some sequences that do not have ID's.

We first need to align samples to be identified with your database by making  a phylogenetic tree.

Your data needs to align you database with your sequences of interest.

Both the sequences of interest and database files should be in fasta format.

```
cat Unique_haplotypes_COI_Dacini.fasta uncertain_ID.fasta > seqs.fasta

muscle -in seqs.fasta -out seqs_aligned.fasta
```
This command will run a fast ML tree search using GTR model.
```
iqtree -s seqs_aligned.fasta -fast -m GTR 
```

Tree can be viewed using figtree.










