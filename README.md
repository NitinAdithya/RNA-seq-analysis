# RNA-seq-analysis
RNA seq analysis using Python and R

We will be using the following data:\
Reference article: [ONECUT2 is a driver of neuroendocrine prostate cancer - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336817/#Sec11)\
Data downloaded from: [GEO Accession viewer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305)\
Github workflow followed: [Bulk RNA-sequencing pipeline and differential gene expression analysis](https://erilu.github.io/bulk-rnaseq-analysis/#Obtaining_raw_data_from_GEO)

The overall workflow which we will follow is as follows:


First, we installed the SRA Toolkit, which provides us the tools 'prefetch' and 'fastq-dump' to download and convert SRA files. 
Next, we downloaded the raw sequence data for each SRR ID using 'prefetch'. This saved the '.sra' files in the default SRA folder and then converted these '.sra' files into compressed FASTQ files using 'fastq-dump' with options that handled paired-end reads (though our data is single-ended) and ensured quality filtering.

To save time and automate this for multiple SRR IDs, I wrote a Python script (which you modified/improved) that sequentially downloads and converts all SRA files to FASTQ format.

Once we had the FASTQ files, we performed quality control using FastQC, which gave us individual reports for each file and then used MultiQC to consolidate these reports into a single HTML file, which helped us decide if trimming was necessary.

For demonstration, we performed optional trimming of one sample using Trimmomatic and then ran FastQC again to confirm that the read quality improved.

Finally, to make our lives easier, we carefully concatenated and renamed the 20 FASTQ files to just 8. After verifying that the files were correct, I safely removed the original '.sra' files. I repeat after verifying that the files were correct...

Suggestion: When sharing results on GitHub, we would just upload processed counts matrices, metadata, and scripts instead of the raw FASTQ files, due to their large size (this is just for those interested in getting registration fee back after making GitHub bulk rnaseq repository worth reproducing!!)
---
Got it! Let me summarize and give you a careful, structured workflow based on your steps so far, along with some tips and safety checks. I’ll also highlight where things could go wrong and how to avoid data loss.

---


### Commands

'''
sudo apt update
sudo apt install sra-toolkit   ##Install SRA Toolkit:

prefetch SRR7179504   ## Download SRA files:

fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR7179504.sra  ## Convert SRA to FASTQ: change paths accordingly

'''

'''
python3 fastq_download.py  ## Automate multiple SRR downloads with Python: You have the .py in email. This makes 20 fastq.gz files, each generated in its own fastq/ directory.
'''

'''
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8  ## FastQC for quality control
multiqc fastqc_results/ -o multiqc_report/   ## MultiQC to get one HTML report. Use it to decide if trimming is necessary.
'''

'''
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 fastq/SRR7179504.fastq.gz fastq/SRR7179504_trimmed.fastq.gz TRAILING:10 -phred33   ## Optional trimming with Trimmomatic: and so try just for any one fastq
'''
After trimming, run FastQC again to compare metrics.

'''
## Concatenate and rename FASTQ files carefully: LNCAP Normoxia replicates | # LNCAP Hypoxia replicates | # PC3 samples (just rename) --> this now makes just 8 files from those 20 files.
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz > LNCAP_Normoxia_S1.fastq.gz
cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz > LNCAP_Normoxia_S2.fastq.gz
cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz > LNCAP_Hypoxia_S1.fastq.gz
cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz > LNCAP_Hypoxia_S2.fastq.gz
mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz
'''

'''
rm -rf SRR*    ## Be extremely careful to avoid overwriting files. Double-check filenames before running cat or mv. Use ls fastq/ to ensure you have the right SRR IDs.
'''

What next??
Now, we will proceed with downloading the HISAT2 prebuilt genome index for the human genome (GRCh38) using wget. This retrieves the compressed tar file containing the genome index that HISAT2 requires for alignment. We will then extracted the tar file using tar -xvzf, which creates the directory containing all necessary genome index files.

```
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz
```

Next, we install HISAT2 and Samtools to do alignment and downstream processing. HISAT2 is used for aligning RNA-seq reads to a reference genome, and Samtools is used to sort and index the resulting BAM files.

```
sudo apt install hisat2
sudo apt install samtools
```

For aligning reads, we run HISAT2 on each FASTQ file. The command uses the -q flag to indicate FASTQ input, -x to specify the genome index directory, and -U for unpaired reads. The alignment output is piped to Samtools to sort the reads into a BAM file using samtools sort, and then indexed with samtools index for efficient access and downstream analysis.

```
hisat2 -q -x grch38/genome -U fastq/sample1.fastq.gz | samtools sort -o sample1.bam | samtools index sample1.bam  ## try this (will take like 20-25 mins)
---
