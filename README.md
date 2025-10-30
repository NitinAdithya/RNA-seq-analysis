# RNA-seq-analysis
RNA seq analysis using Python and R

We will be using the following data:\
Reference article: [ONECUT2 is a driver of neuroendocrine prostate cancer - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336817/#Sec11)\
Data downloaded from: [GEO Accession viewer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305)\
Github workflow followed: [Bulk RNA-sequencing pipeline and differential gene expression analysis](https://erilu.github.io/bulk-rnaseq-analysis/#Obtaining_raw_data_from_GEO)

The overall workflow which we will follow is as follows:
1. Download the raw data from GEO using SRA Toolkit ('prefetch', 'fastqdump')
2. Quality control with FastQC and MultiQC
3. Trimming (optional) with Trimmomatic
4. Alignment with HISAT2
5. Quantification with feature counts
6. Differential gene expression analysis in R (DESeq2)

Note: The first 5 steps are done in python and the scripts are saved in scripts/python/ folder.\
The 6th and last step is done in R and scripts are saved in scripts/R/ folder.

## 1. Downloading the raw data from GEO using SRA Toolkit
As mentioned above, we will be using the dataset from article: [ONECUT2 is a driver of neuroendocrine prostate cancer - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336817/#Sec11). The dataset is available in Gene Expression Omnibus(GEO) under the accession number GSE106305 available here: [GEO Accession viewer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305).

This dataset consists of **RNA-seq samples** from two prostate cancer cell lines:  **LNCaP** and **PC3**. Each cell line is cultured under **two oxygen conditions**: **Normoxia** and **Hypoxia**, with **two biological replicates** per condition.  

The goal of this analysis is to **identify genes that show differential expression under hypoxia** compared to normoxia within each cell line.


### 📊 Sample Metadata

| Sample Name | GSM Identifier | SRA Identifier (SRX) | SRA Runs (SRR, download these) |
|--------------|----------------|----------------------|--------------------------------|
| LNCaP_RNA-Seq_Empty_Vector_Normoxia_rep1 | GSM3145509 | SRX4096735 | SRR7179504, SRR7179505, SRR7179506, SRR7179507 |
| LNCaP_RNA-Seq_Empty_Vector_Normoxia_rep2 | GSM3145510 | SRX4096736 | SRR7179508, SRR7179509, SRR7179510, SRR7179511 |
| LNCaP_RNA-Seq_Empty_Vector_Hypoxia_rep1 | GSM3145513 | SRX4096739 | SRR7179520, SRR7179521, SRR7179522, SRR7179523 |
| LNCaP_RNA-Seq_Empty_Vector_Hypoxia_rep2 | GSM3145514 | SRX4096740 | SRR7179524, SRR7179525, SRR7179526, SRR7179527 |
| PC3_RNA-Seq_siCtrl_Normoxia_rep1 | GSM3145517 | SRX4096743 | SRR7179536 |
| PC3_RNA-Seq_siCtrl_Normoxia_rep2 | GSM3145518 | SRX4096744 | SRR7179537 |
| PC3_RNA-Seq_siCtrl_Hypoxia_rep1 | GSM3145521 | SRX4096747 | SRR7179540 |
| PC3_RNA-Seq_siCtrl_Hypoxia_rep2 | GSM3145522 | SRX4096748 | SRR7179541 |

### 🧾 Explanation of Identifiers

- **GSM Identifier (Gene Expression Omnibus Sample):**  
  A unique ID assigned to each biological sample in the **GEO** database.  
  It represents a specific experimental condition or replicate (e.g., LNCaP normoxia rep1).

- **SRA Identifier (SRX - Sequence Read Experiment):**  
  Corresponds to an **experiment record** in the **Sequence Read Archive (SRA)**,  
  describing how the sequencing was done for that particular sample.

- **SRA Runs (SRR - Sequence Read Runs):**  
  These are the **actual raw sequencing files** generated for an experiment.  
  Each SRX may have one or multiple SRR entries (individual sequencing runs).  
  These are the files you download using tools like `prefetch` or `fastq-dump`.

### Commands

Installing SRA-toolkit:
```
sudo apt update
sudo apt install sra-toolkit   
```
Sample download of one SRA file:\
SRA file SRR7179504 was downloaded using `prefetch`.\
This file was converted into fastq using `fastq-dump`.

```
prefetch SRR7179504   

fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR7179504.sra   

```
Automating the download of all 20 SRA files using python:

In order to automate the download and conversion of SRR files into FASTQ files, we run a python script which is added at scripts/python/fastq_download.py. \
Note: In the python program the sra_path mentioned is according to the my system. You will need to change it according to your system path to run in your systems.

To run the script we use

```
python fastq_download.py
```  
## 2. Quality control with FastQC and MultiQC

### FastQC

Now that we have gotten our sequencing data in fastq format, we need to perform qualtiy control on this sequencing data. For this we use the `fastqc` command:

```
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
```  
FASTQC performs quality control (QC) checks on raw sequencing reads. This generates a quality report in HTML format for each sample with summarized metrics: Per base sequence quality, per sequence quality score, per base sequence content, per sequence GC content, per base N content, Sequence length distribution, Sequence duplication level, Over represented sequences and adapter content.

This step is required to ensure that the sequencing data is reliable before any downstream analysis. Poor quality reads, adapter contamination and sequencing biases can distort alignment and gene quantification leading to false differential expression results. With FastQC we can detect these potential issues before alignment so that:
- **Poor-quality bases can be trimmed (using Trimmomatic).**
- **Contaminated or biased samples can be flagged or excluded.**

Factors mentioned in FASTQC and what they represent
#### **1. Per Base Sequence Quality**
**Represents:**  
The average Phred quality score at each base position across all reads.  
It shows how confident the sequencer was about calling each base.

**Potential Problems:**  
- Quality typically drops toward the 3′ end of reads.  
- A sharp decline below a Phred score of 20 (99% accuracy) suggests unreliable base calls.  

**Solution:**  
- Trim low-quality bases using tools like **Trimmomatic**, **Cutadapt**, or **fastp**.  
- Consider discarding reads below a certain average quality threshold (e.g., Q20 or Q30).

---

#### **2. Per Sequence Quality Scores**
**Represents:**  
Distribution of mean quality scores per read.  
A healthy dataset has a narrow peak toward the high-quality end.

**Potential Problems:**  
- A wider or bimodal distribution indicates that a subset of reads is low quality.  
- May result from instrument malfunction or poor cluster generation.

**Solution:**  
- Filter out low-quality reads before alignment.  
- Re-run sequencing if a large portion of reads is poor.

---

#### **3. Per Base Sequence Content**
**Represents:**  
Proportion of A, T, G, and C bases at each position.  
In random libraries, these proportions should be roughly equal.

**Potential Problems:**  
- Uneven base composition (especially at the start) suggests **bias from random priming**, **adapter remnants**, or **contamination**.  
- In RNA-seq, a small bias near the start of reads is often normal.

**Solution:**  
- Trim adapter or biased regions.  
- Ensure proper library preparation and random priming.

---

#### **4. Per Sequence GC Content**
**Represents:**  
GC distribution across all reads.  
It should approximate a normal distribution for most species.

**Potential Problems:**  
- A shifted or multimodal distribution suggests **contamination** (e.g., bacterial reads) or **bias in amplification**.  

**Solution:**  
- Remove contaminant sequences using **Kraken2**, **FastQ Screen**, or similar tools.  
- Check whether GC bias is biologically expected (e.g., transcriptome of GC-rich species).

---

#### **5. Per Base N Content**
**Represents:**  
Percentage of bases that were called as “N” (uncertain base) at each position.

**Potential Problems:**  
- High or position-specific N content indicates **low signal intensity**, **instrument errors**, or **overexposure**.

**Solution:**  
- Trim or filter reads containing many Ns.  
- If frequent, re-sequencing may be needed.

---

#### **6. Sequence Length Distribution**
**Represents:**  
Distribution of read lengths in the file.

**Potential Problems:**  
- Variable lengths indicate **trimming or adapter removal** has already been performed.  
- Unexpected short reads can bias alignment or indicate incomplete trimming.

**Solution:**  
- Ensure consistent read length after preprocessing.  
- Filter out excessively short reads (e.g., <30 bp for RNA-seq).

---

#### **7. Sequence Duplication Levels**
**Represents:**  
Fraction of duplicate reads in the dataset.  
High duplication can indicate overamplification or low library diversity.

**Potential Problems:**  
- High duplication (>50%) may arise from **PCR bias**, **low input material**, or **over-sequencing of few fragments**.

**Solution:**  
- Remove PCR duplicates after alignment (e.g., **Picard MarkDuplicates**).  
- Improve library preparation to increase complexity.

---

#### **8. Overrepresented Sequences**
**Represents:**  
Sequences that occur more frequently than expected.

**Potential Problems:**  
- Adapter contamination, primer dimers, rRNA reads, or other repetitive contamination.

**Solution:**  
- Identify and remove adapters or contaminants using **Cutadapt**, **Trim Galore**, or **BBduk**.  
- Use rRNA depletion kits in RNA-seq library prep if needed.

---

#### **9. Adapter Content**
**Represents:**  
Proportion of reads containing adapter sequences.

**Potential Problems:**  
- Indicates incomplete trimming or small fragment sizes relative to read length.  
- A high percentage leads to alignment errors and false expression counts.

**Solution:**  
- Trim adapters using **Trimmomatic**, **fastp**, or **Cutadapt** before alignment.  
- Review library prep fragment size distribution.

---

#### ✅ **Summary**

| Category | Common Cause | Impact | Fix |
|-----------|---------------|--------|-----|
| Low-quality bases | Signal decay | Miscalled bases | Trim or filter |
| Uneven GC/base content | Bias or contamination | Mapping errors | Clean and verify library |
| High duplication | Overamplification | Quantification bias | Remove duplicates |
| Adapter contamination | Short inserts | Alignment errors | Trim adapters |

### MultiQC

Often, we deal with a lot of fastq files and it is not possible to go through the fastqc report of each fastq file. MultiQC automatically compiles all these reports, allowing quick visualization of global quality trends, detection of systematic biases, and comparison of sample-level metrics. This provides an efficient overview of sequencing performance and helps identify any outlier or low-quality samples before downstream analysis.

We can get the MultiQC report using `multiqc` command:

```
multiqc fastqc_results/ -o multiqc_report/
```  
This combines all the fastqc results in fastqc_results folder to generate a multiqc_report file.

## 3. Trimming (Optional) with Trimmomatic
If we observe some issues with the sequencing data and need to remove certain parts of reads due to low quality, we use Trimmomatic tool. Using this we can remove low-quality bases and adapter contamination, improving read quality before alignment.

### 🔧 Common Parameters

| Parameter | Description | Example | Notes / When to Use |
|------------|--------------|----------|----------------------|
| **ILLUMINACLIP:** | Removes adapter or primer sequences using a provided adapter file. | `ILLUMINACLIP:TruSeq3-SE.fa:2:30:10` | First argument = adapter file; 2: max mismatches; 30: palindrome clip threshold; 10: simple clip threshold. |
| **LEADING:** | Trims low-quality bases from the start (5′ end) of each read. | `LEADING:20` | Removes poor-quality bases at the beginning. |
| **TRAILING:** | Trims low-quality bases from the end (3′ end) of each read. | `TRAILING:20` | Removes poor-quality trailing bases. |
| **SLIDINGWINDOW:** | Performs dynamic trimming once the average quality in a window falls below a threshold. | `SLIDINGWINDOW:4:20` | Commonly used — balances stringency and read retention. |
| **CROP:** | Keeps only the first *N* bases from each read. | `CROP:100` | Useful if you want to trim reads to a uniform length. |
| **HEADCROP:** | Removes the first *N* bases from each read. | `HEADCROP:10` | Good for removing biased sequence at the start (e.g., random hexamer bias). |
| **MINLEN:** | Discards reads shorter than a given length after trimming. | `MINLEN:36` | Ensures very short reads don’t cause mapping errors. |
| **AVGQUAL:** | Drops reads with an average quality score below a threshold. | `AVGQUAL:20` | Acts as a global quality filter for each read. |
| **TOPHRED33 / TOPHRED64** | Converts quality scores to specified encoding. | `TOPHRED33` | Useful when mixing datasets with different encodings (rare today). |

---

Example command:

```
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 fastq/LNCAP_Hypoxia_S1.fastq.gz fastq/LNCAP_Hypoxia_S1_trimmed.fastq TRAILING:10 -phred33 
``` 
Note: 1. This command removes bases from the 3' end of each read till it encounters a base which has quality score above 10. \
2. Trimmomatic runs on java so make sure java is installed in your system.\

After trimming, rerun FASTQC to confirm improvements in sequencng reads. If FASTQC results are satisfactory, proceed to alignment. If low-quality bases persist, perform additional trimming using more stringent parameters like higher quality threshold. If unavoidable, re-sequencing must be considered.

### Concatenating FASTQ Files
In our data we have some samples having multiple SRA runs. In order to simplify the data, we concatenate all these runs into a single file per sample using `cat`. Samples having single SRA runs are simply renamed using `mv`. 

```
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz > LNCAP_Normoxia_S1.fastq.gz
cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz > LNCAP_Normoxia_S2.fastq.gz
cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz > LNCAP_Hypoxia_S1.fastq.gz
cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz > LNCAP_Hypoxia_S2.fastq.gz
mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz
```

 It is also recommended to remove all the individual SRA-run files using `rm SRR*` command. After concatenation, renaming and removing the SRA files, only final 8 FASTQ files will remain, ready for alignment.


## 4. Alignment with HISAT2

### Setting up the reference genome

Whenever we need to do RNA seq analysis, we need to align our read data to a reference genome. Here we use the GRCh38 human genome and the corresponding annotation file (GTF file). 

#### Brief overview of a GTF file: \
A **GTF (Gene Transfer Format)** file provides **gene and transcript annotations** for a reference genome. 
It describes where genes, exons, CDS (coding sequences), and other genomic features are located, allowing tools like **HISAT2**, **StringTie**, and **featureCounts** to map reads to known genes accurately.

##### Inside a GTF File:

Each line in a GTF file corresponds to a **feature** on the genome (such as a gene, exon, or transcript).  
It contains **9 tab-separated columns**, as shown below:

| Column | Name | Description | Example |
|:-------:|------|--------------|----------|
| 1 | **seqname** | Chromosome or scaffold name | `chr1` |
| 2 | **source** | Program or database that generated the feature | `ENSEMBL` |
| 3 | **feature** | Type of genomic feature | `gene`, `exon`, `CDS`, `transcript` |
| 4 | **start** | Start coordinate of the feature | `11869` |
| 5 | **end** | End coordinate of the feature | `14409` |
| 6 | **score** | Confidence score (often `.` if not applicable) | `.` |
| 7 | **strand** | DNA strand (`+` or `-`) | `+` |
| 8 | **frame** | Reading frame for CDS features (`0`, `1`, or `2`) | `0` |
| 9 | **attributes** | Key-value pairs giving extra info like gene name, transcript ID, gene type | `gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_biotype "transcribed_unprocessed_pseudogene";` |

```
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz
wget ftp://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
gunzip Homo_sapiens.GRCh38.114.gtf.gz
```
### Read Alignment
Using the cleaned FASTQ files we run the script hisat2alignment.sh which is provided under scripts/bash/hisat2alignment.sh. This step produces BAM files along with an index .bai file foe each FASTQ file.

Brief overview of SAM and BAM files:
- SAM (Sequence Alignment/Map):
A plain-text format that stores information about how sequencing reads align to a reference genome. Each line represents one read and includes details such as read name, alignment position, mapping quality, and CIGAR string (which describes mismatches, insertions, or deletions).

- BAM (Binary Alignment/Map):
The binary, compressed version of a SAM file. It contains the same information but takes up less space and can be quickly indexed for random access. BAM files are used for most downstream analyses and visualization in genome browsers like IGV.

#### Sample Commands found inside the hisat2alignment.sh file.
```
# Specify the genome index and input FASTQ file

GENOME_INDEX="grch38/genome"
FASTQ_FILE="fastq/LNCAP_Hypoxia_S1.fastq.gz"

# Define output file name
OUTPUT_BAM="LNCAP_Hypoxia_S1.bam"

# Run HISAT2 alignment and pipe output directly to Samtools for sorting
hisat2 -q -x $GENOME_INDEX -U $FASTQ_FILE | samtools sort -o $OUTPUT_BAM

# Index the sorted BAM file for downstream analysis
samtools index $OUTPUT_BAM

echo "Alignment and indexing complete for $FASTQ_FILE"
```






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
