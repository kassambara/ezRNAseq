ezRNAseq: Analyze easily RNAseq data
==================
    
R Functions
-----------
    
    
(In development, don't install)

Installation
-------------
     
## Pandoc
    
    
required for the function ez_report()
    
https://github.com/rstudio/rmarkdown/blob/master/PANDOC.md
   
   
   
```
cp -r /usr/lib/rstudio-server/bin/pandoc/ ~/bin/pandoc/
cp -r /usr/lib/rstudio-server/bin/pandoc/pandoc-citeproc ~/bin/pandoc/
```
    
    

   
Add this to your profile
     
     
```bash
export PATH=$PATH:~/bin/pandoc
export PATH=$PATH:~/bin/pandoc-citeproc
```
     
     
Download reference genome
-------------------
   
Time ~40 min
  
```
library(ezRNAseq)
install_igenome("hg19", result_dir = "~/genomes")
```
   
   
Create STAR index
-----------------
  
~ 45 min
    
```
fasta_file = "~/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
gtf_file = "~/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
starIndex_dir = "~/genomes/Homo_sapiens/UCSC/hg19/Sequence/StarIndex"

create_star_index(fasta_file, gtf_file, 
    thread = 25, starIndex_dir = starIndex_dir)
```
   
**fasta_file** and **gtf_file** are the path to fasta and gtf files respectively.
  
    
Mapping
--------------
     
     
1. Mapping to reference genome is done using STAR.
2. Required program: **STAR**
3. Required external data: reference genome
   
   
**Outputs**: SAM files

    
    
Organizing SAM/BAM files
-----------------------
   
   
1. Conversion of SAM -> BAM files
2. Sort Reads in BAM files by name and chromosome
3. Create index file for chromosome sorted Bam files
   
   
Required program: **samtools**

   
**Outputs**: /BAM  

    - (*.) bam files
    - (*.) _sorted.bam files: reads sorted by chromosomes
    - (*.)_sorted.bam.bai files: index files
    - (*.) _name_sorted.bam files: reads sorted by name
       
   
Read counting using Bioconductor package
---------------
    
    
1. Required external data: gtf file
2. Required Bioconductor packages: 
    - GenomicFeatures
    - Rsamtools
    - GenomicAlignments
    - BiocParallel
    
    
**Output**:  
   
- SummarizedExperiment
    
    
Functions
-------------------
    
    
- star_align(): Align reads with STAR
- count_reads(): Counting reads using R/Bioconductor
- normalize_counts(): Normalizing count with DESeq2. Creates count.norm.txt and rlog.data.txt.
- check_count_data(): Data quality assessment; 
   
rnaseq_workflow(): Runs the RNAseq workflow: align, count, normalize and check.
    
    

