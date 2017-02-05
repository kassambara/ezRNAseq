
(In development, don't install)

ezRNAseq: Analyze easily RNAseq data
====================================

System requirements
-------------------

The following programs should be installed if you want to align RNAseq data using ezRNAseq:

-   [GNU Parallel](https://www.gnu.org/software/parallel/): GNU parallel is a shell tool for executing jobs in parallel using multiple cores in a computer or multiple computers.
-   [STAR](https://github.com/alexdobin/STAR) software for RNA sequencing data alignment
-   [SAMtools](http://samtools.sourceforge.net/): SAM Tools provide various utilities for manipulating alignments data in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.

Pandoc
------

required for the function ez\_report()

<https://github.com/rstudio/rmarkdown/blob/master/PANDOC.md>

    cp -r /usr/lib/rstudio-server/bin/pandoc/ ~/bin/pandoc/
    cp -r /usr/lib/rstudio-server/bin/pandoc/pandoc-citeproc ~/bin/pandoc/

Add this to your profile

``` bash
export PATH=$PATH:~/bin/pandoc
export PATH=$PATH:~/bin/pandoc-citeproc
```

Download reference genome
-------------------------

Time ~40 min

    library(ezRNAseq)
    install_igenome("hg19", result_dir = "~/genomes")

Create STAR index
-----------------

~ 45 min

    fasta_file = "~/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
    gtf_file = "~/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
    starIndex_dir = "~/genomes/Homo_sapiens/UCSC/hg19/Sequence/StarIndex"

    create_star_index(fasta_file, gtf_file, 
        thread = 25, starIndex_dir = starIndex_dir)

**fasta\_file** and **gtf\_file** are the path to fasta and gtf files respectively.

Mapping
-------

1.  Mapping to reference genome is done using STAR.
2.  Required program: **STAR** & **samtools**
3.  Required external data: reference genome

**Outputs**:

-   Name sorted BAM files
-   Raw count data
-   Normalized count data using DESeq2
