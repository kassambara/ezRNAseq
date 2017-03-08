#' @include utilities.R
NULL
#'Read Aligning
#'@description Read aligning using Bowtie2 software.
#'@param data_dir the directory containing the data. It should contains a
#'  subdirectory called "FASTQ" which contains fastq files. Default value is the
#'  working directory
#'@param samples.annotation the path to "samples.txt" tabulation file containing
#'  sample annotations. This file should contain at least 3 columns: \itemize{
#'  \item name: corresponding to sample name. It will be used to report the
#'  result \item fastq1: the first fastq file name for paired-end sequencing
#'  \item fastq2: the second fastq file name for paired-end sequencing. This is
#'  optional if single-end sequencing is used }
#'@param bowtie2.index reference genome for bowtie2.
#'@param genome.fa path to the fasta file of reference genome. Make sure that
#'  samtools index has been created for the reference genome.
#'@param result.dir The absolute path of result directory.
#'@param keep a character vector specifying the file to be kept. Allowed values
#'  include c("name_sorted_bam", "chr_sorted_bam", "sam" ).
#'@param thread Number of threads to be used. This depends to the available
#'  computer ressources.
#'@param memory a numeric value specifying the memory limit to be used. Default
#'  is 3G.
#'@param rmdup logical. If TRUE, remove PCR duplicates.
#'@param call_variant logical. If TRUE, call variant.
#'@details The function bowtie2_align(), requires bowtie2 and samtools programs
#'to work. Make sure that they are installed. \cr\cr The workflow is as
#'follow:\cr 1. Align all FASTQ files using bowtie2. BAM files are generated.\cr
#'2. Organize BAM files (sorting and indexing). Samtools program required.\cr 3.
#'Count read using Bioconductor packages: "GenomicFeatures", "Rsamtools",
#'"GenomicAlignments" and "BiocParallel" required.\cr
#'
#'@return Three subdirectories are created: \cr \itemize{ \item BAM: containing
#'  unsorted and sorted BAM files, and BAM index file. BAM comtent is sorted by
#'  read names (*_name_sorted.bam) or by chromosome (*_sorted.bam). The index
#'  files are of form *_sorted.bam.bai. *_name_sorted.bam files are used for
#'  read counting using hseq-count or R. *_sorted.bam and *_sorted.bam.bai files
#'  can be used for IGV visualization. \item COUNT: containing read counting
#'  results. It contains the following files:\cr 1. se.RDATA: containing an
#'  object "se" which is an object of class SummarizedExperiment }
#'@name bowtie2_align
#'@rdname bowtie2_align
#'@export
bowtie2_align <- function(data_dir = getwd(), samples.annotation = "samples.txt",
                       bowtie2.index = "/eqmoreaux/genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome",
                       genome.fa = "/eqmoreaux/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa",
                       result.dir = getwd(),
                       keep = c("bam"), thread = 20,
                       memory = 3000000000, rmdup = FALSE, call_variant = FALSE)
{

  # Read samples.txt file
  samples <- read.delim(file=samples.annotation, header = TRUE, row.names = NULL)

  pairedEnd <- !is.null(samples$fastq2)
  file.ext <- tail(unlist(strsplit(as.vector(samples$fastq1)[1], ".", fixed = TRUE)), 1)
  fastq.gz <- file.ext == "gz"

  # Create result dirs
  # ++++++++++++++++++++++++++++++++++++++++
  create_dir(file.path(result.dir, "SAM"))
  create_dir(file.path( result.dir, "BAM"))
  create_dir(file.path( result.dir, "LOG"))


  # Go to the directory containing the FASTQ files
  oldwd <- getwd()
  setwd(file.path(data_dir, "FASTQ"))

  # bowtie2 Alignement
  # ++++++++++++++++++++++++++++++++++++++++
  # Result = BAM file
  message("\n- Starting Alignment...\n")

  # Option to read compressed FASTQ files
  if(fastq.gz) read_file_command <- "--readFilesCommand gunzip -c"
  else read_file_command <- ""

  if(pairedEnd)
    cmd <- with(samples, paste("(bowtie2 -p", thread, "-x ", bowtie2.index,
                               "-1 ", fastq1, "-2 ", fastq2,
                               "-S", file.path(result.dir, "SAM", paste0(name, ".sam")),
                               ") 2>", file.path(result.dir, "LOG", paste0(name, ".log")),
                               sep=" " ))
  else cmd <- with(samples, paste("(bowtie2 -p", thread, "-x ", bowtie2.index,
                                  fastq1,
                                  "-S", file.path(result.dir, "SAM", paste0(name, ".sam")),
                                  ") 2>", file.path(result.dir, "LOG", paste0(name, ".log")),
                                  sep=" " ))

  for(c in cmd) {
    cat(c, "\n--------------------------------------\n")
    system(c) # Run alignment
  }

  # Write alignment params
  # +++++++++++++++++++++
  sink(file.path(result.dir, "LOG", "_alignment.params"))
  for(c in cmd) cat(c, "\n\n")
  sink()

  # Convert SAM -> BAM file
  # ++++++++++++++++++++++++++++++++++++++++
  #  ls *.sam | parallel "samtools view -bS {} > ../BAM/{.}.bam"
  # message("- Converting SAM to BAM Files...\n")
   sam_dir <- file.path(result.dir, "SAM")
   bam_dir <- file.path(result.dir, "BAM/{.}.bam")
   setwd(sam_dir) # Go to SAM dir
   cmd <- paste( 'ls',  '*.sam | parallel -t',
                  '"samtools view -bS {} >', bam_dir, '"', sep = " ")
   system(cmd)

  # Remove or keep SAM file
  # +++++++++++++++++++++++++++++++
   if(!("sam" %in% keep)){
    unlink(file.path(result.dir, "SAM"),
            recursive = TRUE, force = TRUE)
   }

   # Organizing BAM files
   #++++++++++++++++++++++++++++++++++++++++++
   # - sort by chromosome and create index for IGV
   #  ls *.bam | parallel "samtools sort file.bam  file_name_sorted.bam"
   message("- Organizing BAM files... \n")
   setwd(file.path(result.dir, "BAM"))
   message("- Sorting reads by chromosome... \n")
   create_dir("chr_sorted")
   system(paste('ls', '*.bam | parallel -t',
                '"samtools sort {} chr_sorted/{.}_chrsorted -m', memory, '"', sep = " ")) # by chromosome
   setwd("chr_sorted")
   if(!rmdup) system(paste('ls', '*.bam | parallel -t',
                '"samtools index {} {}"', sep = " ")) #index
   system("rm ../*.bam") # Removing unsorted BAM
   system("mv *.bam* ../") # moving chr_sorted content to BAM
   system("rm -r chr_sorted")

   # Removing PCR duplicates
   #++++++++++++++++++++++++++++++++++++++
   if(rmdup){

     setwd(file.path(result.dir, "BAM"))
     create_dir(file.path(result.dir, "RMDUP"))
     system('ls *.bam | parallel -t "samtools rmdup {} ../RMDUP/{.}_rmdup.bam"')

     setwd(file.path(result.dir, "RMDUP"))
     system(paste('ls', '*.bam | parallel -t',
                  '"samtools index {} {}"', sep = " "))

     setwd(result.dir)
     system('rm -r BAM')
     system('mv RMDUP BAM')

   }

   # Calling variants
   #++++++++++++++++++++++++++++++++++++++
   if(call_variant){
     create_dir(file.path(result.dir, "VCF"))
     setwd(file.path(result.dir, "BAM"))

     cmd <- paste( 'ls *.bam | parallel -t',
                   'samtools mpileup -Buf', genome.fa, '{} |',
                   'bcftools view - -vcg > ../VCF/{.}.vcf',
                    sep = "")
     system(cmd)
   }

  setwd(oldwd)
}
