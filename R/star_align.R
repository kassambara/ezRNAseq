#' @include utilities.R
NULL
#'Read Aligning
#'@description
#'Read aligning using STAR software.
#'@param data_dir the directory containing the data. It should contains a subdirectory called "FASTQ" which
#'contains fastq files. Default value is the working directory
#'@param samples.annotation the path to "samples.txt" tabulation file containing sample annotations.
#'This file should contain at least 3 columns:
#'\itemize{
#'\item name: corresponding to sample name. It will be used to report the result
#'\item fastq1: the first fastq file name for paired-end sequencing
#'\item fastq2: the second fastq file name for paired-end sequencing.
#'This is optional if single-end sequencing is used
#'}
#'@param star.index reference genome for STAR.
#'@param result.dir The absolute path of result directory.
#'@param keep a character vector specifying the file to be kept.
#' Allowed values include c("name_sorted_bam", "chr_sorted_bam", "sam" ).
#'@param thread Number of threads to be used. This depends to the available computer ressources.
#'@param pairedEnd specify if the data are paired-end sequencing data. Default is TRUE.
#'@param fastq.gz specify if the FASTQ files are compressed. Default is TRUE.
#'@details
#' The function star_align(), requires STAR and samtools programs to work.
#' Make sure that they are installed. \cr\cr
#' The workflow is as follow:\cr
#' 1. Align all FASTQ files using STAR. BAM files are generated.\cr
#' 2. Organize BAM files (sorting and indexing). Samtools program required.\cr
#' 3. Count read using Bioconductor packages: "GenomicFeatures", "Rsamtools", "GenomicAlignments" and "BiocParallel" required.\cr
#'
#' @return Three subdirectories are created: \cr
#' \itemize{
#' \item BAM: containing unsorted and sorted BAM files, and BAM index file.
#' BAM comtent is sorted by read names (*_name_sorted.bam) or by chromosome (*_sorted.bam).
#' The index files are of form *_sorted.bam.bai. *_name_sorted.bam files are used for read counting
#' using hseq-count or R. *_sorted.bam and *_sorted.bam.bai files can be used for IGV visualization.
#' \item COUNT: containing read counting results. It contains the following files:\cr
#' 1. se.RDATA: containing an object "se" which is an object of class SummarizedExperiment
#' }
#' @name star_align
#' @rdname star_align
#'@export
star_align <- function(data_dir = getwd(), samples.annotation = "samples.txt",
                       star.index = "/eqmoreaux/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/StarIndex/",
                       result.dir = getwd(),
                       keep = c("name_sorted_bam"),
                       pairedEnd = TRUE, fastq.gz = TRUE, thread = 10)
{

  # Read samples.txt file
  samples <- read.delim(file=samples.annotation, header = TRUE, row.names = NULL)

  # Go to the directory containing the FASTQ files
  oldwd <- getwd()
  setwd(file.path(data_dir, "FASTQ"))

  # Create result dirs
  # ++++++++++++++++++++++++++++++++++++++++
  # create_dir(file.path(result.dir, "SAM"))
  create_dir(file.path( result.dir, "BAM"))

  # STAR Alignement
  # ++++++++++++++++++++++++++++++++++++++++
  # Result = BAM file
  message("\n- Starting Alignment...\n")

  # Option to read compressed FASTQ files
  if(fastq.gz) read_file_command <- "--readFilesCommand gunzip -c"
  else read_file_command <- ""

  if(pairedEnd)
    cmd <- with(samples, paste("STAR --genomeDir", star.index,
                               read_file_command,
                               "--readFilesIn", fastq1, fastq2,
                               "--outSAMstrandField intronMotif",
                               "--outFileNamePrefix", file.path(result.dir, "BAM", paste0(name, "_")),
                               "--outSAMtype BAM Unsorted",
                               "--runThreadN", thread, sep=" " ))
  else cmd <- with(samples, paste("STAR --genomeDir", star.index,
                                  read_file_command,
                                  "--readFilesIn", fastq1,
                                  "--outSAMstrandField intronMotif",
                                  "--outFileNamePrefix", file.path(result.dir, "BAM", paste0(name, "_")),
                                  "--outSAMtype BAM Unsorted",
                                  "--runThreadN", thread, sep=" " ))

  for(c in cmd) {
    cat(c, "\n--------------------------------------\n")
    system(c) # Run alignment
  }


  # Renaming BAM file
  # ++++++++++++++++++++++++++++++++++++++
  message("- Renaming BAM Files....\n")
  initial_name <- file.path(result.dir, "BAM", paste0(samples$name, "_Aligned.out.bam"))
  final_name <- file.path(result.dir, "BAM", paste0(samples$name, ".bam"))
  for(i in 1:length(samples$name)) file.rename(initial_name[i], final_name[i])

  # Convert SAM -> BAM file
  # ++++++++++++++++++++++++++++++++++++++++
  #  ls *.sam | parallel "samtools view -bS {} > ../BAM/{.}.bam"
  # message("- Converting SAM to BAM Files...\n")
  # sam_dir <- file.path(result.dir, "SAM")
  # bam_dir <- file.path(result.dir, "BAM/{.}.bam")
  # setwd(sam_dir) # Go to SAM dir
  # cmd <- paste( 'ls',  '*.sam | parallel -t -j', thread,
  #               '"samtools view -bS {} >', bam_dir, '"', sep = " ")
  # system(cmd)

  # Remove or keep SAM file
  # +++++++++++++++++++++++++++++++
  # if(!("sam" %in% keep)){
  #   unlink(file.path(result.dir, "SAM"),
  #          recursive = TRUE, force = TRUE)
  # }


  # Organizing BAM files
  #++++++++++++++++++++++++++++++++++++++++++
  # - sort by name for htseq-count
  # - sort by chromosome and create index for IGV
  #  ls *.bam | parallel "samtools sort file.bam  file_name_sorted.bam"
  message("- Organizing BAM files... \n")
  setwd(file.path(result.dir, "BAM"))

  if("name_sorted_bam" %in% keep) {
    message("- Sorting reads by name... \n")
    create_dir("name_sorted")
    cmd <- paste('ls', '*.bam | parallel -t -j', thread,
                 '"samtools sort -n {} name_sorted/{.}_name_sorted"', sep = " ")
    system(cmd)
  }

  if("chr_sorted_bam" %in% keep){
    message("- Sorting reads by chromosome... \n")
    create_dir("chr_sorted")
    system(paste('ls', '*.bam | parallel -t -j', thread,
                 '"samtools sort {} chr_sorted/{.}_sorted"', sep = " ")) # by chromosome
    setwd("chr_sorted")
    system(paste('ls', '*.bam | parallel -t -j', thread,
                 '"samtools index {} {}"', sep = " ")) #index
  }

  # Removing unsorted BAM
  # +++++++++++++++++++++++++++++
  system("rm *.bam")

  setwd(oldwd)
}
